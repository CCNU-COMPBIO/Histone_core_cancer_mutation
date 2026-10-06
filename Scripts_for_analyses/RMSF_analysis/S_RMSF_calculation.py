#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSF-based Pseudo-symmetry Deviation Score (S_RMSF) - Global (octamer) level
Strategy: compute per dyad chain pair, then pool and summarize by residue count.
Output: mean + SEM + per-run S_RMSF.
Alignment: use CHAIN_START to map residues to the canonical histone numbering, align by canonical numbering.
Trimming: symmetrically remove N_TRIM residues from both ends on the shared canonical numbering.

One run produces two variants (the R2-2 excluded/no-exclusion comparison):
  Variant A "excluded": all systems x run1-run3, reads regular files,
                 outputs SRMSF_PseudoSymmetry_Deviation_Trim{N}_PerRunDetail.csv
  Variant B "no-exclusion": only processes systems that have _noexcl files (others are skipped, no duplicate computation),
                 reads _noexcl files, and recomputes WT using the same run subset,
                 outputs ..._PerRunDetail_noexcl.csv
  Also outputs SRMSF_octamer_compare_excl_vs_noexcl.csv for filling in Table S17.

★ Folder names are already unified as <histone>_<mutation>_rw, e.g. H2A_R29Q_rw, H2B_F70L_rw, H3_E50K_rw.
  Therefore the MUTANTS keys use that naming directly (with the "_rw" suffix removed).
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
WT_SYSTEM = "WT_rw"
NOEXCL = "_noexcl"
SEQUENCE_CSV = "sequence.csv"

# Two variants
VARIANTS = [
    {"suffix": "",     "label": "excluded segments"},
    {"suffix": NOEXCL, "label": "no excluded segments"},
]

N_TRIM = 5

CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

DYAD_PAIRS = [
    {"mut": {"chain": "A", "histone": "H3"},  "non": {"chain": "E", "histone": "H3"}},
    {"mut": {"chain": "B", "histone": "H4"},  "non": {"chain": "F", "histone": "H4"}},
    {"mut": {"chain": "C", "histone": "H2A"}, "non": {"chain": "G", "histone": "H2A"}},
    {"mut": {"chain": "D", "histone": "H2B"}, "non": {"chain": "H", "histone": "H2B"}},
]

# ★ Keys must match the folder names (without "_rw") ★
#   Folder lookup: f"{mutant_key}_rw" -> H2A_R29Q_rw / H2B_F70L_rw / H3_E50K_rw ...
MUTANTS = {
    "H2A_R29Q": {"histone_mut": "H2A", "mut_chain": "C", "resid": 29},
    "H2A_E56K": {"histone_mut": "H2A", "mut_chain": "C", "resid": 56},
    "H2A_K74N": {"histone_mut": "H2A", "mut_chain": "C", "resid": 74},
    "H2B_F70L": {"histone_mut": "H2B", "mut_chain": "D", "resid": 70},
    "H2B_E105K": {"histone_mut": "H2B", "mut_chain": "D", "resid": 105},
    "H3_E50K":  {"histone_mut": "H3",  "mut_chain": "A", "resid": 50},
    "H3_E73K":  {"histone_mut": "H3",  "mut_chain": "A", "resid": 73},
    "H3_E105K": {"histone_mut": "H3",  "mut_chain": "A", "resid": 105},
}


# ========== Helper functions ==========
def read_sequence_csv(path: Path) -> pd.DataFrame:
    for enc in ["utf-8", "utf-8-sig", "gbk", "cp936", "latin1"]:
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError:
            continue
    else:
        raise UnicodeDecodeError(f"Failed to decode {path}")

    rename = {}
    if "#Res" in df.columns:
        rename["#Res"] = "Res"
    if "chain" in df.columns:
        rename["chain"] = "Chain"
    df = df.rename(columns=rename)
    df["Res"] = pd.to_numeric(df["Res"], errors="coerce")
    df = df.dropna(subset=["Res"]).copy()
    df["Res"] = df["Res"].astype(int)
    df["Chain"] = df["Chain"].astype(str).str.strip()
    return df


def build_res_to_std(seq_df: pd.DataFrame, chain_id: str, chain_start: int) -> dict:
    sub = seq_df[seq_df["Chain"] == chain_id]
    if sub.empty:
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv")
    first_res = sub["Res"].min()
    return {res: chain_start + (res - first_res) for res in sub["Res"]}


def scan_noexcl() -> dict:
    """Scan which systems/runs have _noexcl files (probe using H3 chain A)."""
    out = {}
    for mutant_key in MUTANTS:
        sys_name = f"{mutant_key}_rw"
        runs = [r for r in RUNS
                if (BASE_DIR / sys_name / r / f"rmsf-H3_chain_A{NOEXCL}.dat").exists()]
        if runs:
            out[sys_name] = runs
    return out


def load_single_run_rmsf(system: str, chain: str, histone: str,
                         res_to_std: dict, run: str, suffix: str = "") -> dict:
    """Read a single chain / single run RMSF, returns {canonical std: RMSF}; no trimming."""
    f = BASE_DIR / system / run / f"rmsf-{histone}_chain_{chain}{suffix}.dat"
    if not f.exists():
        return {}

    try:
        df = pd.read_csv(f, comment="#", sep=r'\s+', header=None,
                         engine='python', names=["Res", "RMSF"])
    except Exception:
        return {}

    df["Res"] = df["Res"].round().astype(int)
    result = {}
    for _, row in df.iterrows():
        std = res_to_std.get(row["Res"])
        if std is not None and not np.isnan(row["RMSF"]):
            result[std] = float(row["RMSF"])
    return result


def calc_s_rmsf_single_run(system: str, seq_df: pd.DataFrame, run: str,
                           n_trim: int = 5, suffix: str = "") -> float:
    """Single-run octamer S_RMSF: pool all four dyad pairs (accumulate numerator and denominator separately)."""
    total_abs_diff_sum = 0.0
    total_n = 0

    for pair in DYAD_PAIRS:
        mut_info = pair["mut"]
        non_info = pair["non"]

        try:
            mut_res_to_std = build_res_to_std(
                seq_df, mut_info["chain"], CHAIN_START[mut_info["chain"]])
            non_res_to_std = build_res_to_std(
                seq_df, non_info["chain"], CHAIN_START[non_info["chain"]])
        except ValueError:
            continue

        rmsf_mut = load_single_run_rmsf(system, mut_info["chain"], mut_info["histone"],
                                        mut_res_to_std, run, suffix)
        rmsf_non = load_single_run_rmsf(system, non_info["chain"], non_info["histone"],
                                        non_res_to_std, run, suffix)

        common_std = sorted(set(rmsf_mut.keys()) & set(rmsf_non.keys()))
        if len(common_std) > 2 * n_trim:
            common_std = common_std[n_trim:-n_trim]
        else:
            print(f"      [WARN] {mut_info['chain']}-{non_info['chain']} only "
                  f"{len(common_std)} common std positions, cannot trim, skip")
            continue

        for std in common_std:
            total_abs_diff_sum += abs(rmsf_mut[std] - rmsf_non[std])
        total_n += len(common_std)

    return total_abs_diff_sum / total_n if total_n > 0 else np.nan


def calc_s_rmsf_aggregate(system: str, seq_df: pd.DataFrame, runs,
                          n_trim: int = 5, suffix: str = ""):
    """Aggregate across runs (mean + SEM); returns (mean, sem, per_run_dict)."""
    per_run = {run: calc_s_rmsf_single_run(system, seq_df, run, n_trim=n_trim, suffix=suffix)
               for run in runs}
    valid_vals = [v for v in per_run.values() if not np.isnan(v)]
    if not valid_vals:
        return np.nan, np.nan, per_run
    mean_val = float(np.mean(valid_vals))
    sem_val = (float(np.std(valid_vals, ddof=1) / np.sqrt(len(valid_vals)))
               if len(valid_vals) > 1 else 0.0)
    return mean_val, sem_val, per_run


# ========== Single-variant computation ==========
def run_variant(variant: dict, noexcl_map: dict, wt_runs: list):
    suffix, label = variant["suffix"], variant["label"]

    print(f"\n{'#'*80}")
    print(f"#  VARIANT: {label}   (suffix = '{suffix or 'none'}')")
    if suffix:
        print(f"#  systems analysed: {sorted(noexcl_map.keys())}")
    print(f"#  WT runs: {wt_runs}")
    print(f"#  Fixed truncation: {N_TRIM} residues from each chain end")
    print(f"{'#'*80}")

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)

    # --- WT ---
    s_wt, sem_wt, wt_per_run = calc_s_rmsf_aggregate(
        WT_SYSTEM, seq_df, wt_runs, n_trim=N_TRIM, suffix=suffix)
    print(f">>> WT  S_RMSF = {s_wt:.4f} ± {sem_wt:.4f} Å   per-run: "
          + " | ".join(f"{r}: {v:.4f}" if not np.isnan(v) else f"{r}: N/A"
                       for r, v in wt_per_run.items()))

    all_results = []
    for mutant_key, info in MUTANTS.items():
        sys_name = f"{mutant_key}_rw"

        # ★ The no-exclusion variant only processes systems with _noexcl; others are skipped (no duplicate computation) ★
        if suffix and sys_name not in noexcl_map:
            print(f">>> {mutant_key}: no _noexcl files -> skipped in this variant")
            continue

        use_runs = noexcl_map.get(sys_name, RUNS) if suffix else list(RUNS)
        resid, h_mut, c_mut = info["resid"], info["histone_mut"], info["mut_chain"]

        s_mut, sem_mut, mut_per_run = calc_s_rmsf_aggregate(
            sys_name, seq_df, use_runs, n_trim=N_TRIM, suffix=suffix)

        delta_s = s_mut - s_wt if not (np.isnan(s_mut) or np.isnan(s_wt)) else np.nan
        delta_sem = (np.sqrt(sem_mut ** 2 + sem_wt ** 2)
                     if not (np.isnan(sem_mut) or np.isnan(sem_wt)) else np.nan)

        print(f">>> {mutant_key} ({h_mut} {c_mut}{resid})  runs={use_runs}")
        print(f"    Mutant S_RMSF = {s_mut:.4f} ± {sem_mut:.4f} Å   per-run: "
              + " | ".join(f"{r}: {v:.4f}" if not np.isnan(v) else f"{r}: N/A"
                           for r, v in mut_per_run.items()))
        if not np.isnan(delta_s):
            print(f"    ΔS_RMSF = {delta_s:+.4f} ± {delta_sem:.4f} Å")

        record = {
            "Variant": mutant_key,
            "Histone_Subunit": h_mut,
            "Mutated_Chain": c_mut,
            "Residue_Number": resid,
            "SRMSF_WT_Mean_A": s_wt,
            "SRMSF_WT_SEM_A": sem_wt,
            "SRMSF_Variant_Mean_A": s_mut,
            "SRMSF_Variant_SEM_A": sem_mut,
            "Delta_SRMSF_Mean_A": delta_s,
            "Delta_SRMSF_SEM_A": delta_sem,
            "Runs_Used": ",".join(use_runs),
        }
        for i, run in enumerate(RUNS, start=1):
            wt_v = wt_per_run.get(run, np.nan)
            mut_v = mut_per_run.get(run, np.nan)
            record[f"SRMSF_WT_Run{i}_A"] = wt_v
            record[f"SRMSF_Variant_Run{i}_A"] = mut_v
            record[f"Delta_SRMSF_Run{i}_A"] = (mut_v - wt_v
                                               if not (np.isnan(wt_v) or np.isnan(mut_v))
                                               else np.nan)

        all_results.append(record)

    # ---- Save ----
    out_csv = BASE_DIR / f"SRMSF_PseudoSymmetry_Deviation_Trim{N_TRIM}_PerRunDetail{suffix}.csv"
    df_out = pd.DataFrame(all_results)

    if len(df_out):
        df_out.to_csv(out_csv, index=False, float_format="%.6f")
        print(f"\n[OK] Saved to: {out_csv.resolve()}")
        summary_cols = ["Variant", "Histone_Subunit",
                        "SRMSF_WT_Mean_A", "SRMSF_Variant_Mean_A",
                        "Delta_SRMSF_Mean_A", "Delta_SRMSF_SEM_A"]
        print(df_out[summary_cols].round(4).to_string(index=False))
    else:
        print(f"\n[WARN] no records for variant '{label}' -> nothing written")
        df_out.to_csv(out_csv, index=False)

    return s_wt, sem_wt, df_out


# ========== Main ==========
def main():
    noexcl_map = scan_noexcl()
    wt_runs_b = sorted({r for rs in noexcl_map.values() for r in rs})

    print("=" * 80)
    print("[Variant B scan] Which systems/runs have _noexcl files (only these use Variant B)")
    if noexcl_map:
        for sys_name, rs in noexcl_map.items():
            print(f"  {sys_name:12s} -> {rs}")
        print(f"  WT runs used in Variant B = {wt_runs_b}")
    else:
        print("  No _noexcl files found; Variant B will be empty")
    print("=" * 80)

    # Variant A: excluded
    wt_a, _, df_a = run_variant({"suffix": "", "label": "excluded segments"},
                                noexcl_map, list(RUNS))

    # Variant B: no-exclusion (only systems with _noexcl; WT uses the same run subset)
    if noexcl_map:
        wt_b, _, df_b = run_variant({"suffix": NOEXCL, "label": "no excluded segments"},
                                    noexcl_map, wt_runs_b)
    else:
        wt_b, df_b = np.nan, pd.DataFrame()

    # ---- Comparison of the two variants (for filling in Table S17) ----
    if len(df_a) and len(df_b):
        cols = ["Variant", "SRMSF_Variant_Mean_A", "Delta_SRMSF_Mean_A"]
        m = df_a[cols].merge(df_b[cols], on="Variant", suffixes=("_excl", "_noexcl"))
        m["dSRMSF_diff"] = m["Delta_SRMSF_Mean_A_noexcl"] - m["Delta_SRMSF_Mean_A_excl"]
        cmp_csv = BASE_DIR / "SRMSF_octamer_compare_excl_vs_noexcl.csv"
        m.to_csv(cmp_csv, index=False, float_format="%.6f")

        print("\n" + "=" * 96)
        print("[Comparison] Octamer S_RMSF: excluded vs no-exclusion")
        print(f"  WT S_RMSF  excl = {wt_a:.4f}   noexcl = {wt_b:.4f}")
        print("=" * 96)
        print(m.round(4).to_string(index=False))
        print(f"\n[OK] {cmp_csv.name}")


if __name__ == "__main__":
    main()