#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSF-based Pseudo-symmetry Deviation Score (S_RMSF) - Full Per-Run Output
Strategy: Calculate separately for each dyad chain pair, then aggregate.
Output: Mean + SEM + per-run S_RMSF.
Alignment: Use CHAIN_START to map residues to canonical histone numbering, then align by canonical numbering.
Trimming: After sorting each chain by canonical numbering, remove a fixed N_TRIM residues from both ends.
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

N_TRIM = 5

# Canonical histone numbering corresponding to the first residue of each chain
CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

DYAD_PAIRS = [
    {"mut": {"chain": "A", "histone": "H3"}, "non": {"chain": "E", "histone": "H3"}},
    {"mut": {"chain": "B", "histone": "H4"}, "non": {"chain": "F", "histone": "H4"}},
    {"mut": {"chain": "C", "histone": "H2A"}, "non": {"chain": "G", "histone": "H2A"}},
    {"mut": {"chain": "D", "histone": "H2B"}, "non": {"chain": "H", "histone": "H2B"}},
]

MUTANTS = {
    "R29Q": {"histone_mut": "H2A", "mut_chain": "C", "resid": 29},
    "E56K": {"histone_mut": "H2A", "mut_chain": "C", "resid": 56},
    "K74N": {"histone_mut": "H2A", "mut_chain": "C", "resid": 74},
    "F70L": {"histone_mut": "H2B", "mut_chain": "D", "resid": 70},
    "E50K": {"histone_mut": "H3", "mut_chain": "A", "resid": 50},
    "E73K": {"histone_mut": "H3", "mut_chain": "A", "resid": 73},
    "E105K": {"histone_mut": "H3", "mut_chain": "A", "resid": 105},
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
    """
    Build the {original Res: canonical std} mapping for this chain.
    Canonical numbering = chain_start + (Res - first Res of this chain).
    """
    sub = seq_df[seq_df["Chain"] == chain_id]
    if sub.empty:
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv")
    first_res = sub["Res"].min()
    res_to_std = {}
    for res in sub["Res"]:
        std = chain_start + (res - first_res)
        res_to_std[res] = std
    return res_to_std


def trim_fixed_ends_std(std_rmsf: dict, n_trim: int = 5) -> dict:
    """Fixed removal of n_trim residues from each end by canonical numbering."""
    if not std_rmsf:
        return {}
    stds = sorted(std_rmsf.keys())
    if len(stds) <= 2 * n_trim:
        print(f"      [WARN] Chain too short ({len(stds)} residues) "
              f"for trimming {n_trim} from each end, returning empty")
        return {}
    trimmed_stds = stds[n_trim:-n_trim]
    return {s: std_rmsf[s] for s in trimmed_stds}


def load_single_run_rmsf(system: str, chain: str, histone: str,
                         res_to_std: dict, run: str) -> dict:
    """
    Load RMSF for a single chain and single run; return {canonical std: RMSF}.
    Note: do not trim ends here; trimming is deferred until common_std is determined uniformly.
    """
    f = BASE_DIR / system / run / RMSF_FILE_TEMPLATE.format(
        histone=histone, chain=chain)
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
    return result                      # Do not call trim_fixed_ends_std here


def calc_s_rmsf_single_run(system: str, seq_df: pd.DataFrame,
                           run: str, n_trim: int = 5) -> float:
    """
    Compute the full S_RMSF independently for a single run.
    After alignment, symmetrically trim both ends only on shared canonical positions.
    """
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

        rmsf_mut = load_single_run_rmsf(
            system, mut_info["chain"], mut_info["histone"],
            mut_res_to_std, run)
        rmsf_non = load_single_run_rmsf(
            system, non_info["chain"], non_info["histone"],
            non_res_to_std, run)

        common_std = sorted(set(rmsf_mut.keys()) & set(rmsf_non.keys()))

        # Symmetrically truncate at shared homologous positions
        if len(common_std) > 2 * n_trim:
            common_std = common_std[n_trim:-n_trim]
        else:
            # If there are too few shared positions to trim n_trim from each end, skip this chain pair
            print(f"      [WARN] {mut_info['chain']}-{non_info['chain']} "
                  f"only {len(common_std)} common std positions, "
                  f"cannot trim {n_trim} from each end, skip this pair")
            continue

        for std in common_std:
            total_abs_diff_sum += abs(rmsf_mut[std] - rmsf_non[std])
        total_n += len(common_std)

    return total_abs_diff_sum / total_n if total_n > 0 else np.nan


def calc_s_rmsf_aggregate(system: str, seq_df: pd.DataFrame, n_trim: int = 5):
    """
    Aggregate across runs (mean + SEM), and also return the independent S_RMSF for each run.
    Returns: (mean, sem, per_run_dict)
    """
    per_run = {}
    for run in RUNS:
        per_run[run] = calc_s_rmsf_single_run(system, seq_df, run, n_trim=n_trim)

    valid_vals = [v for v in per_run.values() if not np.isnan(v)]

    if not valid_vals:
        return np.nan, np.nan, per_run

    mean_val = float(np.mean(valid_vals))
    if len(valid_vals) > 1:
        sem_val = float(np.std(valid_vals, ddof=1) / np.sqrt(len(valid_vals)))
    else:
        sem_val = 0.0

    return mean_val, sem_val, per_run


# ========== Main program ==========
def main():
    print("=" * 80)
    print("Global Nucleosome S_RMSF (Per-Dyad-Pair, Standard Residue Alignment)")
    print(f"Fixed truncation: {N_TRIM} residues from each chain end (by standard numbering)")
    print(f"Runs: {RUNS}")
    print("=" * 80)

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)
    print(f"[OK] Loaded sequence file\n")

    all_results = []

    # --- WT ---
    print(">>> WT Reference")
    s_wt, sem_wt, wt_per_run = calc_s_rmsf_aggregate("WT_rw", seq_df, n_trim=N_TRIM)
    run_str = " | ".join([f"{r}: {v:.4f}" if not np.isnan(v) else f"{r}: N/A"
                          for r, v in wt_per_run.items()])
    print(f"    WT S_RMSF = {s_wt:.4f} ± {sem_wt:.4f} Å")
    print(f"    Per-run:   {run_str}\n")

    # --- Mutants ---
    for mutant_key, info in MUTANTS.items():
        resid = info["resid"]
        h_mut = info["histone_mut"]
        c_mut = info["mut_chain"]
        sys_name = f"{mutant_key}_rw"

        print(f">>> Mutant: {mutant_key} ({h_mut} {c_mut}{resid})")
        s_mut, sem_mut, mut_per_run = calc_s_rmsf_aggregate(
            sys_name, seq_df, n_trim=N_TRIM)

        delta_s = s_mut - s_wt if not (np.isnan(s_mut) or np.isnan(s_wt)) else np.nan
        delta_sem = np.sqrt(sem_mut ** 2 + sem_wt ** 2) if not (np.isnan(sem_mut) or np.isnan(sem_wt)) else np.nan

        run_str = " | ".join([f"{r}: {v:.4f}" if not np.isnan(v) else f"{r}: N/A"
                              for r, v in mut_per_run.items()])
        print(f"    Mutant S_RMSF = {s_mut:.4f} ± {sem_mut:.4f} Å")
        print(f"    Per-run:        {run_str}")

        if not np.isnan(delta_s):
            sign = "+" if delta_s > 0 else ""
            tag = "↑ Increased asymmetry" if delta_s > 0.05 else (
                "↓ Enhanced symmetry" if delta_s < -0.05 else "≈ Minimal change")
            print(f"    ΔS_RMSF = {sign}{delta_s:.4f} ± {delta_sem:.4f} Å  {tag}")
        print()

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
        }

        for i, run in enumerate(RUNS, start=1):
            wt_v = wt_per_run.get(run, np.nan)
            mut_v = mut_per_run.get(run, np.nan)
            record[f"SRMSF_WT_Run{i}_A"] = wt_v
            record[f"SRMSF_Variant_Run{i}_A"] = mut_v
            record[f"Delta_SRMSF_Run{i}_A"] = (mut_v - wt_v) if not (np.isnan(wt_v) or np.isnan(mut_v)) else np.nan

        all_results.append(record)

    # ========== Save and output ==========
    print(f"\n🔍 DEBUG: all_results contains {len(all_results)} records")

    if all_results:
        df_out = pd.DataFrame(all_results)

        OUT_DIR = BASE_DIR
        OUT_DIR.mkdir(exist_ok=True)
        out_csv = OUT_DIR / f"SRMSF_PseudoSymmetry_Deviation_Trim{N_TRIM}_PerRunDetail.csv"

        df_out.to_csv(out_csv, index=False, float_format="%.6f")
        print(f"[OK] Saved to: {out_csv.resolve()}")

        print("\n" + "=" * 80)
        print("SUMMARY TABLE")
        summary_cols = [
            "Variant", "Histone_Subunit",
            "SRMSF_WT_Mean_A", "SRMSF_WT_SEM_A",
            "SRMSF_Variant_Mean_A", "SRMSF_Variant_SEM_A",
            "Delta_SRMSF_Mean_A", "Delta_SRMSF_SEM_A"
        ]
        print(df_out[summary_cols].round(4).to_string(index=False))

        print("\nPER-RUN DETAIL")
        run_cols = ["Variant"] + [c for c in df_out.columns
                                  if c.startswith(("SRMSF_WT_Run", "SRMSF_Variant_Run", "Delta_SRMSF_Run"))]
        print(df_out[run_cols].round(4).to_string(index=False))
        print("=" * 80)
    else:
        print("\n[ERROR] all_results is empty! Check indentation of all_results.append(record)")


if __name__ == "__main__":
    main()