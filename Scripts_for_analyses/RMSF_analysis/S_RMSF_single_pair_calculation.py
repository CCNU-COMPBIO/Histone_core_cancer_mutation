#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-Pair RMSF-based Pseudo-symmetry Deviation Score (S_RMSF)
Computed for a single histone dimer chain pair (mutant chain vs its copy chain).
Output: mean + SEM + per-run values.
Alignment: use CHAIN_START to map to canonical histone numbering, align by std.
Trimming: symmetrically remove N_TRIM residues from both ends on the shared std positions.

Folder names are already unified as <histone>_<mutation>_rw, e.g. H2A_R29Q_rw, H2B_E105K_rw, H3_E50K_rw.
Therefore the MUTANTS keys use that naming directly (with the "_rw" suffix removed).
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

N_TRIM = 5   # number of residues to trim from each end

WT_SYSTEM = "WT_rw"
NOEXCL = "_noexcl"
VARIANTS = [
    {"suffix": "",     "label": "excluded segments"},
    {"suffix": NOEXCL, "label": "no excluded segments"},
]

# Canonical histone numbering assigned to the first residue of each chain
CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

# Chain pair definitions (mut chain is the mutant chain, non chain is its copy)
DYAD_PAIRS = [
    {"mut": {"chain": "A", "histone": "H3"},  "non": {"chain": "E", "histone": "H3"}},
    {"mut": {"chain": "B", "histone": "H4"},  "non": {"chain": "F", "histone": "H4"}},
    {"mut": {"chain": "C", "histone": "H2A"}, "non": {"chain": "G", "histone": "H2A"}},
    {"mut": {"chain": "D", "histone": "H2B"}, "non": {"chain": "H", "histone": "H2B"}},
]

# Keys match the folder names (without "_rw")
#   Folder lookup: f"{mutant_key}_rw" -> H2A_R29Q_rw / H2B_E105K_rw / H3_E50K_rw ...
MUTANTS = {
    "H2A_R29Q":  {"histone_mut": "H2A", "mut_chain": "C", "resid": 29},
    "H2A_E56K":  {"histone_mut": "H2A", "mut_chain": "C", "resid": 56},
    "H2A_K74N":  {"histone_mut": "H2A", "mut_chain": "C", "resid": 74},
    "H2B_E105K": {"histone_mut": "H2B", "mut_chain": "D", "resid": 105},   # newly added
    "H2B_F70L":  {"histone_mut": "H2B", "mut_chain": "D", "resid": 70},
    "H3_E50K":   {"histone_mut": "H3",  "mut_chain": "A", "resid": 50},
    "H3_E73K":   {"histone_mut": "H3",  "mut_chain": "A", "resid": 73},
    "H3_E105K":  {"histone_mut": "H3",  "mut_chain": "A", "resid": 105},
}


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
    Build the {original Res: canonical std numbering} mapping for the chain.
    std = chain_start + (Res - first Res of that chain)
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


def load_single_run_rmsf(system: str, chain: str, histone: str,
                         res_to_std: dict, run: str, suffix: str = "") -> dict:
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


# ---------- Single chain-pair computation ----------
def calc_s_rmsf_single_pair(system: str, seq_df: pd.DataFrame, pair: dict,
                            run: str, n_trim: int = 5, suffix: str = "") -> float:
    mut_info, non_info = pair["mut"], pair["non"]
    try:
        mut_map = build_res_to_std(seq_df, mut_info["chain"], CHAIN_START[mut_info["chain"]])
        non_map = build_res_to_std(seq_df, non_info["chain"], CHAIN_START[non_info["chain"]])
    except ValueError:
        return np.nan

    rmsf_mut = load_single_run_rmsf(system, mut_info["chain"], mut_info["histone"],
                                    mut_map, run, suffix)
    rmsf_non = load_single_run_rmsf(system, non_info["chain"], non_info["histone"],
                                    non_map, run, suffix)

    common_std = sorted(set(rmsf_mut.keys()) & set(rmsf_non.keys()))
    if not common_std:
        return np.nan
    if len(common_std) > 2 * n_trim:
        common_std = common_std[n_trim:-n_trim]
    else:
        print(f"      [WARN] {mut_info['chain']}-{non_info['chain']} only "
              f"{len(common_std)} common std positions, cannot trim, skip")
        return np.nan

    return sum(abs(rmsf_mut[s] - rmsf_non[s]) for s in common_std) / len(common_std)


def calc_s_rmsf_aggregate_pair(system: str, seq_df: pd.DataFrame, pair: dict,
                               runs, n_trim: int = 5, suffix: str = ""):
    per_run = {run: calc_s_rmsf_single_pair(system, seq_df, pair, run,
                                            n_trim=n_trim, suffix=suffix)
               for run in runs}
    valid = [v for v in per_run.values() if not np.isnan(v)]
    if not valid:
        return np.nan, np.nan, per_run
    mean_val = float(np.mean(valid))
    sem_val = (float(np.std(valid, ddof=1) / np.sqrt(len(valid))) if len(valid) > 1 else 0.0)
    return mean_val, sem_val, per_run


# ========== Main ==========
def main():
    noexcl_map = scan_noexcl()
    wt_runs_b = sorted({r for rs in noexcl_map.values() for r in rs})

    print("=" * 80)
    print("[Variant B scan] Which systems/runs have _noexcl files")
    for s, rs in (noexcl_map.items() if noexcl_map else []):
        print(f"  {s:12s} -> {rs}")
    if noexcl_map:
        print(f"  WT runs used in Variant B = {wt_runs_b}")
    print("=" * 80)

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)

    chain_to_pair = {p["mut"]["chain"]: p for p in DYAD_PAIRS}
    results = {}

    for variant in VARIANTS:
        suffix, label = variant["suffix"], variant["label"]
        print(f"\n{'#'*80}\n#  VARIANT: {label}   (suffix = '{suffix or 'none'}')\n{'#'*80}")

        rows = []
        for mutant_key, info in MUTANTS.items():
            sys_mut = f"{mutant_key}_rw"          # e.g. H2B_E105K_rw

            if suffix and sys_mut not in noexcl_map:
                print(f">>> {mutant_key}: no _noexcl files -> skipped in this variant")
                continue

            mut_chain = info["mut_chain"]
            pair = chain_to_pair.get(mut_chain)
            if pair is None:
                continue

            use_runs = noexcl_map.get(sys_mut, RUNS) if suffix else list(RUNS)
            wt_runs = wt_runs_b if suffix else list(RUNS)

            wt_mean, wt_sem, wt_per_run = calc_s_rmsf_aggregate_pair(
                WT_SYSTEM, seq_df, pair, wt_runs, n_trim=N_TRIM, suffix=suffix)
            mut_mean, mut_sem, mut_per_run = calc_s_rmsf_aggregate_pair(
                sys_mut, seq_df, pair, use_runs, n_trim=N_TRIM, suffix=suffix)

            delta_mean = (mut_mean - wt_mean
                          if not (np.isnan(wt_mean) or np.isnan(mut_mean)) else np.nan)
            delta_sem = (np.sqrt(wt_sem ** 2 + mut_sem ** 2)
                         if not (np.isnan(wt_sem) or np.isnan(mut_sem)) else np.nan)

            print(f">>> {mutant_key} ({info['histone_mut']} "
                  f"{pair['mut']['chain']}/{pair['non']['chain']}) runs={use_runs}")
            print(f"    WT S_RMSF = {wt_mean:.4f} ± {wt_sem:.4f} | "
                  f"Mut S_RMSF = {mut_mean:.4f} ± {mut_sem:.4f} | "
                  f"Δ = {delta_mean:+.4f}")

            record = {
                "Variant": mutant_key,
                "Histone": info["histone_mut"],
                "Chain_Mut": pair["mut"]["chain"],
                "Chain_Copy": pair["non"]["chain"],
                "Residue": info["resid"],
                "SRMSF_WT_Mean_A": wt_mean,
                "SRMSF_WT_SEM_A": wt_sem,
                "SRMSF_Mut_Mean_A": mut_mean,
                "SRMSF_Mut_SEM_A": mut_sem,
                "Delta_SRMSF_Mean_A": delta_mean,
                "Delta_SRMSF_SEM_A": delta_sem,
                "Runs_Used": ",".join(use_runs),
            }
            for i, run in enumerate(RUNS, start=1):
                wt_v = wt_per_run.get(run, np.nan)
                mut_v = mut_per_run.get(run, np.nan)
                record[f"SRMSF_WT_Run{i}_A"] = wt_v
                record[f"SRMSF_Mut_Run{i}_A"] = mut_v
                record[f"Delta_SRMSF_Run{i}_A"] = (mut_v - wt_v
                                                   if not (np.isnan(wt_v) or np.isnan(mut_v))
                                                   else np.nan)
            rows.append(record)

        df = pd.DataFrame(rows)
        out_csv = BASE_DIR / f"SRMSF_SinglePair_Trim{N_TRIM}_PerRunDetail{suffix}.csv"
        df.to_csv(out_csv, index=False, float_format="%.6f")
        print(f"[OK] {out_csv.name}  rows={len(df)}")
        results[suffix] = df

    # ---- Comparison of the two variants (for filling in the local ΔS_RMSF column of Table S17) ----
    if "" in results and NOEXCL in results and len(results[""]) and len(results[NOEXCL]):
        cols = ["Variant", "SRMSF_WT_Mean_A", "SRMSF_Mut_Mean_A", "Delta_SRMSF_Mean_A"]
        m = results[""][cols].merge(results[NOEXCL][cols], on="Variant",
                                    suffixes=("_excl", "_noexcl"))
        m["dSRMSF_diff"] = (m["Delta_SRMSF_Mean_A_noexcl"]
                            - m["Delta_SRMSF_Mean_A_excl"])
        cmp_csv = BASE_DIR / "SRMSF_singlepair_compare_excl_vs_noexcl.csv"
        m.to_csv(cmp_csv, index=False, float_format="%.6f")
        print("\n" + "=" * 96)
        print("[Comparison] Chain-pair-level S_RMSF: excluded vs no-exclusion")
        print("=" * 96)
        print(m.round(4).to_string(index=False))
        print(f"\n[OK] {cmp_csv.name}")


if __name__ == "__main__":
    main()