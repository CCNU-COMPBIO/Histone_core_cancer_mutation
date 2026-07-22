#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSF Correlation Analysis – Per‑Run Version
Only symmetry chain pairs (C vs G, D vs H, A vs E, B vs F) are considered.
Removes all self‑comparison metrics (mutant chain vs itself).
Column names are updated to publication‑ready naming.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

CHAIN_START = {
    "A": 38, "E": 37, "B": 25, "F": 21,
    "C": 14, "G": 15, "D": 30, "H": 33,
}

MUTANTS = {
    "R29Q":  {"histone": "H2A", "mut_chain": "C", "ref_chain": "G", "resid": 29},
    "E56K":  {"histone": "H2A", "mut_chain": "C", "ref_chain": "G", "resid": 56},
    "K74N":  {"histone": "H2A", "mut_chain": "C", "ref_chain": "G", "resid": 74},
    "F70L":  {"histone": "H2B", "mut_chain": "D", "ref_chain": "H", "resid": 70},
    "E50K":  {"histone": "H3",  "mut_chain": "A", "ref_chain": "E", "resid": 50},
    "E73K":  {"histone": "H3",  "mut_chain": "A", "ref_chain": "E", "resid": 73},
    "E105K": {"histone": "H3",  "mut_chain": "A", "ref_chain": "E", "resid": 105},
}


# ========== Helper Functions ==========
def read_sequence_csv(path: Path) -> pd.DataFrame:
    """Read the sequence CSV file with automatic encoding detection."""
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


def build_chain_res_to_pos(seq_df: pd.DataFrame, chain_id: str) -> pd.DataFrame:
    """
    Map residue numbers to sequential positions (1‑based) for a given chain.
    """
    sub = seq_df[seq_df["Chain"] == chain_id].copy()
    if sub.empty:
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv")
    sub = sub.sort_values("Res").reset_index(drop=True)
    sub["Pos"] = np.arange(1, len(sub) + 1)
    return sub[["Res", "Pos"]]


def load_single_run_rmsf(system, histone, chain, res_to_pos, run):
    """
    Load RMSF data for a given system, histone, chain, and simulation run.
    Returns (X, RMSF) arrays aligned to the global index (X = chain_start + Pos - 1).
    """
    f = BASE_DIR / system / run / RMSF_FILE_TEMPLATE.format(histone=histone, chain=chain)
    if not f.exists():
        return None, None
    try:
        df = pd.read_csv(f, comment="#", delim_whitespace=True, header=None,
                         names=["Res", "RMSF"], dtype={"Res": float, "RMSF": float})
    except Exception:
        return None, None

    df["Res"] = df["Res"].round().astype(int)
    merged = df.merge(res_to_pos, on="Res", how="inner")
    if merged.empty:
        return None, None
    merged["X"] = CHAIN_START[chain] + merged["Pos"] - 1
    merged = merged.sort_values("X")
    return merged["X"].to_numpy(dtype=int), merged["RMSF"].to_numpy(dtype=float)


def calc_correlation(x1, y1, x2, y2):
    """Compute Pearson and Spearman correlations on common X positions."""
    if x1 is None or x2 is None:
        return np.nan, np.nan, 0
    common = np.intersect1d(x1, x2)
    if len(common) == 0:
        return np.nan, np.nan, 0

    idx1 = {int(v): i for i, v in enumerate(x1)}
    idx2 = {int(v): i for i, v in enumerate(x2)}
    y1c = np.array([y1[idx1[int(v)]] for v in common])
    y2c = np.array([y2[idx2[int(v)]] for v in common])

    r_p, _ = pearsonr(y1c, y2c)
    r_s, _ = spearmanr(y1c, y2c)
    return r_p, r_s, len(common)


def mean_sem(values):
    """Compute mean and standard error of the mean (SEM) from a list."""
    valid = [v for v in values if not np.isnan(v)]
    if not valid:
        return np.nan, np.nan
    m = float(np.mean(valid))
    s = float(np.std(valid, ddof=1) / np.sqrt(len(valid))) if len(valid) > 1 else 0.0
    return m, s


def calc_for_single_run(run, mutant_key, info, res_maps):
    """
    Calculate all correlation metrics for one run and one mutant.
    Returns a dictionary with WT/Mutant cross‑chain correlations and deltas.
    """
    h = info["histone"]
    c_mut = info["mut_chain"]
    c_ref = info["ref_chain"]

    sys_wt = "WT_rw"
    sys_mut = f"{mutant_key}_rw"

    x_wt_m, y_wt_m = load_single_run_rmsf(sys_wt, h, c_mut, res_maps[c_mut], run)
    x_wt_r, y_wt_r = load_single_run_rmsf(sys_wt, h, c_ref, res_maps[c_ref], run)
    x_mt_m, y_mt_m = load_single_run_rmsf(sys_mut, h, c_mut, res_maps[c_mut], run)
    x_mt_r, y_mt_r = load_single_run_rmsf(sys_mut, h, c_ref, res_maps[c_ref], run)

    # WT symmetry pair
    r_xwt_p, r_xwt_s, n = calc_correlation(x_wt_m, y_wt_m, x_wt_r, y_wt_r)
    # Mutant symmetry pair
    r_xmt_p, r_xmt_s, _ = calc_correlation(x_mt_m, y_mt_m, x_mt_r, y_mt_r)

    delta_p = r_xwt_p - r_xmt_p if not (np.isnan(r_xwt_p) or np.isnan(r_xmt_p)) else np.nan
    delta_s = r_xwt_s - r_xmt_s if not (np.isnan(r_xwt_s) or np.isnan(r_xmt_s)) else np.nan

    return {
        "Cross_WT_Pearson": r_xwt_p,
        "Cross_WT_Spearman": r_xwt_s,
        "Cross_Mut_Pearson": r_xmt_p,
        "Cross_Mut_Spearman": r_xmt_s,
        "Delta_Pearson": delta_p,
        "Delta_Spearman": delta_s,
        "N_Residues": n,
    }


# ========== Main Script ==========
def main():
    print("=" * 80)
    print("RMSF Correlation — Symmetry Pairs Only (Per‑Run)")
    print(f"Runs: {RUNS}")
    print("=" * 80)

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)
    print("[OK] Loaded sequence file\n")

    needed_chains = set()
    for info in MUTANTS.values():
        needed_chains.update([info["mut_chain"], info["ref_chain"]])

    res_maps = {}
    for ch in needed_chains:
        try:
            res_maps[ch] = build_chain_res_to_pos(seq_df, ch)
        except ValueError as e:
            print(f"[WARN] {e}")

    all_results = []

    # Mapping from internal names to publication‑ready column names
    name_map = {
        "Cross_WT_Pearson": "Symmetry_WT_Pearson",
        "Cross_WT_Spearman": "Symmetry_WT_Spearman",
        "Cross_WT_SEM_P": "Symmetry_WT_Pearson_SEM",
        "Cross_WT_SEM_S": "Symmetry_WT_Spearman_SEM",
        "Cross_Mut_Pearson": "Symmetry_Mutant_Pearson",
        "Cross_Mut_Spearman": "Symmetry_Mutant_Spearman",
        "Cross_Mut_SEM_P": "Symmetry_Mutant_Pearson_SEM",
        "Cross_Mut_SEM_S": "Symmetry_Mutant_Spearman_SEM",
        "Delta_Pearson": "Symmetry_Breaking_Pearson",
        "Delta_Spearman": "Symmetry_Breaking_Spearman",
        "Delta_SEM_P": "Symmetry_Breaking_Pearson_SEM",
        "Delta_SEM_S": "Symmetry_Breaking_Spearman_SEM",
    }

    for mutant_key, info in MUTANTS.items():
        pair_label = f"{info['mut_chain']} vs {info['ref_chain']}"
        print(f">>> {mutant_key} ({info['histone']} {pair_label})")

        per_run = {run: calc_for_single_run(run, mutant_key, info, res_maps) for run in RUNS}

        metric_keys = [k for k in per_run[RUNS[0]].keys() if k != "N_Residues"]
        summary = {}
        for mk in metric_keys:
            vals = [per_run[r][mk] for r in RUNS]
            m, s = mean_sem(vals)
            summary[mk] = m
            sem_key = mk.replace("_Pearson", "_SEM_P").replace("_Spearman", "_SEM_S")
            summary[sem_key] = s

        print(f"    Cross‑WT:        {summary['Cross_WT_Pearson']:.4f} ± {summary['Cross_WT_SEM_P']:.4f}")
        print(f"    Cross‑Mutant:    {summary['Cross_Mut_Pearson']:.4f} ± {summary['Cross_Mut_SEM_P']:.4f}")
        print(f"    ΔSymmetry:       {summary['Delta_Pearson']:.4f} ± {summary['Delta_SEM_P']:.4f}")
        print()

        record = {
            "Mutant": mutant_key,
            "Histone": info["histone"],
            "Chain_Pair": pair_label,
        }

        # Add mean values and SEMs using new column names
        for old_key, val in summary.items():
            if old_key in name_map:
                record[name_map[old_key]] = val

        # Add per‑run values for non‑delta metrics
        for run in RUNS:
            for old_key in metric_keys:
                if not old_key.startswith("Delta"):
                    new_key = name_map[old_key]  # must exist
                    record[f"{new_key}_{run}"] = per_run[run][old_key]

        # Number of common residues (same for all runs)
        record["N_Common_Residues"] = per_run[RUNS[0]]["N_Residues"]

        all_results.append(record)

    if all_results:
        df_out = pd.DataFrame(all_results)
        out_csv = BASE_DIR / "RMSF_correlation_symmetry_pairs.csv"
        df_out.to_csv(out_csv, index=False, float_format="%.6f", encoding="utf-8-sig")
        print(f"[OK] Saved to: {out_csv}\n")

        print("=" * 80)
        print("SUMMARY (Mean Pearson r ± SEM)")
        cols = ["Mutant", "Chain_Pair",
                "Symmetry_WT_Pearson", "Symmetry_WT_Pearson_SEM",
                "Symmetry_Mutant_Pearson", "Symmetry_Mutant_Pearson_SEM",
                "Symmetry_Breaking_Pearson", "Symmetry_Breaking_Pearson_SEM"]
        print(df_out[cols].round(4).to_string(index=False))
        print("=" * 80)
    else:
        print("[WARN] No results collected.")


if __name__ == "__main__":
    main()