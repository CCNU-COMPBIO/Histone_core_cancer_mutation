#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-Pair RMSF-based Pseudo-symmetry Deviation Score (S_RMSF)
Calculated for a single histone dimer chain pair (mutant chain vs. its copy).
Output: mean + SEM + per-run independent values.
Alignment: use CHAIN_START to map residues to canonical histone numbering, align by std.
Trimming: symmetrically remove N_TRIM residues from both ends at shared std positions.
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

N_TRIM = 5   # Number of residues to trim from each end

# Canonical histone numbering corresponding to the first residue of each chain
CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

# Chain pair definitions (mut chain is the mutant chain, non chain is its copy)
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
    std = chain_start + (Res - first Res of this chain)
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
                         res_to_std: dict, run: str) -> dict:
    """
    Load RMSF for a single chain and single run; return {canonical std: RMSF}.
    Note: do not trim here; trimming is performed uniformly on common_std.
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
    return result


# ---------- Calculation for a single chain pair ----------
def calc_s_rmsf_single_pair(system: str, seq_df: pd.DataFrame,
                            pair: dict, run: str, n_trim: int = 5) -> float:
    """
    Compute S_RMSF for a single chain pair (one simulation run).
    First take shared std positions, then symmetrically trim both ends at shared positions.
    """
    mut_info = pair["mut"]
    non_info = pair["non"]

    try:
        mut_map = build_res_to_std(
            seq_df, mut_info["chain"], CHAIN_START[mut_info["chain"]])
        non_map = build_res_to_std(
            seq_df, non_info["chain"], CHAIN_START[non_info["chain"]])
    except ValueError:
        return np.nan

    rmsf_mut = load_single_run_rmsf(
        system, mut_info["chain"], mut_info["histone"], mut_map, run)
    rmsf_non = load_single_run_rmsf(
        system, non_info["chain"], non_info["histone"], non_map, run)

    common_std = sorted(set(rmsf_mut.keys()) & set(rmsf_non.keys()))
    if not common_std:
        return np.nan

    # Symmetrically truncate at shared homologous positions
    if len(common_std) > 2 * n_trim:
        common_std = common_std[n_trim:-n_trim]
    else:
        print(f"      [WARN] {mut_info['chain']}-{non_info['chain']} "
              f"only {len(common_std)} common std positions, "
              f"cannot trim {n_trim} from each end, skip this pair")
        return np.nan

    abs_diff_sum = sum(abs(rmsf_mut[s] - rmsf_non[s]) for s in common_std)
    return abs_diff_sum / len(common_std)


def calc_s_rmsf_aggregate_pair(system: str, seq_df: pd.DataFrame,
                               pair: dict, n_trim: int = 5):
    """
    Aggregate across runs for a single chain pair (mean + SEM), and also return
    per-run independent values.
    Returns: (mean, sem, per_run_dict)
    """
    per_run = {}
    for run in RUNS:
        per_run[run] = calc_s_rmsf_single_pair(system, seq_df, pair, run, n_trim=n_trim)

    valid_vals = [v for v in per_run.values() if not np.isnan(v)]
    if not valid_vals:
        return np.nan, np.nan, per_run

    mean_val = float(np.mean(valid_vals))
    sem_val = float(np.std(valid_vals, ddof=1) / np.sqrt(len(valid_vals))) if len(valid_vals) > 1 else 0.0
    return mean_val, sem_val, per_run


# ========== Main program ==========
def main():
    print("=" * 80)
    print("Single-Pair S_RMSF (Mutant chain vs. its copy, Standard Residue Alignment)")
    print(f"Fixed truncation: {N_TRIM} residues from each end (by shared standard numbering)")
    print(f"Runs: {RUNS}")
    print("=" * 80)

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)
    print(f"[OK] Loaded sequence file\n")

    chain_to_pair = {}
    for pair in DYAD_PAIRS:
        mut_chain = pair["mut"]["chain"]
        chain_to_pair[mut_chain] = pair

    all_records = []

    for mutant_key, info in MUTANTS.items():
        mut_chain = info["mut_chain"]
        histone = info["histone_mut"]
        resid = info["resid"]
        sys_mut = f"{mutant_key}_rw"
        sys_wt = "WT_rw"

        pair = chain_to_pair.get(mut_chain)
        if pair is None:
            print(f"[WARN] No chain pair found for chain {mut_chain}, skipping {mutant_key}")
            continue

        mut_chain_label = pair["mut"]["chain"]
        non_chain_label = pair["non"]["chain"]
        print(f">>> Mutant: {mutant_key} ({histone} {mut_chain_label}/{non_chain_label})")

        wt_mean, wt_sem, wt_per_run = calc_s_rmsf_aggregate_pair(
            sys_wt, seq_df, pair, n_trim=N_TRIM)

        mut_mean, mut_sem, mut_per_run = calc_s_rmsf_aggregate_pair(
            sys_mut, seq_df, pair, n_trim=N_TRIM)

        if not (np.isnan(wt_mean) or np.isnan(mut_mean)):
            delta_mean = mut_mean - wt_mean
            delta_sem = np.sqrt(wt_sem**2 + mut_sem**2) if (wt_sem is not None and mut_sem is not None) else np.nan
        else:
            delta_mean, delta_sem = np.nan, np.nan

        print(f"    WT S_RMSF     = {wt_mean:.4f} ± {wt_sem:.4f} Å")
        print(f"    Mutant S_RMSF = {mut_mean:.4f} ± {mut_sem:.4f} Å")
        if not np.isnan(delta_mean):
            print(f"    ΔS_RMSF       = {delta_mean:+.4f} ± {delta_sem:.4f} Å")
        print()

        record = {
            "Variant": mutant_key,
            "Histone": histone,
            "Chain_Mut": mut_chain_label,
            "Chain_Copy": non_chain_label,
            "Residue": resid,
            "SRMSF_WT_Mean_A": wt_mean,
            "SRMSF_WT_SEM_A": wt_sem,
            "SRMSF_Mut_Mean_A": mut_mean,
            "SRMSF_Mut_SEM_A": mut_sem,
            "Delta_SRMSF_Mean_A": delta_mean,
            "Delta_SRMSF_SEM_A": delta_sem,
        }
        for i, run in enumerate(RUNS, start=1):
            record[f"SRMSF_WT_Run{i}_A"] = wt_per_run.get(run, np.nan)
            record[f"SRMSF_Mut_Run{i}_A"] = mut_per_run.get(run, np.nan)
            wt_v = wt_per_run.get(run, np.nan)
            mut_v = mut_per_run.get(run, np.nan)
            record[f"Delta_SRMSF_Run{i}_A"] = (mut_v - wt_v) if not (np.isnan(wt_v) or np.isnan(mut_v)) else np.nan

        all_records.append(record)

    if all_records:
        df_out = pd.DataFrame(all_records)
        out_csv = BASE_DIR / f"SRMSF_SinglePair_Trim{N_TRIM}_PerRunDetail.csv"
        df_out.to_csv(out_csv, index=False, float_format="%.6f")
        print(f"[OK] Saved to: {out_csv.resolve()}")

        print("\n" + "=" * 80)
        print("SUMMARY TABLE")
        summary_cols = [
            "Variant", "Histone", "Chain_Mut", "Chain_Copy",
            "SRMSF_WT_Mean_A", "SRMSF_WT_SEM_A",
            "SRMSF_Mut_Mean_A", "SRMSF_Mut_SEM_A",
            "Delta_SRMSF_Mean_A", "Delta_SRMSF_SEM_A"
        ]
        print(df_out[summary_cols].round(4).to_string(index=False))

        print("\nPER-RUN DETAIL")
        run_cols = ["Variant"] + [c for c in df_out.columns
                                  if c.startswith(("SRMSF_WT_Run", "SRMSF_Mut_Run", "Delta_SRMSF_Run"))]
        print(df_out[run_cols].round(4).to_string(index=False))
        print("=" * 80)
    else:
        print("[ERROR] No records generated.")


if __name__ == "__main__":
    main()