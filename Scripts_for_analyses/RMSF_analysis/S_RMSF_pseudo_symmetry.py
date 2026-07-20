#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
srmsf_pseudo_symmetry.py

Compute RMSF-based pseudo-symmetry deviation score (S_RMSF) for nucleosome systems.
The score measures the mean absolute RMSF difference between dyad-related homologous
residues across the two halves of the nucleosome. Both mean and per-run values are
reported, along with delta (mutant – WT) and standard error of the mean (SEM).
N-terminal and C-terminal residues are trimmed to avoid tail-dominated signals.
"""

from pathlib import Path
import numpy as np
import pandas as pd

# ========== CONFIGURATION ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

N_TRIM = 5  # Number of residues trimmed from each end of each chain

# Dyad-related chain pairs (mutated half chain → non‑mutated half chain)
DYAD_PAIRS = [
    {"mut": {"chain": "A", "histone": "H3"}, "non": {"chain": "E", "histone": "H3"}},
    {"mut": {"chain": "B", "histone": "H4"}, "non": {"chain": "F", "histone": "H4"}},
    {"mut": {"chain": "C", "histone": "H2A"}, "non": {"chain": "G", "histone": "H2A"}},
    {"mut": {"chain": "D", "histone": "H2B"}, "non": {"chain": "H", "histone": "H2B"}},
]

# Mutant definitions
MUTANTS = {
    "R29Q": {"histone_mut": "H2A", "mut_chain": "C", "resid": 29},
    "E56K": {"histone_mut": "H2A", "mut_chain": "C", "resid": 56},
    "K74N": {"histone_mut": "H2A", "mut_chain": "C", "resid": 74},
    "F70L": {"histone_mut": "H2B", "mut_chain": "D", "resid": 70},
    "E50K": {"histone_mut": "H3", "mut_chain": "A", "resid": 50},
    "E73K": {"histone_mut": "H3", "mut_chain": "A", "resid": 73},
    "E105K": {"histone_mut": "H3", "mut_chain": "A", "resid": 105},
}


# ========== HELPER FUNCTIONS ==========
def read_sequence_csv(path: Path) -> pd.DataFrame:
    """Read sequence CSV file, handling various encodings and column names."""
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


def build_res_to_pos(seq_df: pd.DataFrame, chain_id: str) -> dict:
    """
    Build a mapping from residue number (Res) to sequential position (Pos)
    for a given chain, using the order in the sequence file.
    """
    sub = seq_df[seq_df["Chain"] == chain_id].sort_values("Res").reset_index(drop=True)
    if sub.empty:
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv")
    return dict(zip(sub["Res"].astype(int), range(1, len(sub) + 1)))


def trim_fixed_ends(pos_rmsf: dict, n_trim: int = 5) -> dict:
    """Remove the first and last n_trim positions from the RMSF dictionary."""
    if not pos_rmsf:
        return {}
    positions = sorted(pos_rmsf.keys())
    if len(positions) <= 2 * n_trim:
        print(f"      [WARN] Chain too short ({len(positions)} residues) "
              f"for trimming {n_trim} from each end, returning empty")
        return {}
    trimmed_positions = positions[n_trim:-n_trim]
    return {p: pos_rmsf[p] for p in trimmed_positions}


def load_single_run_rmsf(system: str, chain: str, histone: str,
                         res_to_pos: dict, run: str, n_trim: int = 5) -> dict:
    """
    Load RMSF data for a single chain and run, convert to positions,
    trim ends, and return as {Pos: RMSF}.
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
        pos = res_to_pos.get(row["Res"])
        if pos is not None and not np.isnan(row["RMSF"]):
            result[pos] = float(row["RMSF"])

    return trim_fixed_ends(result, n_trim=n_trim)


def calc_s_rmsf_single_run(system: str, seq_df: pd.DataFrame,
                           run: str, n_trim: int = 5) -> float:
    """
    Compute S_RMSF for a single system and run.
    S_RMSF = (1/N) * sum_i | RMSF_i(mutated half) - RMSF_m(i)(non-mutated half) |
    where i runs over all matched residue positions across all dyad pairs.
    Returns NaN if no valid residue pairs are found.
    """
    total_abs_diff_sum = 0.0
    total_n = 0

    for pair in DYAD_PAIRS:
        mut_info = pair["mut"]
        non_info = pair["non"]

        try:
            mut_map = build_res_to_pos(seq_df, mut_info["chain"])
            non_map = build_res_to_pos(seq_df, non_info["chain"])
        except ValueError:
            continue

        rmsf_mut = load_single_run_rmsf(
            system, mut_info["chain"], mut_info["histone"],
            mut_map, run, n_trim=n_trim)
        rmsf_non = load_single_run_rmsf(
            system, non_info["chain"], non_info["histone"],
            non_map, run, n_trim=n_trim)

        common_pos = set(rmsf_mut.keys()) & set(rmsf_non.keys())
        if common_pos:
            for p in common_pos:
                total_abs_diff_sum += abs(rmsf_mut[p] - rmsf_non[p])
            total_n += len(common_pos)

    return total_abs_diff_sum / total_n if total_n > 0 else np.nan


def calc_s_rmsf_aggregate(system: str, seq_df: pd.DataFrame, n_trim: int = 5):
    """
    Compute S_RMSF for all runs of a given system.
    Returns (mean, SEM, per_run_dict) where per_run_dict maps run name to S_RMSF.
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


# ========== MAIN ==========
def main():
    print("=" * 80)
    print("Global Nucleosome S_RMSF (Dyad-Pair Based, Per-Run Output)")
    print(f"Fixed truncation: {N_TRIM} residues from each chain end")
    print(f"Runs: {RUNS}")
    print("=" * 80)

    seq_df = read_sequence_csv(BASE_DIR / SEQUENCE_CSV)
    print(f"[OK] Loaded sequence file\n")

    all_results = []

    # --- Wild-type reference ---
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

        # Build record with descriptive column names
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

        # Add per-run values
        for i, run in enumerate(RUNS, start=1):
            wt_v = wt_per_run.get(run, np.nan)
            mut_v = mut_per_run.get(run, np.nan)
            record[f"SRMSF_WT_Run{i}_A"] = wt_v
            record[f"SRMSF_Variant_Run{i}_A"] = mut_v
            record[f"Delta_SRMSF_Run{i}_A"] = (mut_v - wt_v) if not (np.isnan(wt_v) or np.isnan(mut_v)) else np.nan

        all_results.append(record)

    # ========== SAVE AND SUMMARY ==========
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
        print("\n[ERROR] No results collected.")


if __name__ == "__main__":
    main()