#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare RMSF profiles between mutant systems and wild-type for specific histone chains.
Generates line plots with standard deviation bands and exports delta tables.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import seaborn as sns
sns.set_theme(style="whitegrid", context="paper")

import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# ========== PATHS AND RUN PARAMETERS ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"
WT_SYSTEM = "WT_rw"

CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

# ========== MUTANT SYSTEM MAPPING ==========
MUTATION_MAP = {
    "R29Q":  {"histone": "H2A", "chain": "C", "x_resnum": 29},
    "E56K":  {"histone": "H2A", "chain": "C", "x_resnum": 56},
    "K74N":  {"histone": "H2A", "chain": "C", "x_resnum": 74},
    "F70L":  {"histone": "H2B", "chain": "D", "x_resnum": 70},
    "E50K":  {"histone": "H3",  "chain": "A", "x_resnum": 50},
    "E73K":  {"histone": "H3",  "chain": "A", "x_resnum": 73},
    "E105K": {"histone": "H3",  "chain": "A", "x_resnum": 105},
}

# Fixed chain pairs: H3 mutants plot A+B, H2A/H2B mutants plot C+D
HISTONE_CHAIN_PAIRS = {
    "H3":  ("A", "B"),
    "H2A": ("C", "D"),
    "H2B": ("C", "D"),
}

# Actual histone type for each chain (used to build RMSF filenames)
CHAIN_TO_HISTONE = {
    "A": "H3", "B": "H4",
    "C": "H2A", "D": "H2B",
    "E": "H3", "F": "H4",
    "G": "H2A", "H": "H2B",
}

PLOT_STD_BAND = True
DPI = 600
DIFF_THRESHOLD = 0.5


# ========== UTILITY FUNCTIONS ==========
def read_sequence_csv(path: Path) -> pd.DataFrame:
    encodings = ["utf-8", "utf-8-sig", "gbk", "cp936", "latin1"]
    last = None
    for enc in encodings:
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError as e:
            last = e
            df = None
    if df is None:
        raise UnicodeDecodeError("unknown", b"", 0, 1, f"Failed to decode {path}; last={last}")

    rename = {}
    if "#Res" in df.columns: rename["#Res"] = "Res"
    if "chain" in df.columns: rename["chain"] = "Chain"
    if "residue" in df.columns: rename["residue"] = "ResidueName"
    if "#Orig" in df.columns: rename["#Orig"] = "Orig"
    df = df.rename(columns=rename)

    if "Res" not in df.columns or "Chain" not in df.columns:
        raise ValueError(f"{path} must contain #Res/Res and chain/Chain. Got: {list(df.columns)}")

    df["Res"] = pd.to_numeric(df["Res"], errors="coerce")
    df = df.dropna(subset=["Res"]).copy()
    df["Res"] = df["Res"].astype(int)
    df["Chain"] = df["Chain"].astype(str).str.strip()
    return df


def build_chain_res_to_pos(seq_df: pd.DataFrame, chain_id: str) -> pd.DataFrame:
    sub = seq_df[seq_df["Chain"] == chain_id].copy()
    if sub.empty:
        chains = sorted(seq_df["Chain"].unique().tolist())
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv. Available: {chains}")
    sub = sub.sort_values("Res").reset_index(drop=True)
    sub["Pos"] = np.arange(1, len(sub) + 1)
    return sub[["Res", "Pos"]]


def read_cpptraj_rmsf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", delim_whitespace=True, header=None,
                     names=["Res", "RMSF"], dtype={"Res": float, "RMSF": float})
    if df.empty:
        raise ValueError(f"{path} parsed as empty dataframe.")
    df["Res"] = df["Res"].round().astype(int)
    return df


def load_system_chain_runs(system, file_histone, chain, res_to_pos, start_index):
    """
    file_histone: histone name used to construct the RMSF filename
                  (may differ from the mutation's histone)
    """
    per_run = []
    for run in RUNS:
        f = BASE_DIR / system / run / RMSF_FILE_TEMPLATE.format(histone=file_histone, chain=chain)
        if not f.exists():
            raise FileNotFoundError(f"Missing file: {f}")
        df = read_cpptraj_rmsf(f)
        merged = df.merge(res_to_pos, on="Res", how="inner")
        if merged.empty:
            raise ValueError(f"No residues matched between {f} and sequence.csv for chain={chain}.")
        merged["X"] = start_index + merged["Pos"] - 1
        merged = merged.sort_values("X")
        per_run.append((merged["X"].to_numpy(dtype=int), merged["RMSF"].to_numpy(dtype=float)))
    return per_run


def align_and_stack(series_list):
    sets = [set(x.tolist()) for x, _ in series_list]
    common = sorted(set.intersection(*sets))
    if not common:
        raise ValueError("No common X positions across runs.")
    idx = {v: i for i, v in enumerate(common)}
    Y = np.full((len(series_list), len(common)), np.nan, dtype=float)
    for i, (x, y) in enumerate(series_list):
        for xv, yv in zip(x, y):
            j = idx.get(int(xv))
            if j is not None:
                Y[i, j] = float(yv)
    valid = ~np.isnan(Y).any(axis=0)
    return np.array(common)[valid], Y[:, valid]


def summarize_runs(per_run):
    common_x, Y = align_and_stack(per_run)
    mean = Y.mean(axis=0)
    std = Y.std(axis=0, ddof=1) if Y.shape[0] > 1 else np.zeros_like(mean)
    return common_x, mean, std


def compute_system_mean_for_chain(system, chain_id, res_to_pos):
    """Determine the histone type for the RMSF filename based on the chain."""
    file_histone = CHAIN_TO_HISTONE[chain_id]
    per_run = load_system_chain_runs(system, file_histone, chain_id, res_to_pos, CHAIN_START[chain_id])
    return summarize_runs(per_run)


def x_to_pos(chain_start, x_resnum):
    return int(x_resnum - chain_start + 1)


def annotate_single_mutation(ax, mut_name, mut_info, chain_id, res_to_pos, y_map):
    if mut_info["chain"] != chain_id:
        return
    chain_start = CHAIN_START[chain_id]
    x_mut = int(mut_info["x_resnum"])
    pos_mut = x_to_pos(chain_start, x_mut)

    if pos_mut < 1 or not (res_to_pos["Pos"] == pos_mut).any():
        print(f"[WARN] Mutation {mut_name} (X={x_mut}) not in chain {chain_id}")
        return
    if x_mut not in y_map:
        print(f"[WARN] Mutation {mut_name} X={x_mut} not found on mean curve for chain {chain_id}")
        return

    ax.plot([x_mut], [y_map[x_mut]],
            marker="o", markersize=10,
            markerfacecolor="crimson", markeredgecolor="white",
            markeredgewidth=1.2, linestyle="None", zorder=10)


def diff_table_for_chain(chain_id, res_to_pos, sys_a, sys_b):
    x_a, mean_a, _ = compute_system_mean_for_chain(sys_a, chain_id, res_to_pos)
    x_b, mean_b, _ = compute_system_mean_for_chain(sys_b, chain_id, res_to_pos)
    common = np.intersect1d(x_a, x_b)
    if common.size == 0:
        raise ValueError(f"No common X between {sys_a} and {sys_b} for chain {chain_id}")

    idx_a = {int(x): i for i, x in enumerate(x_a)}
    idx_b = {int(x): i for i, x in enumerate(x_b)}
    rows = []
    for x in common:
        wt = float(mean_a[idx_a[int(x)]])
        mut = float(mean_b[idx_b[int(x)]])
        diff = mut - wt
        rows.append({
            "Chain": chain_id, "X": int(x),
            f"{sys_a}_mean": wt, f"{sys_b}_mean": mut,
            "diff": diff, "abs_diff": abs(diff),
        })
    return pd.DataFrame(rows).sort_values("X").reset_index(drop=True)


def print_diff_hits(df, threshold, label):
    hits = df[df["RMSF_delta_abs"] >= float(threshold)].copy()
    if hits.empty:
        print(f"[INFO] [{label}] No sites with |RMSF_delta| >= {threshold}")
        return
    hits = hits.sort_values("RMSF_delta_abs", ascending=False)
    print(f"\n[DIFF] [{label}] Sites with |RMSF_delta| >= {threshold}")
    for _, r in hits.iterrows():
        print(f"  Chain {r['chain']}  resid={int(r['resid'])}  RMSF_delta={r['RMSF_delta']:.3f}  (abs={r['RMSF_delta_abs']:.3f})")


# ========== MAIN PIPELINE ==========
def main():
    seq_path = BASE_DIR / SEQUENCE_CSV
    if not seq_path.exists():
        raise FileNotFoundError(f"Cannot find: {seq_path}")
    seq_df = read_sequence_csv(seq_path)

    pal = sns.color_palette("deep", 2)
    color_wt, color_mut = pal[0], pal[1]

    for mut_name, mut_info in MUTATION_MAP.items():
        histone = mut_info["histone"]
        mut_system = f"{mut_name}_rw"
        display_label = f"{histone}{mut_name}"

        if histone not in HISTONE_CHAIN_PAIRS:
            print(f"[WARN] No chain pair defined for {histone}, skipping {mut_name}.")
            continue

        # Fixed chain pairs: H3->(A,B), H2A/H2B->(C,D)
        chain1, chain2 = HISTONE_CHAIN_PAIRS[histone]
        print(f"\n{'='*60}")
        print(f"Processing: {display_label} ({mut_system}) | chains: {chain1}, {chain2}")
        print(f"  (File histone mapping: {chain1}->{CHAIN_TO_HISTONE[chain1]}, {chain2}->{CHAIN_TO_HISTONE[chain2]})")
        print(f"{'='*60}")

        chain_to_res_to_pos = {}
        for cid in (chain1, chain2):
            try:
                chain_to_res_to_pos[cid] = build_chain_res_to_pos(seq_df, cid)
            except ValueError as e:
                print(f"[ERROR] {e}")

        for chain_id in (chain1, chain2):
            if chain_id not in chain_to_res_to_pos:
                print(f"[SKIP] Chain {chain_id} not available")
                continue

            res_to_pos = chain_to_res_to_pos[chain_id]
            file_histone = CHAIN_TO_HISTONE[chain_id]

            # --- Difference statistics ---
            try:
                df_diff = diff_table_for_chain(chain_id, res_to_pos,
                                               sys_a=WT_SYSTEM, sys_b=mut_system)

                # Rename columns to standard names
                df_diff.rename(columns={
                    "Chain": "chain",
                    "X": "resid",
                    f"{WT_SYSTEM}_mean": "RMSF_WT",
                    f"{mut_system}_mean": f"RMSF_{histone}{mut_name}",
                    "diff": "RMSF_delta",
                    "abs_diff": "RMSF_delta_abs"
                }, inplace=True)

                # Print delta hits
                print_diff_hits(df_diff, DIFF_THRESHOLD, display_label)

                # Save CSV (4 decimals)
                out_csv = BASE_DIR / f"RMSF_delta_{histone}{mut_name}_chain{chain_id}.csv"
                df_diff.to_csv(out_csv, index=False, float_format='%.4f')
                print(f"[OK] Saved diff table: {out_csv}")
            except Exception as e:
                print(f"[ERROR] Diff failed for {display_label} chain {chain_id}: {e}")
                continue

            # --- Plotting ---
            fig, ax = plt.subplots(figsize=(8, 4))
            mean_by_system = {}

            for system, color, ls, legend_label in [
                (WT_SYSTEM, color_wt, "-", "WT"),
                (mut_system, color_mut, "--", display_label),
            ]:
                try:
                    x, mean, std = compute_system_mean_for_chain(
                        system, chain_id, res_to_pos)
                except Exception as e:
                    print(f"[ERROR] Cannot load {system} chain {chain_id} (file: {file_histone}): {e}")
                    continue

                mean_by_system[system] = (x, mean)
                ax.plot(x, mean, color=color, linestyle=ls, linewidth=2.0,
                        label=legend_label)

                if PLOT_STD_BAND and len(RUNS) > 1:
                    ax.fill_between(x, mean - std, mean + std,
                                    color=color, alpha=0.20, linewidth=0)

            if mut_system in mean_by_system:
                x_curve, y_curve = mean_by_system[mut_system]
                y_map = {int(x): float(y) for x, y in zip(x_curve, y_curve)}
                annotate_single_mutation(ax, mut_name, mut_info, chain_id,
                                         res_to_pos, y_map)

            ax = plt.gca()
            for spine in ax.spines.values():
                spine.set_linewidth(2)
                spine.set_color('black')

            ax.set_xlabel("Residue index", fontsize=20, weight='bold')
            ax.set_ylabel("RMSF (Å)", fontsize=20, weight='bold')
            ax.set_title(f"WT vs {display_label} — {CHAIN_TO_HISTONE[chain_id]} (chain {chain_id})",
                         fontsize=22, pad=10, weight='bold')
            ax.grid(True, alpha=0.25)
            ax.legend(frameon=False, ncol=2, fontsize=20, loc="upper center")

            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(18)

            fig.tight_layout()
            savefig = f"RMSF_WT_vs_{mut_system}_{histone}_chain_{chain_id}_seaborn.png"
            out = BASE_DIR / savefig
            fig.savefig(out, bbox_inches='tight', dpi=DPI)
            print(f"[OK] Saved: {out}")
            plt.show()


if __name__ == "__main__":
    main()