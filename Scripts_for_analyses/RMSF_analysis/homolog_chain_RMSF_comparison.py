#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compare RMSF profiles between two homologous chains (e.g., chain A vs chain E for H3)
within the same nucleosome system. Generates overlay plots and delta tables.
"""

from pathlib import Path
import numpy as np
import pandas as pd

import seaborn as sns
sns.set_theme(style="whitegrid", context="paper")

import matplotlib.pyplot as plt

# ===== GLOBAL FONT AND FRAME SETTINGS =====
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.edgecolor'] = 'black'

# ========== PATHS AND BASIC CONFIGURATION ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"

# Mapping histone types to their two chain IDs (homologous copies)
HISTONE_PAIR = {
    "H3":  ("A", "E"),
    "H4":  ("B", "F"),
    "H2A": ("C", "G"),
    "H2B": ("D", "H"),
}

CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

RMSF_FILE_TEMPLATE = "rmsf-{histone}_chain_{chain}.dat"

PLOT_STD_BAND = True
DPI = 300
DIFF_THRESHOLD = 0.5
OUTPUT_CSV = True

# Which systems and histones to analyze
SYSTEMS_TO_ANALYZE = {
    "WT_rw":   ["H2A", "H2B", "H3"],          # H4 excluded
    "R29Q_rw": ["H2A"],
    "E56K_rw": ["H2A"],
    "K74N_rw": ["H2A"],
    "F70L_rw": ["H2B"],
    "E50K_rw": ["H3"],
    "E73K_rw": ["H3"],
    "E105K_rw": ["H3"],
}

MUTATION_MAP = {
    "R29Q":  {"histone": "H2A", "chain": "C", "x_resnum": 29},
    "E56K":  {"histone": "H2A", "chain": "C", "x_resnum": 56},
    "K74N":  {"histone": "H2A", "chain": "C", "x_resnum": 74},
    "F70L":  {"histone": "H2B", "chain": "D", "x_resnum": 70},
    "E50K":  {"histone": "H3",  "chain": "A", "x_resnum": 50},
    "E73K":  {"histone": "H3",  "chain": "A", "x_resnum": 73},
    "E105K": {"histone": "H3",  "chain": "A", "x_resnum": 105},
}


# -------------------------------------------------------------
def get_system_label(system_name: str) -> str:
    """
    Convert system directory name to a display label.
    WT_rw      -> WT
    R29Q_rw    -> H2AR29Q (histone prefix from MUTATION_MAP)
    Others     -> remove _rw suffix.
    """
    if system_name == "WT_rw":
        return "WT"

    mut_key = system_name.replace("_rw", "")
    if mut_key in MUTATION_MAP:
        histone = MUTATION_MAP[mut_key]["histone"]
        return f"{histone}{mut_key}"

    return mut_key


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
    df = pd.read_csv(
        path, comment="#", delim_whitespace=True, header=None,
        names=["Res", "RMSF"], dtype={"Res": float, "RMSF": float},
    )
    if df.empty:
        raise ValueError(f"{path} parsed as empty dataframe.")
    df["Res"] = df["Res"].round().astype(int)
    return df


def load_system_chain_runs(system: str, histone: str, chain: str,
                           res_to_pos: pd.DataFrame, start_index: int):
    per_run = []
    for run in RUNS:
        f = BASE_DIR / system / run / RMSF_FILE_TEMPLATE.format(histone=histone, chain=chain)
        if not f.exists():
            print(f"[WARN] Missing file: {f}, skipping this run")
            continue
        df = read_cpptraj_rmsf(f)
        merged = df.merge(res_to_pos, on="Res", how="inner")
        if merged.empty:
            raise ValueError(
                f"No residues matched between {f} and sequence.csv for chain={chain}. "
                f"Maybe RMSF uses #Orig instead of #Res?"
            )
        merged["X"] = start_index + merged["Pos"] - 1
        merged = merged.sort_values("X")
        per_run.append((merged["X"].to_numpy(dtype=int), merged["RMSF"].to_numpy(dtype=float)))
    if not per_run:
        raise RuntimeError(f"No valid runs found for {system} {histone} chain {chain}")
    return per_run


def align_and_stack(series_list):
    sets = [set(x.tolist()) for x, _ in series_list]
    common = sorted(set.intersection(*sets))
    if not common:
        raise ValueError("No common X positions across runs (different residue sets).")
    idx = {v: i for i, v in enumerate(common)}
    Y = np.full((len(series_list), len(common)), np.nan, dtype=float)
    for i, (x, y) in enumerate(series_list):
        for xv, yv in zip(x, y):
            j = idx.get(int(xv))
            if j is not None:
                Y[i, j] = float(yv)
    valid = ~np.isnan(Y).any(axis=0)
    common_x = np.array(common)[valid]
    Y = Y[:, valid]
    return common_x, Y


def summarize_runs(per_run):
    common_x, Y = align_and_stack(per_run)
    mean = Y.mean(axis=0)
    std = Y.std(axis=0, ddof=1) if Y.shape[0] > 1 else np.zeros_like(mean)
    return common_x, mean, std


def mark_mutation_on_chain(ax, mutation: dict, chain_id: str, res_to_pos: pd.DataFrame,
                           x_curve, y_curve):
    if mutation["chain"] != chain_id:
        return
    x_curve = np.asarray(x_curve, dtype=int)
    y_curve = np.asarray(y_curve, dtype=float)
    y_map = {int(x): float(y) for x, y in zip(x_curve, y_curve)}
    chain_start = CHAIN_START[chain_id]

    x_mut = int(mutation["x_resnum"])
    pos_mut = x_mut - chain_start + 1
    if pos_mut < 1 or not (res_to_pos["Pos"] == pos_mut).any():
        print(f"[WARN] Mutation (X={x_mut}) not found in chain {chain_id} sequence")
        return
    if x_mut not in y_map:
        print(f"[WARN] Mutation X={x_mut} not found on RMSF curve for chain {chain_id}")
        return
    y_mut = y_map[x_mut]
    ax.plot(x_mut, y_mut, marker="o", markersize=7.5,
            markerfacecolor="crimson", markeredgecolor="white",
            markeredgewidth=1.2, linestyle="None", zorder=10)


def compute_chain_mean_for_system(system: str, chain_id: str, histone: str,
                                  res_to_pos: pd.DataFrame):
    per_run = load_system_chain_runs(
        system=system, histone=histone, chain=chain_id,
        res_to_pos=res_to_pos, start_index=CHAIN_START[chain_id],
    )
    return summarize_runs(per_run)


def diff_between_two_chains(chain1_id, chain2_id, histone,
                            res_to_pos1, res_to_pos2, system):
    x1, mean1, _ = compute_chain_mean_for_system(system, chain1_id, histone, res_to_pos1)
    x2, mean2, _ = compute_chain_mean_for_system(system, chain2_id, histone, res_to_pos2)

    common_x = np.intersect1d(x1, x2)
    if common_x.size == 0:
        raise ValueError(f"No common X between chains {chain1_id} and {chain2_id}")

    idx1 = {int(x): i for i, x in enumerate(x1)}
    idx2 = {int(x): i for i, x in enumerate(x2)}
    rows = []
    for x in common_x:
        i1 = idx1[int(x)]
        i2 = idx2[int(x)]
        m1 = float(mean1[i1])
        m2 = float(mean2[i2])
        diff = m2 - m1
        rows.append({
            "Histone": histone, "Chain1": chain1_id, "Chain2": chain2_id,
            "X": int(x), f"{chain1_id}_mean": m1, f"{chain2_id}_mean": m2,
            "diff": diff, "abs_diff": abs(diff),
        })
    return pd.DataFrame(rows).sort_values("X").reset_index(drop=True)


def print_diff_hits(df, threshold):
    hits = df[df["abs_diff"] >= float(threshold)].copy()
    if hits.empty:
        print(f"  No sites with |diff| >= {threshold}")
        return
    hits = hits.sort_values("abs_diff", ascending=False)
    print(f"  Sites with |Chain2 - Chain1| >= {threshold}:")
    for _, r in hits.iterrows():
        print(f"    X={int(r['X'])}  diff={r['diff']:.3f}  (abs={r['abs_diff']:.3f})")


def main():
    seq_path = BASE_DIR / SEQUENCE_CSV
    if not seq_path.exists():
        raise FileNotFoundError(f"Cannot find: {seq_path}")
    seq_df = read_sequence_csv(seq_path)

    all_chains = set()
    for ch1, ch2 in HISTONE_PAIR.values():
        all_chains.update([ch1, ch2])
    chain_res_map = {chain: build_chain_res_to_pos(seq_df, chain) for chain in all_chains}

    pal = sns.color_palette("deep", 2)
    color_chain1 = pal[0]
    color_chain2 = pal[1]

    for system_name, histones in SYSTEMS_TO_ANALYZE.items():
        sys_dir = BASE_DIR / system_name
        if not sys_dir.exists():
            print(f"[ERROR] System directory not found: {sys_dir}, skipping")
            continue

        title_label = get_system_label(system_name)

        for histone in histones:
            if histone not in HISTONE_PAIR:
                print(f"[WARN] {histone} not in HISTONE_PAIR, skipping")
                continue
            chain1, chain2 = HISTONE_PAIR[histone]
            print(f"\n=== {system_name} | {histone}: {chain1} vs {chain2} ===")

            x1, mean1, std1 = compute_chain_mean_for_system(
                system_name, chain1, histone, chain_res_map[chain1])
            x2, mean2, std2 = compute_chain_mean_for_system(
                system_name, chain2, histone, chain_res_map[chain2])

            fig, ax = plt.subplots(figsize=(8, 4))

            ax.plot(x1, mean1, color=color_chain1, linestyle="-", linewidth=2.0,
                    label=f"Chain {chain1} (n={len(RUNS)})")
            if PLOT_STD_BAND and len(RUNS) > 1:
                ax.fill_between(x1, mean1 - std1, mean1 + std1,
                                color=color_chain1, alpha=0.20, linewidth=0)

            ax.plot(x2, mean2, color=color_chain2, linestyle="--", linewidth=2.0,
                    label=f"Chain {chain2} (n={len(RUNS)})")
            if PLOT_STD_BAND and len(RUNS) > 1:
                ax.fill_between(x2, mean2 - std2, mean2 + std2,
                                color=color_chain2, alpha=0.20, linewidth=0)

            # Mark mutation if applicable
            mut_key = system_name.replace("_rw", "")
            if mut_key in MUTATION_MAP:
                mut_info = MUTATION_MAP[mut_key]
                mark_mutation_on_chain(ax, mut_info, chain1, chain_res_map[chain1], x1, mean1)
                mark_mutation_on_chain(ax, mut_info, chain2, chain_res_map[chain2], x2, mean2)

            ax = plt.gca()
            for spine in ax.spines.values():
                spine.set_linewidth(2)
                spine.set_color('black')

            ax.set_title(f"{title_label} | {chain1} vs {chain2}",
                         fontsize=20, pad=10, weight='bold')
            ax.set_xlabel("Residue index", fontsize=18, weight='bold')
            ax.set_ylabel("RMSF (Å)", fontsize=18, weight='bold')

            ax.grid(True, alpha=0.25)
            ax.legend(frameon=False, fontsize=18, loc="upper center", ncol=2)

            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontsize(18)

            fig.tight_layout()

            label = get_system_label(system_name)
            if label == "WT":
                name_part = f"{histone}_WT"
            else:
                name_part = label

            out_png = BASE_DIR / f"RMSF_{name_part}_{chain1}vs{chain2}.png"
            fig.savefig(out_png, dpi=DPI)
            print(f"[OK] Figure saved: {out_png}")
            plt.close(fig)

            df_diff = diff_between_two_chains(chain1, chain2, histone,
                                              chain_res_map[chain1], chain_res_map[chain2],
                                              system=system_name)
            if OUTPUT_CSV:
                csv_out = BASE_DIR / f"RMSF_delta_{name_part}_{chain1}-{chain2}.csv"
                df_diff.to_csv(csv_out, index=False)

            print_diff_hits(df_diff, DIFF_THRESHOLD)

    print("\nAll done.")


if __name__ == "__main__":
    main()