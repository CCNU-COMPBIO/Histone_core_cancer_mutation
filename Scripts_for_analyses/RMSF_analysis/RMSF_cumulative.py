#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S25 - time convergence of RMSF.
Changes: line width 1.2; wider color-scale span; 0-100 on top; SEM band removed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid", context="paper")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
WINDOWS = [100, 200, 300, 400, 500]
DPI = 600
FIGSIZE = (6.0, 4.0)

# Sequential color scale: light -> dark maps to 100 -> 500 ns
WIN_COLORS = ["#9ECAE1", "#6BAED6", "#2171B5", "#08519C", "#08306B"]

PANELS = [
    ("WT_rw",        "A", "WT (H3)",   "WT_H3"),
    ("WT_rw",        "C", "WT (H2A)",  "WT_H2A"),
    ("WT_rw",        "D", "WT (H2B)",  "WT_H2B"),
    ("H3_E50K_rw",   "A", "H3E50K",    "H3E50K"),
    ("H3_E73K_rw",   "A", "H3E73K",    "H3E73K"),
    ("H3_E105K_rw",  "A", "H3E105K",   "H3E105K"),
    ("H2A_R29Q_rw",  "C", "H2AR29Q",   "H2AR29Q"),
    ("H2A_E56K_rw",  "C", "H2AE56K",   "H2AE56K"),
    ("H2A_K74N_rw",  "C", "H2AK74N",   "H2AK74N"),
    ("H2B_F70L_rw",  "D", "H2BF70L",   "H2BF70L"),
    ("H2B_E105K_rw", "D", "H2BE105K",  "H2BE105K"),
]

CHAIN_START = {"A": 38, "E": 37, "B": 25, "F": 21,
               "C": 14, "G": 15, "D": 30, "H": 33}
CHAIN_TO_HISTONE = {"A": "H3", "B": "H4", "C": "H2A", "D": "H2B",
                    "E": "H3", "F": "H4", "G": "H2A", "H": "H2B"}


def read_sequence_csv(path):
    for enc in ["utf-8", "utf-8-sig", "gbk", "cp936", "latin1"]:
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError:
            df = None
    ren = {}
    if "#Res" in df.columns: ren["#Res"] = "Res"
    if "chain" in df.columns: ren["chain"] = "Chain"
    df = df.rename(columns=ren)
    df["Res"] = pd.to_numeric(df["Res"], errors="coerce")
    df = df.dropna(subset=["Res"]).copy()
    df["Res"] = df["Res"].astype(int)
    df["Chain"] = df["Chain"].astype(str).str.strip()
    return df


def build_chain_res_to_pos(seq_df, chain_id):
    sub = seq_df[seq_df["Chain"] == chain_id].sort_values("Res").reset_index(drop=True)
    sub["Pos"] = np.arange(1, len(sub) + 1)
    return sub[["Res", "Pos"]]


def read_cpptraj_rmsf(path):
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None,
                     names=["Res", "RMSF"])
    df["Res"] = df["Res"].round().astype(int)
    return df


def compute_window(system, hist, chain, res_to_pos, start_index, T):
    per_run = []
    for run in RUNS:
        f = BASE_DIR / system / run / f"rmsf-{hist}_chain_{chain}_cum{T}ns.dat"
        if not f.exists():
            raise FileNotFoundError(f"missing: {f}")
        df = read_cpptraj_rmsf(f).merge(res_to_pos, on="Res", how="inner")
        df["X"] = start_index + df["Pos"] - 1
        df = df.sort_values("X")
        per_run.append((df["X"].to_numpy(int), df["RMSF"].to_numpy(float)))
    common = sorted(set(per_run[0][0]).intersection(
        *[set(x) for x, _ in per_run[1:]]))
    idx = {v: i for i, v in enumerate(common)}
    Y = np.full((len(per_run), len(common)), np.nan)
    for i, (x, y) in enumerate(per_run):
        for xv, yv in zip(x, y):
            if int(xv) in idx:
                Y[i, idx[int(xv)]] = yv
    mean = np.nanmean(Y, axis=0)
    return np.array(common), mean


def draw_panel(system, chain, title, fname, seq):
    hist = CHAIN_TO_HISTONE[chain]
    res_to_pos = build_chain_res_to_pos(seq, chain)
    start = CHAIN_START[chain]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    x_ref = None
    for i, T in enumerate(WINDOWS):
        x, mean = compute_window(system, hist, chain, res_to_pos, start, T)
        x_ref = x
        # Draw 0-100 on top
        z = 10 - i
        ax.plot(x, mean, color=WIN_COLORS[i], linewidth=1.4, zorder=z)

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color("black")
    ax.set_xlabel("Residue index", fontsize=20, weight="bold")
    ax.set_ylabel("RMSF (Å)", fontsize=20, weight="bold")
    ax.set_title(title, fontsize=22, pad=10, weight="bold")
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlim(x_ref.min(), x_ref.max())
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(18)

    fig.tight_layout()
    out = BASE_DIR / f"FigS25_{fname}.png"
    fig.savefig(out, bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print("  [OK]", out.name)


def draw_legend_panel():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axis("off")
    handles = [Line2D([0], [0], color=WIN_COLORS[i], lw=2.0,
                      label=f"0–{T} ns") for i, T in enumerate(WINDOWS)]
    ax.legend(handles=handles, loc="center", frameon=False,
              fontsize=22, handlelength=2.4, labelspacing=0.8)
    out = BASE_DIR / "FigS25_legend.png"
    fig.savefig(out, bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print("  [OK]", out.name)


def main():
    seq = read_sequence_csv(BASE_DIR / "sequence.csv")
    for system, chain, title, fname in PANELS:
        draw_panel(system, chain, title, fname, seq)
    draw_legend_panel()


if __name__ == "__main__":
    main()