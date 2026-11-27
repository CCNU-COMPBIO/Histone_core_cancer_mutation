#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. S24 - per-replicate RMSF: 11 panels + 1 standalone legend panel.
Changes: darker/thicker lines with distinct line styles; removed red dots;
legend shows only run 1/2/3 + mean.
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
DPI = 600
FIGSIZE = (6.0, 4.0)

# Per-run line style: three saturated colors, all solid
RUN_STYLE = [
    dict(color="#0072B2", ls="-", lw=1.8),   # run 1, blue
    dict(color="#D55E00", ls="-", lw=1.8),   # run 2, orange-red
    dict(color="#009E73", ls="-", lw=1.8),   # run 3, teal
]

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
    df = pd.read_csv(path, comment="#", delim_whitespace=True, header=None,
                     names=["Res", "RMSF"])
    df["Res"] = df["Res"].round().astype(int)
    return df


def compute_per_run(system, hist, chain, res_to_pos, start_index):
    per_run = []
    for run in RUNS:
        f = BASE_DIR / system / run / f"rmsf-{hist}_chain_{chain}.dat"
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
    sem = np.nanstd(Y, axis=0, ddof=1) / np.sqrt(Y.shape[0])
    return np.array(common), Y, mean, sem


def draw_panel(system, chain, title, fname, seq):
    hist = CHAIN_TO_HISTONE[chain]
    res_to_pos = build_chain_res_to_pos(seq, chain)
    x, Y, mean, sem = compute_per_run(system, hist, chain, res_to_pos,
                                      CHAIN_START[chain])

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # Semi-transparent SEM band (annotated in caption, not shown in legend)
    ax.fill_between(x, mean - sem, mean + sem, color="black",
                    alpha=0.12, linewidth=0, zorder=1)

    for i in range(Y.shape[0]):
        ax.plot(x, Y[i], zorder=2 + i, **RUN_STYLE[i])
    ax.plot(x, mean, color="black", linewidth=2.4, linestyle=(0, (5, 2)), zorder=6)

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color("black")

    ax.set_xlabel("Residue index", fontsize=20, weight="bold")
    ax.set_ylabel("RMSF (Å)", fontsize=20, weight="bold")
    ax.set_title(title, fontsize=22, pad=10, weight="bold")
    ax.grid(True, alpha=0.25)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    for lab in ax.get_xticklabels() + ax.get_yticklabels():
        lab.set_fontsize(18)

    fig.tight_layout()
    out = BASE_DIR / f"FigS24_{fname}.png"
    fig.savefig(out, bbox_inches="tight", dpi=DPI)
    plt.close(fig)
    print("  [OK]", out.name)


def draw_legend_panel():
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.axis("off")
    handles = [
        Line2D([0], [0], label="run 1", **RUN_STYLE[0]),
        Line2D([0], [0], label="run 2", **RUN_STYLE[1]),
        Line2D([0], [0], label="run 3", **RUN_STYLE[2]),
        Line2D([0], [0], color="black", lw=2.4,
               linestyle=(0, (5, 2)), label="mean"),
    ]
    ax.legend(handles=handles, loc="center", frameon=False,
              fontsize=22, handlelength=2.4, labelspacing=0.9)
    out = BASE_DIR / "FigS24_legend.png"
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