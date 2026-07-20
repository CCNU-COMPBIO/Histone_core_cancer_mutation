"""
Plot MMGBSA binding free energy from Excel data.
Only mutated side chains are plotted.
Individual run values are overlaid as scatter points.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---- Input / Output ----
EXCEL_PATH = "MMGBSA_result.xlsx"
SHEET_NAME = "mmgbsa"
MUT_CHAINS = ["chainA-chainB", "chainC-chainD", "chainA-DNA", "chainC-DNA"]
OUTPUT_PNG = "MMGBSA_binding_free_energy.png"
DPI = 300

# ---- Bar plot settings ----
CHAIN_GROUP_GAP = 2.2          # spacing between chain groups
BAR_WIDTH = 0.28
EDGE_LW = 1.0
CAPSIZE = 4

# ---- System label (below bars) ----
SYSTEM_LABEL_Y = -3.0
SYSTEM_LABEL_ROT = 90

# ---- Scatter overlay for individual runs ----
SHOW_RUN_POINTS = True
POINT_SIZE = 35
POINT_ALPHA = 0.95
POINT_EDGE = "black"
POINT_LW = 1.0
POINT_JITTER = 0.06

RUN_COLS = ["run1", "run2", "run3"]

# ---- Font settings (Arial) ----
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False

TITLE_FONTSIZE = 20
YLABEL_FONTSIZE = 18
XTICK_FONTSIZE = 18
YTICK_FONTSIZE = 18
SYSTEM_LABEL_FONTSIZE = 18

# ---- Read and preprocess data ----
df = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_NAME)
required_cols = {"chain", "system", *RUN_COLS}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing columns: {missing}. Actual: {list(df.columns)}")

df = df.dropna(subset=["chain", "system"]).copy()
df["chain"] = df["chain"].astype(str).str.strip()
df["system"] = df["system"].astype(str).str.strip()
for c in RUN_COLS:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Compute mean and SEM across runs
runs_mat = df[RUN_COLS].to_numpy(dtype=float)
df["mean_calc"] = np.nanmean(runs_mat, axis=1)
n = np.sum(~np.isnan(runs_mat), axis=1)
std = np.nanstd(runs_mat, axis=1, ddof=1)
sem = std / np.sqrt(n)
sem[n < 2] = np.nan
df["SEM_calc"] = sem
df = df.dropna(subset=["mean_calc"])

# ---- Color mapping for chains ----
def chain_color_map(chains):
    cmap_primary = plt.get_cmap("Set2")
    cmap_fallback = plt.get_cmap("Pastel1")
    colors = {}
    for i, ch in enumerate(chains):
        colors[ch] = cmap_primary(i) if i < cmap_primary.N else cmap_fallback(i % cmap_fallback.N)
    return colors

# ---- Main plotting function ----
def plot_mmgbsa(data, chains_order, title, out_png):
    data = data[data["chain"].isin(chains_order)].copy()
    if data.empty:
        raise ValueError(f"No data for chains: {chains_order}")

    data["chain"] = pd.Categorical(data["chain"], categories=chains_order, ordered=True)
    data = data.sort_values(["chain", "system"])

    chains_present = [c for c in chains_order if c in set(data["chain"].astype(str))]
    colors = chain_color_map(chains_present)
    group_centers = np.arange(len(chains_present)) * CHAIN_GROUP_GAP

    fig, ax = plt.subplots(figsize=(16, 5), constrained_layout=True)

    for gi, ch in enumerate(chains_present):
        sub = data[data["chain"] == ch].copy()
        if sub.empty:
            continue

        # Ensure no duplicate (chain, system) entries
        dup = sub.duplicated(subset=["chain", "system"], keep=False)
        if dup.any():
            dups = sub.loc[dup, ["chain", "system"]].drop_duplicates()
            raise ValueError(
                "Duplicate (chain, system) rows. Please aggregate first.\n"
                f"Examples:\n{dups.to_string(index=False)}"
            )

        systems = sorted(sub["system"].unique().tolist())
        n_sys = len(systems)
        offsets = (np.arange(n_sys) - (n_sys - 1) / 2.0) * BAR_WIDTH

        for si, sys in enumerate(systems):
            row = sub[sub["system"] == sys].iloc[0]
            mean = float(row["mean_calc"])
            sem_val = row["SEM_calc"]
            x = group_centers[gi] + offsets[si]

            # Bar with error bar
            yerr = None if pd.isna(sem_val) else sem_val
            ax.bar(
                x, mean,
                width=BAR_WIDTH,
                color=colors[ch],
                edgecolor="black",
                linewidth=EDGE_LW,
                yerr=yerr,
                capsize=CAPSIZE,
                zorder=1
            )

            # Scatter overlay for individual runs
            if SHOW_RUN_POINTS:
                run_vals = row[RUN_COLS].to_numpy(dtype=float)
                run_vals = run_vals[~np.isnan(run_vals)]
                if run_vals.size:
                    jit = np.linspace(-POINT_JITTER, POINT_JITTER, run_vals.size) if run_vals.size > 1 else np.array([0.0])
                    ax.scatter(
                        np.full(run_vals.size, x) + jit,
                        run_vals,
                        s=POINT_SIZE,
                        alpha=POINT_ALPHA,
                        facecolor="white",
                        edgecolor=POINT_EDGE,
                        linewidth=POINT_LW,
                        zorder=3
                    )

            # System name label below the bar
            ax.text(
                x, SYSTEM_LABEL_Y, sys,
                rotation=SYSTEM_LABEL_ROT,
                ha="center", va="top",
                fontsize=SYSTEM_LABEL_FONTSIZE,
                color="black",
                zorder=4
            )

    # Axis styling
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    ax.set_title(title, fontsize=TITLE_FONTSIZE, weight='bold', pad=12)
    ax.set_ylabel("Binding free energy (kcal/mol)", fontsize=YLABEL_FONTSIZE, weight='bold', labelpad=10)
    ax.set_xticks(group_centers)
    ax.set_xticklabels(chains_present, fontsize=XTICK_FONTSIZE, weight='bold')
    ax.tick_params(axis='y', labelsize=YTICK_FONTSIZE)
    ax.grid(axis="y", linestyle="--", alpha=0.35, zorder=0)

    # Keep bottom space for system labels
    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min, 0)

    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

# ---- Execute ----
plot_mmgbsa(
    df,
    MUT_CHAINS,
    "MMGBSA Binding Free Energy",
    OUTPUT_PNG
)

print(f"Done. Saved: {OUTPUT_PNG}")