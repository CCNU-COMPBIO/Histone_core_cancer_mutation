from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Parameters
excel_path = "MMGBSA_result.xlsx"
sheet_name = "mmgbsa"

mut_side_chains = ["chainA-chainB", "chainC-chainD", "chainA-DNA", "chainC-DNA"]
OUT_MUT = "mmgbsa_mutated_side.png"
DPI = 600

CHAIN_GROUP_GAP = 2.2
BAR_WIDTH = 0.28
EDGE_LW = 1.0
CAPSIZE = 4

SYSTEM_LABEL_Y = -3.0
SYSTEM_LABEL_ROT = 90

SHOW_RUN_POINTS = True
POINT_SIZE = 35
POINT_ALPHA = 0.95
POINT_EDGE = "black"
POINT_LW = 1.0
POINT_JITTER = 0.06

DRAW_LEGEND = False
RUN_COLS = ["run1", "run2", "run3"]

# Font settings (Arial)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.unicode_minus"] = False

TITLE_FONTSIZE = 20
YLABEL_FONTSIZE = 20
XTICK_FONTSIZE = 20
YTICK_FONTSIZE = 20
SYSTEM_LABEL_FONTSIZE = 20
LEGEND_FONTSIZE = 20
LEGEND_TITLE_FONTSIZE = 20

# Read data and compute mean/SEM from runs
df = pd.read_excel(excel_path, sheet_name=sheet_name)

required_cols = {"chain", "system", *RUN_COLS}
missing = required_cols - set(df.columns)
if missing:
    raise ValueError(f"Excel missing columns: {missing}. Actual columns: {list(df.columns)}")

df = df.dropna(subset=["chain", "system"]).copy()
df["chain"] = df["chain"].astype(str).str.strip()
df["system"] = df["system"].astype(str).str.strip()

for c in RUN_COLS:
    df[c] = pd.to_numeric(df[c], errors="coerce")

runs_mat = df[RUN_COLS].to_numpy(dtype=float)
df["mean_calc"] = np.nanmean(runs_mat, axis=1)

n = np.sum(~np.isnan(runs_mat), axis=1)
std = np.nanstd(runs_mat, axis=1, ddof=1)
sem = std / np.sqrt(n)
sem[n < 2] = np.nan
df["SEM_calc"] = sem

df = df.dropna(subset=["mean_calc"])

# =========================
# Color scheme
# =========================
def chain_color_map(chains):
    cmap_primary = plt.get_cmap("Set2")
    cmap_fallback = plt.get_cmap("Pastel1")
    colors = {}
    for i, ch in enumerate(chains):
        if i < cmap_primary.N:
            colors[ch] = cmap_primary(i)
        else:
            colors[ch] = cmap_fallback(i % cmap_fallback.N)
    return colors

# Chain display name mapping (supports composite names)
def map_display_name(raw_name: str) -> str:
    """
    Convert chainA->H3, chainB->H4, chainC->H2A, chainD->H2B,
    also handles composite names like "chainA-chainB" -> "H3-H4".
    """
    mapping = {
        "chainA": "H3",
        "chainB": "H4",
        "chainC": "H2A",
        "chainD": "H2B",
    }
    result = raw_name
    for old, new in mapping.items():
        result = result.replace(old, new)
    return result

# Plotting function
def plot_side(data: pd.DataFrame, chains_order, title: str, out_png: str):
    data = data[data["chain"].isin(chains_order)].copy()
    if data.empty:
        raise ValueError(f"No data found for chains: {chains_order}")

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

        dup = sub.duplicated(subset=["chain", "system"], keep=False)
        if dup.any():
            dups = sub.loc[dup, ["chain", "system"]].drop_duplicates()
            raise ValueError(
                "Duplicate (chain, system) entries found; aggregation required.\n"
                f"Example duplicates:\n{dups.to_string(index=False)}"
            )

        systems = sorted(sub["system"].unique().tolist())
        n_sys = len(systems)
        offsets = (np.arange(n_sys) - (n_sys - 1) / 2.0) * BAR_WIDTH

        for si, sys in enumerate(systems):
            row = sub[sub["system"] == sys].iloc[0]
            mean = float(row["mean_calc"])
            sem = row["SEM_calc"]
            x = group_centers[gi] + offsets[si]

            yerr = None if pd.isna(sem) else sem
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

            # System label below bar
            ax.text(
                x, SYSTEM_LABEL_Y, sys,
                rotation=SYSTEM_LABEL_ROT,
                ha="center", va="top",
                fontsize=SYSTEM_LABEL_FONTSIZE,
                color="black",
                zorder=4
            )

    # Thicken spines
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    ax.set_title(title, fontsize=TITLE_FONTSIZE, weight='bold', pad=15)
    ax.set_ylabel("Binding free energy (kcal/mol)", fontsize=YLABEL_FONTSIZE, weight='bold', labelpad=10)

    # Replace chain names using map_display_name
    display_labels = [map_display_name(ch) for ch in chains_present]
    ax.set_xticks(group_centers)
    ax.set_xticklabels(display_labels, fontsize=XTICK_FONTSIZE, weight='bold')

    ax.tick_params(axis='y', labelsize=YTICK_FONTSIZE)
    ax.grid(axis="y", linestyle="--", alpha=0.35, zorder=0)

    y_min, y_max = ax.get_ylim()
    ax.set_ylim(y_min, 0)

    if DRAW_LEGEND:
        handles = [Patch(facecolor=colors[ch], edgecolor="black", label=ch) for ch in chains_present]
        ax.legend(
            handles=handles,
            title="chain",
            loc="lower left",
            bbox_to_anchor=(1.01, 0.0),
            frameon=False,
            fontsize=LEGEND_FONTSIZE,
            title_fontsize=LEGEND_TITLE_FONTSIZE
        )

    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
    plt.close(fig)

# Generate plot (mutated side only)
plot_side(df, mut_side_chains, "MMGBSA Free Energy (Mutated side)", OUT_MUT)
print(f"Done. Saved:\n- {OUT_MUT}")
