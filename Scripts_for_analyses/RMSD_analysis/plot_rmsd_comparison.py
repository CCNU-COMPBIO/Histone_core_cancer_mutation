#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ===== User parameters =====
# Mutants and their full display names
MUTANTS = {
    "R29Q": "H2AR29Q",
    "E56K": "H2AE56K",
    "K74N": "H2AK74N",
    "F70L": "H2BF70L",
    "E50K": "H3E50K",
    "E73K": "H3E73K",
    "E105K": "H3E105K",
}

# Chain information for each mutant (used for chain‑specific RMSD)
MUTATION_MAP = {
    "R29Q": {"histone": "H2A", "chain": "C", "x_resnum": 29},
    "E56K": {"histone": "H2A", "chain": "C", "x_resnum": 56},
    "K74N": {"histone": "H2A", "chain": "C", "x_resnum": 74},
    "F70L": {"histone": "H2B", "chain": "D", "x_resnum": 70},
    "E50K": {"histone": "H3", "chain": "A", "x_resnum": 50},
    "E73K": {"histone": "H3", "chain": "A", "x_resnum": 73},
    "E105K": {"histone": "H3", "chain": "A", "x_resnum": 105},
}

RUNS = ["run1", "run2", "run3"]

# Time conversion: ns per frame. Set to None to plot frame numbers.
DT_NS = 0.02  # 0.02 ns/frame
MAX_NS = 500  # Only plot the first 500 ns

# Plotting options
PLOT_SEM_BAND = True          # Shade mean ± standard error
SHOW_INDIVIDUAL_RUNS = False  # Show individual run traces as faint lines
DPI = 600

# Global font settings (Arial)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False

TITLE_FONTSIZE = 22
AXIS_LABEL_FONTSIZE = 20
TICK_FONTSIZE = 18
LEGEND_FONTSIZE = 20

# Fixed colors: WT blue, mutant red
COLOR_WT = "#1f77b4"
COLOR_MUT = "#d62728"


# ===== Helper functions =====
def read_cpptraj_rmsd(path: Path) -> pd.DataFrame:
    """
    Read a cpptraj‑generated RMSD file.
    Expected format:
        #Frame    RMSD
        1         0.0000
    """
    df = pd.read_csv(
        path,
        comment="#",
        delim_whitespace=True,
        header=None,
        names=["Frame", "RMSD"],
        dtype={"Frame": int, "RMSD": float},
    )
    if df.empty:
        raise ValueError(f"{path} parsed as empty dataframe.")
    return df


def align_and_stack(series_list):
    """
    Align multiple (x, y) series to a common x‑axis (intersection).
    Returns the common x‑coordinates and a matrix Y of shape [n_runs, n_points].
    """
    x_sets = [set(x.tolist()) for x, _ in series_list]
    common_x = sorted(set.intersection(*x_sets))
    if len(common_x) == 0:
        raise ValueError("No common x-axis points across runs.")

    x_index = {v: i for i, v in enumerate(common_x)}
    Y = np.full((len(series_list), len(common_x)), np.nan, dtype=float)

    for i, (x, y) in enumerate(series_list):
        for xv, yv in zip(x, y):
            j = x_index.get(int(xv))
            if j is not None:
                Y[i, j] = float(yv)

    # Remove columns with any NaN (frames missing in some runs)
    valid_cols = ~np.isnan(Y).any(axis=0)
    common_x = np.array(common_x)[valid_cols]
    Y = Y[:, valid_cols]
    return common_x, Y


def plot_comparison(ax, systems, runs, base_path, display_names, filename):
    """
    Plot RMSD comparison between two systems on the given Axes object.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to draw on.
    systems : list of str
        Two system names, e.g. ["WT_rw", "R29Q_rw"].
    runs : list of str
        List of run identifiers, e.g. ["run1", "run2", "run3"].
    base_path : Path
        Base directory where the system/run folders are located.
    display_names : dict
        Mapping from system name to display label, e.g. {"WT_rw":"WT", "R29Q_rw":"H2AR29Q"}.
    filename : str
        Name of the RMSD data file to read (e.g. "rmsd_histone.dat" or "rmsd_chain_C.dat").

    Returns
    -------
    bool
        True if both systems had valid data, False otherwise.
    """
    colors = {systems[0]: COLOR_WT, systems[1]: COLOR_MUT}
    max_frame = int(MAX_NS / DT_NS) if DT_NS is not None else None

    for system in systems:
        per_run = []
        for run in runs:
            f = base_path / system / run / filename
            try:
                if not f.exists():
                    raise FileNotFoundError(f"File not found: {f}")
                df = read_cpptraj_rmsd(f)

                # Time truncation
                if DT_NS is not None:
                    df = df[df["Frame"] * DT_NS <= MAX_NS]

                # Downsample to 1 ns spacing (if DT_NS is defined)
                if DT_NS is not None:
                    sampling_step = int(round(1.0 / DT_NS))  # e.g. 50 for 0.02 ns/frame
                    df = df[df["Frame"] % sampling_step == 0]
                    if df.empty:
                        print(f"Warning: {f} has no data after downsampling to 1 ns, skipping.")
                        continue

                x = df["Frame"].to_numpy()
                y = df["RMSD"].to_numpy()
                per_run.append((x, y))

                # Plot individual runs as faint lines (if enabled)
                if SHOW_INDIVIDUAL_RUNS:
                    ax.plot(
                        x if DT_NS is None else x * DT_NS,
                        y,
                        color=colors[system],
                        alpha=0.25,
                        linewidth=1.0,
                        label=None,
                    )
            except Exception as e:
                print(f"Warning: Failed to read {f}: {e}. Skipping this run.")
                continue

        # If no valid runs for this system, skip entirely
        if not per_run:
            print(f"Warning: No valid data for system {system} in file {filename}, skipping this mutant.")
            return False

        # Align runs and compute mean and standard deviation
        common_x, Y = align_and_stack(per_run)
        mean = Y.mean(axis=0)
        std = Y.std(axis=0, ddof=1) if Y.shape[0] > 1 else np.zeros_like(mean)

        # Standard error of the mean
        n_runs = Y.shape[0]
        sem = std / np.sqrt(n_runs)

        plot_x = common_x if DT_NS is None else common_x * DT_NS
        display_name = display_names.get(system, system.replace("_rw", ""))

        ax.plot(
            plot_x,
            mean,
            color=colors[system],
            linewidth=2.5,
            label=f"{display_name} (n={len(per_run)})",
        )

        if PLOT_SEM_BAND and len(per_run) > 1:
            ax.fill_between(
                plot_x,
                mean - sem,
                mean + sem,
                color=colors[system],
                alpha=0.4,
                linewidth=0,
                label=None,
            )

    return True


def main():
    base = Path(__file__).resolve().parent

    for mutant, full_name in MUTANTS.items():
        systems = ["WT_rw", f"{mutant}_rw"]
        display_names = {
            "WT_rw": "WT",
            f"{mutant}_rw": full_name,
        }
        print(f"\nProcessing: WT vs {full_name}")

        # ----- 1. Histone‑wide RMSD (rmsd_histone.dat) -----
        fig1, ax1 = plt.subplots(figsize=(8, 4))
        success1 = plot_comparison(ax1, systems, RUNS, base, display_names, "rmsd_histone.dat")

        if success1:
            for spine in ax1.spines.values():
                spine.set_linewidth(2)
                spine.set_color('black')
            ax1.set_xlabel("Time (ns)" if DT_NS is not None else "Frame",
                           fontsize=AXIS_LABEL_FONTSIZE, weight='bold')
            ax1.set_ylabel("RMSD (Å)", fontsize=AXIS_LABEL_FONTSIZE, weight='bold')
            ax1.set_title(f"RMSD comparison: WT vs {full_name}",
                          fontsize=TITLE_FONTSIZE, pad=10, weight='bold')
            ax1.grid(True, alpha=0.3)
            ax1.legend(frameon=False, fontsize=LEGEND_FONTSIZE, loc='lower right')
            ax1.tick_params(axis='both', labelsize=TICK_FONTSIZE)

            fig1.tight_layout()
            out_name1 = f"rmsd_histone_WT-{mutant}.png"
            fig1.savefig(out_name1, bbox_inches='tight', dpi=DPI)
            plt.close(fig1)
            print(f"  -> Histone RMSD saved: {out_name1}")
        else:
            plt.close(fig1)
            print(f"  -> Skipped histone RMSD for {full_name} (no data)")

        # ----- 2. Chain‑specific RMSD (only if chain info is defined) -----
        if mutant in MUTATION_MAP:
            chain = MUTATION_MAP[mutant]["chain"]
            filename_chain = f"rmsd_chain_{chain}.dat"
            fig2, ax2 = plt.subplots(figsize=(8, 4))
            success2 = plot_comparison(ax2, systems, RUNS, base, display_names, filename_chain)

            if success2:
                for spine in ax2.spines.values():
                    spine.set_linewidth(2)
                    spine.set_color('black')
                ax2.set_xlabel("Time (ns)" if DT_NS is not None else "Frame",
                               fontsize=AXIS_LABEL_FONTSIZE, weight='bold')
                ax2.set_ylabel("RMSD (Å)", fontsize=AXIS_LABEL_FONTSIZE, weight='bold')
                ax2.set_title(f"Chain {chain} RMSD comparison: WT vs {full_name}",
                              fontsize=TITLE_FONTSIZE, pad=10, weight='bold')
                ax2.grid(True, alpha=0.3)
                ax2.legend(frameon=False, fontsize=LEGEND_FONTSIZE, loc='lower right')
                ax2.tick_params(axis='both', labelsize=TICK_FONTSIZE)

                fig2.tight_layout()
                out_name2 = f"rmsd_chain_{chain}_WT-{mutant}.png"
                fig2.savefig(out_name2, bbox_inches='tight', dpi=DPI)
                plt.close(fig2)
                print(f"  -> Chain {chain} RMSD saved: {out_name2}")
            else:
                plt.close(fig2)
                print(f"  -> Skipped chain {chain} RMSD for {full_name} (no data)")
        else:
            print(f"  -> No chain info for {full_name}, chain RMSD skipped.")

    print("\nAll figures generated successfully.")


if __name__ == "__main__":
    main()
