import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
from matplotlib.patches import Patch

# =============================================================
# Global font and size settings
# =============================================================
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['axes.titlesize'] = 22
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['legend.fontsize'] = 18

# =============================================================
# System definitions, SHL pairs, and colors
# =============================================================
systems = [
    ("WT_rw", "WT", "WT"),
    ("R29Q_rw", "H2AR29Q", "H2A"),
    ("E56K_rw", "H2AE56K", "H2A"),
    ("K74N_rw", "H2AK74N", "H2A"),
    ("F70L_rw", "H2BF70L", "H2B"),
    ("E50K_rw", "H3E50K", "H3"),
    ("E73K_rw", "H3E73K", "H3"),
    ("E105K_rw", "H3E105K", "H3"),
]

# Symmetric SHL pairs to be plotted together
shl_pairs = [("2_6", "6_2"), ("3_5", "5_3")]
pair_suffix = ["pair1", "pair2"]          # used for output file naming
shl_labels_all = ["2_6", "6_2", "3_5", "5_3"]

runs = ["run1", "run2", "run3"]

# Color scheme
WT_COLOR = "#86A8E7"
MUT_COLOR = "#D95F4C"

fig_dir = "symmetry_boxplot_figures"
os.makedirs(fig_dir, exist_ok=True)

mutants = [sys for sys in systems if sys[0] != "WT_rw"]

# =============================================================
# Helper: load data for a given system and SHL label
# =============================================================
def load_data_for_system(system_folder, shl_label):
    system_short = system_folder.replace("_rw", "")
    all_vals = []
    for run in runs:
        filepath = os.path.join(system_folder, "distance_results",
                                f"{system_short}_{run}_SHL_{shl_label}_gap.dat")
        if os.path.exists(filepath):
            data = np.loadtxt(filepath, comments='#')
            if data.ndim == 1:
                vals = np.array([data[1]])
            else:
                vals = data[:, 1]
            all_vals.extend(vals)
        else:
            print(f"Warning: file not found: {filepath}")
    return np.array(all_vals)

# =============================================================
# Generate boxplots for each mutant
# =============================================================
for mutant_sys in mutants:
    mut_folder, mut_label, mut_group = mutant_sys
    wt_folder = "WT_rw"

    # Load data for all SHL labels
    all_shl_vals = {}
    skip = False
    for shl in shl_labels_all:
        wt_vals = load_data_for_system(wt_folder, shl)
        mut_vals = load_data_for_system(mut_folder, shl)
        if len(wt_vals) == 0 or len(mut_vals) == 0:
            print(f"Warning: missing data for {mut_label} at SHL {shl}, skipping mutant.")
            skip = True
            break
        all_shl_vals[shl] = {'WT': wt_vals, 'Mut': mut_vals}
    if skip:
        continue

    # Generate a separate figure for each symmetric pair
    for pair_idx, (shl_a, shl_b) in enumerate(shl_pairs):
        fig, ax = plt.subplots(figsize=(7, 5))

        wt_a = all_shl_vals[shl_a]['WT']
        wt_b = all_shl_vals[shl_b]['WT']
        mut_a = all_shl_vals[shl_a]['Mut']
        mut_b = all_shl_vals[shl_b]['Mut']

        box_data = [wt_a, wt_b, mut_a, mut_b]
        positions = [0, 1, 2, 3]
        colors = [WT_COLOR, WT_COLOR, MUT_COLOR, MUT_COLOR]

        # Box plots
        bp = ax.boxplot(box_data, positions=positions, widths=0.4,
                        patch_artist=True, showmeans=False, showfliers=False,
                        medianprops={'color': 'black', 'linewidth': 2.5,
                                     'solid_capstyle': 'butt'},
                        boxprops={'linewidth': 1.5, 'edgecolor': 'black'},
                        whiskerprops={'linewidth': 1.2, 'color': 'black'},
                        capprops={'linewidth': 1.2, 'color': 'black'})
        for patch, col in zip(bp['boxes'], colors):
            patch.set_facecolor(col)
            patch.set_edgecolor('black')

        # Scatter (jittered)
        jitter = 0.08
        for i, data in enumerate(box_data):
            x_jitter = np.random.uniform(-jitter, jitter, size=len(data))
            ax.scatter(positions[i] + x_jitter, data,
                       alpha=0.3, s=20, edgecolors='black',
                       c=colors[i], zorder=4)

        # X‑axis labels: add minus signs
        label_a = f"-{shl_a.replace('_','/+')}"
        label_b = f"-{shl_b.replace('_','/+')}"
        ax.set_xticks(positions)
        ax.set_xticklabels([f'WT\n{label_a}', f'WT\n{label_b}',
                            f'{mut_label}\n{label_a}', f'{mut_label}\n{label_b}'],
                           fontsize=18)
        ax.set_ylabel('Gapping Distance (Å)', fontweight='bold')

        # Thicken all four spines
        for spine in ax.spines.values():
            spine.set_linewidth(2)
            spine.set_color('black')

        ax.yaxis.grid(True, linestyle='--', alpha=0.4, color='gray')
        ax.set_axisbelow(True)
        ax.tick_params(direction='out', length=4, width=1.2)

        # Adjust y‑axis limits
        all_vals = np.concatenate(box_data)
        y_min, y_max = np.min(all_vals), np.max(all_vals)
        y_range = y_max - y_min
        if y_range == 0:
            y_range = 1.0
        y_max_plot = y_max + 0.25 * y_range
        y_min_plot = y_min - 0.05 * y_range
        ax.set_ylim(y_min_plot, y_max_plot)

        from matplotlib.ticker import MaxNLocator
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        # ---- Significance markers within each group ----
        # WT group: compare shl_a vs shl_b
        p_wt = stats.mannwhitneyu(wt_a, wt_b, alternative='two-sided').pvalue
        x1_wt, x2_wt = positions[0], positions[1]
        y_top_wt = max(np.max(wt_a), np.max(wt_b))
        y_line_wt = y_top_wt + 0.03 * y_range
        ax.plot([x1_wt, x2_wt], [y_line_wt, y_line_wt], color='k', linewidth=1.2)
        ax.plot([x1_wt, x1_wt], [y_line_wt, y_line_wt - 0.02 * y_range], color='k', linewidth=1.2)
        ax.plot([x2_wt, x2_wt], [y_line_wt, y_line_wt - 0.02 * y_range], color='k', linewidth=1.2)
        star_wt = '***' if p_wt < 0.001 else '**' if p_wt < 0.01 else '*' if p_wt < 0.05 else 'ns'
        ax.text((x1_wt + x2_wt)/2, y_line_wt + 0.02 * y_range, star_wt,
                ha='center', va='bottom', fontsize=18, fontweight='bold')

        # Mutant group: compare shl_a vs shl_b
        p_mut = stats.mannwhitneyu(mut_a, mut_b, alternative='two-sided').pvalue
        x1_mut, x2_mut = positions[2], positions[3]
        y_top_mut = max(np.max(mut_a), np.max(mut_b))
        y_line_mut = y_top_mut + 0.03 * y_range
        ax.plot([x1_mut, x2_mut], [y_line_mut, y_line_mut], color='k', linewidth=1.2)
        ax.plot([x1_mut, x1_mut], [y_line_mut, y_line_mut - 0.02 * y_range], color='k', linewidth=1.2)
        ax.plot([x2_mut, x2_mut], [y_line_mut, y_line_mut - 0.02 * y_range], color='k', linewidth=1.2)
        star_mut = '***' if p_mut < 0.001 else '**' if p_mut < 0.01 else '*' if p_mut < 0.05 else 'ns'
        ax.text((x1_mut + x2_mut)/2, y_line_mut + 0.02 * y_range, star_mut,
                ha='center', va='bottom', fontsize=16, fontweight='bold')

        # ---- Legend (top, horizontal) ----
        legend_elements = [
            Patch(facecolor=WT_COLOR, edgecolor='black', label='WT'),
            Patch(facecolor=MUT_COLOR, edgecolor='black', label=mut_label)
        ]
        ax.legend(handles=legend_elements, loc='upper center',
                  bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False)

        plt.subplots_adjust(top=0.9)   # leave room for the legend

        # Save
        outname = f"symmetry_boxes_{mut_label}_vs_WT_{pair_suffix[pair_idx]}.png"
        outpath = os.path.join(fig_dir, outname)
        fig.savefig(outpath, dpi=600, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {outpath}")

print(f"\nAll done. Figures saved to {fig_dir}/")