import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ========== 1. Data loading ==========
csv_path = "SRMSF_PseudoSymmetry_Deviation_Trim5_PerRunDetail.csv"
df_raw = pd.read_csv(csv_path)

# Build labels: WT + mutant labels (e.g., H2AR29Q)
labels = ["WT"] + [f"{row['Histone_Subunit']}{row['Variant']}" for _, row in df_raw.iterrows()]

# Extract means and SEMs
means = [df_raw["SRMSF_WT_Mean_A"].iloc[0]] + df_raw["SRMSF_Variant_Mean_A"].tolist()
sems  = [df_raw["SRMSF_WT_SEM_A"].iloc[0]]  + df_raw["SRMSF_Variant_SEM_A"].tolist()

S_WT = means[0]

# ========== 2. Color assignment by histone subunit ==========
COLOR_MAP = {
    "H2A": "#D16BA5",
    "H2B": "#86A8E7",
    "H3":  "#ff9671",
    "WT":  "#808080",
}
bar_colors = [COLOR_MAP["WT"]]
for _, row in df_raw.iterrows():
    histone = row['Histone_Subunit']
    bar_colors.append(COLOR_MAP.get(histone, "#808080"))

# ========== 3. Global font settings ==========
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.unicode_minus': False,
    'font.size': 18,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
})

# ========== 4. Plotting (vertical bar chart with alpha=0.7) ==========
fig, ax = plt.subplots(figsize=(16, 5))

x_pos = np.arange(len(labels))

bars = ax.bar(
    x_pos, means, yerr=sems, width=0.55,
    color=bar_colors,
    edgecolor='black', linewidth=1.5,
    alpha=0.7,
    error_kw=dict(ecolor='#333333', elinewidth=1.5, capsize=4, capthick=1.5),
    zorder=2
)

# WT reference line (horizontal dashed)
ax.axhline(y=S_WT, color='#666666', linewidth=1.0, linestyle='--', zorder=1)

# ========== 5. Axes and labels ==========
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, rotation=45, ha='right', va='top')
ax.set_ylabel(r"$\mathbf{S_{RMSF}}$ (Å)", fontweight='bold', fontsize=20)

y_max = max(np.nanmax(means) + np.nanmax(sems), S_WT) * 1.15
ax.set_ylim(0, y_max)

for spine in ax.spines.values():
    spine.set_linewidth(2.0)
    spine.set_edgecolor('black')

# ========== 6. Legend (with alpha matching bars) ==========
legend_elements = [
    Patch(facecolor=COLOR_MAP["H2A"], edgecolor='black', label='H2A mutant', alpha=0.7),
    Patch(facecolor=COLOR_MAP["H2B"], edgecolor='black', label='H2B mutant', alpha=0.7),
    Patch(facecolor=COLOR_MAP["H3"], edgecolor='black', label='H3 mutant', alpha=0.7),
    Patch(facecolor=COLOR_MAP["WT"], edgecolor='black', label='WT', alpha=0.7),
]
ax.legend(handles=legend_elements, bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig("s_rmsf_bar_vertical_colored_alpha.png", dpi=600, bbox_inches='tight')
plt.show()
