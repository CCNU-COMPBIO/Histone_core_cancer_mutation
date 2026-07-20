import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ========== 1. Data loading and preparation ==========
csv_path = "rmsf_symmetry_deviation_summary_trim5.csv"
df_raw = pd.read_csv(csv_path)

# Build labels: WT + (Histone + Mutant)
labels = ["WT"] + [f"{row['Histone']}{row['Mutant']}" for _, row in df_raw.iterrows()]

# Extract means and SEMs
means = [df_raw["S_WT"].iloc[0]] + df_raw["S_Mutant"].tolist()
sems = [df_raw["SEM_WT"].iloc[0]] + df_raw["SEM_Mutant"].tolist()

S_WT = means[0]
y_pos = np.arange(len(labels))

# ========== 2. Unified colors (no gradient) ==========
WT_COLOR = "#E8E8E8"       # Light gray for WT
MUT_COLOR = "#7EA6D9"      # Uniform soft blue for mutants
bar_colors = [WT_COLOR] + [MUT_COLOR] * len(df_raw)

# ========== 3. Global font settings ==========
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.unicode_minus': False,
    'font.size': 20,
    'axes.labelsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
})

# ========== 4. Plotting ==========
fig, ax = plt.subplots(figsize=(7, 6))

# Bar plot with SEM error bars
bars = ax.barh(
    y_pos, means, height=0.55, xerr=sems,
    color=bar_colors,
    edgecolor="black", linewidth=1.5,
    error_kw=dict(ecolor='#333333', elinewidth=1.5, capsize=4, capthick=1.5),
    zorder=2
)

# WT reference line
ax.axvline(x=S_WT, color='#666666', linewidth=1.0, linestyle='--', zorder=1)

# ========== 5. Axes and borders ==========
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, va='center')
ax.set_xlabel(r"$\mathbf{S_{RMSF}}$ (Å)", fontweight='bold', fontsize=20)

# X-axis limit: accommodate mean + SEM and WT
x_max = max(np.nanmax(means) + max(sems), S_WT) * 1.15
ax.set_xlim(0, x_max)

for spine in ax.spines.values():
    spine.set_linewidth(2.0)
    spine.set_edgecolor('black')

plt.tight_layout()
plt.savefig("srmsf_barplot.png", dpi=300, bbox_inches='tight')
plt.show()