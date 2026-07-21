import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the input CSV with standardized column names
csv_path = "SRMSF_PseudoSymmetry_Deviation_Trim5_PerRunDetail.csv"
df_raw = pd.read_csv(csv_path)

# Build y-axis labels: WT + (Histone_Subunit + Variant)
labels = ["WT"] + [f"{row['Histone_Subunit']}{row['Variant']}" for _, row in df_raw.iterrows()]

# Extract mean and SEM values using new column names
means = [df_raw["SRMSF_WT_Mean_A"].iloc[0]] + df_raw["SRMSF_Variant_Mean_A"].tolist()
sems  = [df_raw["SRMSF_WT_SEM_A"].iloc[0]]  + df_raw["SRMSF_Variant_SEM_A"].tolist()

S_WT = means[0]
y_pos = np.arange(len(labels))

# Color scheme: light gray for WT, soft blue for mutants
WT_COLOR = "#E8E8E8"
MUT_COLOR = "#7EA6D9"
bar_colors = [WT_COLOR] + [MUT_COLOR] * len(df_raw)

# Global font settings (Arial)
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.unicode_minus': False,
    'font.size': 20,
    'axes.labelsize': 20,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
})

# Create horizontal bar plot with error bars (SEM)
fig, ax = plt.subplots(figsize=(7, 6))
bars = ax.barh(
    y_pos, means, height=0.55, xerr=sems,
    color=bar_colors,
    edgecolor="black", linewidth=1.5,
    error_kw=dict(ecolor='#333333', elinewidth=1.5, capsize=4, capthick=1.5),
    zorder=2
)

# Dashed vertical line for WT reference
ax.axvline(x=S_WT, color='#666666', linewidth=1.0, linestyle='--', zorder=1)

# Axis labels and limits
ax.set_yticks(y_pos)
ax.set_yticklabels(labels, va='center')
ax.set_xlabel(r"$\mathbf{S_{RMSF}}$ (Å)", fontweight='bold', fontsize=20)

x_max = max(np.nanmax(means) + max(sems), S_WT) * 1.15
ax.set_xlim(0, x_max)

# Thick black spines
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
    spine.set_edgecolor('black')

plt.tight_layout()
plt.savefig("s_rmsf_bar_sem_only.png", dpi=300, bbox_inches='tight')
plt.show()