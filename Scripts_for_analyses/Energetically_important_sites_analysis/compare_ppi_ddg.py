# Script name: ppi_violin_compare.py
# Compare two |ΔΔG_PPI| distributions with a single violin plot.
# Group 1 — all cancer mutations from mutation_hotspot_ddG.xlsx (PPI sheet, no weighting)
# Group 2 — 3afa hotspot sites with |ΔΔG_PPI| >= 1.5
# Both groups are plotted as absolute values.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Global style
# ----------------------------
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Arial'
plt.rcParams['mathtext.it'] = 'Arial:italic'
plt.rcParams['mathtext.bf'] = 'Arial:bold'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['xtick.major.width'] = 2.5
plt.rcParams['ytick.major.width'] = 2.5
plt.rcParams['xtick.major.size'] = 8
plt.rcParams['ytick.major.size'] = 8
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['text.color'] = 'black'
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'

THRESHOLD = 1.5

# ----------------------------
# 1. Group 1 — all cancer mutations, absolute value
# ----------------------------
df_ppi = pd.read_excel('mutation_hotspot_ddG.xlsx', sheet_name='PPI')
df_ppi.columns = df_ppi.columns.str.strip()
df_ppi['PPI_val'] = pd.to_numeric(df_ppi['PPI'], errors='coerce')

signed1 = df_ppi['PPI_val'].dropna()
group1 = signed1.abs().to_numpy()

# ----------------------------
# 2. Group 2 — 3afa |ΔΔG_PPI| >= 1.5, absolute value
# ----------------------------
histone_sheets = ['H2A', 'H2B', 'H3', 'H4']
df_list = []
for sheet in histone_sheets:
    d = pd.read_excel('3afa.xlsx', sheet_name=sheet)
    d['histone_type'] = sheet
    df_list.append(d)

df_3afa = pd.concat(df_list, ignore_index=True)
df_3afa.columns = df_3afa.columns.str.strip()
df_3afa['PPI_val'] = pd.to_numeric(df_3afa['PPI'], errors='coerce')

signed2 = df_3afa['PPI_val'].dropna()
group2 = signed2[signed2.abs() >= THRESHOLD].abs().to_numpy()

# Check for empty data
if len(group1) == 0 or len(group2) == 0:
    raise ValueError("At least one group has no usable data; check the Excel files and threshold.")

# ----------------------------
# 3. Print counts and stats
# ----------------------------
print(f"Group 1 — all cancer mutations |ΔΔG_PPI| : n = {len(group1)}")
print(f"Group 2 — 3afa |ΔΔG_PPI| >= {THRESHOLD} : n = {len(group2)}")
print(f"  group 1: mean = {group1.mean():.3f}, median = {np.median(group1):.3f}")
print(f"  group 2: mean = {group2.mean():.3f}, median = {np.median(group2):.3f}")

# Count positive/negative values in the original signed PPI of cancer mutations
n_ge = int(np.sum(signed1 >= THRESHOLD))
n_le = int(np.sum(signed1 <= -THRESHOLD))
print(f"\nCancer mutations with ΔΔG_PPI >= +{THRESHOLD}: {n_ge} / {len(signed1)} "
      f"({100.0 * n_ge / len(signed1):.1f}%)")
print(f"Cancer mutations with ΔΔG_PPI <= -{THRESHOLD}: {n_le} / {len(signed1)} "
      f"({100.0 * n_le / len(signed1):.1f}%)")

# Count positive/negative values in the original signed PPI of 3afa
n2_ge = int(np.sum(signed2 >= THRESHOLD))
n2_le = int(np.sum(signed2 <= -THRESHOLD))
print(f"3afa mutations with ΔΔG_PPI >= +{THRESHOLD}: {n2_ge} / {len(signed2)} "
      f"({100.0 * n2_ge / len(signed2):.1f}%)")
print(f"3afa mutations with ΔΔG_PPI <= -{THRESHOLD}: {n2_le} / {len(signed2)} "
      f"({100.0 * n2_le / len(signed2):.1f}%)")

# ----------------------------
# 4. Violin plot
# ----------------------------
fig, ax = plt.subplots(figsize=(10, 7))

data = [group1, group2]
positions = [1, 2]

# NPG palette
colors = ['#4DBBD5', '#E64B35']

vp = ax.violinplot(
    data, positions=positions, widths=0.75,
    showmeans=False, showmedians=True, showextrema=False
)

for body, c in zip(vp['bodies'], colors):
    body.set_facecolor(c)
    body.set_edgecolor('black')
    body.set_linewidth(1.8)
    body.set_alpha(0.80)

vp['cmedians'].set_color('black')
vp['cmedians'].set_linewidth(2.2)

# Reference lines: 0 (gray dotted) and +1.5 (red dashed)
ax.axhline(0, color='gray', linestyle=':', linewidth=1.4, alpha=0.8, zorder=1)
ax.axhline(THRESHOLD, color='#B22222', linestyle='--',
           linewidth=1.6, alpha=0.75, zorder=1)

# x-axis labels
ax.set_xticks(positions)
ax.set_xticklabels([
    'Cancer mutations',
    r'$|\Delta\Delta G_{PPI}| \geq 1.5$'
], fontsize=15)

# y-axis label
ax.set_ylabel(
    r'$\mathbf{|\Delta\Delta G_{PPI}|}$' + ' ' + r'$\mathbf{(kcal/mol)}$',
    fontsize=16, labelpad=10
)

# Title
ax.set_title(
    r'$\mathbf{|\Delta\Delta G_{PPI}|}$ distribution comparison',
    fontsize=16, fontweight='bold', pad=14
)

ax.tick_params(axis='y', labelsize=14, width=2.5, length=8)
ax.set_xlim(0.4, 2.6)

for label in ax.get_yticklabels():
    label.set_fontname('Arial')

for spine in ax.spines.values():
    spine.set_linewidth(2.5)
    spine.set_color('black')

plt.tight_layout()
plt.savefig('ppi_violin_compare.png', dpi=600, bbox_inches='tight')
plt.show()
