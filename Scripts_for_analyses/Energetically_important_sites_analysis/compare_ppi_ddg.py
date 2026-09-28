# Script name: ppi_ddg_weighted_pcc_scatter.py
# Description: Weighted-average (by mutation frequency) ΔΔG_PPI of cancer mutations
#              vs ΔΔG_PPI of alanine mutations at the same hotspot residues.
#              Computes Pearson correlation coefficient (PCC) and p-value,
#              then draws a vivid scatter plot with |ΔΔG_PPI| = 1.5 threshold
#              lines on both axes.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

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


# ----------------------------
# 1. Read AlanineScan_DDG_Nucleosome_Averages.xlsx
#    使用 DDG_PPI_avg 作为丙氨酸突变的 ΔΔG_PPI
# ----------------------------
file_ala = 'AlanineScan_DDG_Nucleosome_Averages.xlsx'
histone_sheets = ['H2A', 'H2B', 'H3', 'H4']

df_list = []
for sheet in histone_sheets:
    df = pd.read_excel(file_ala, sheet_name=sheet)
    df['histone_type'] = sheet
    df_list.append(df)

df_ala = pd.concat(df_list, ignore_index=True)
df_ala.columns = df_ala.columns.str.strip()

df_ala['frequency'] = pd.to_numeric(df_ala['frequency'], errors='coerce')
df_ala = df_ala[df_ala['frequency'] >= 8].copy()

# 直接用已计算好的平均 PPI ΔΔG
df_ala['PPI_ala'] = pd.to_numeric(df_ala['DDG_PPI_avg'], errors='coerce')

df_ala['site_str'] = df_ala['site'].apply(
    lambda x: str(int(x)) if pd.notnull(x) else ''
)

df_ala['key'] = (
    df_ala['histone_type'] + '_' +
    df_ala['wild'].astype(str).str.strip() +
    df_ala['site_str']
)

# ----------------------------
# 2. Read PPI sheet — cancer mutations with frequency
#    此文件本次未提供，仍保留原文件名
# ----------------------------
df_ppi = pd.read_excel('mutation_hotspot_ddG.xlsx', sheet_name='PPI')
df_ppi.columns = df_ppi.columns.str.strip()

df_ppi['PPI_val'] = pd.to_numeric(df_ppi['PPI'], errors='coerce')
df_ppi['frequency'] = pd.to_numeric(df_ppi['frequency'], errors='coerce')
df_ppi['mutation_prefix'] = df_ppi['mutation'].astype(str).str.strip().str[:-1]
df_ppi['key'] = (
    df_ppi['histone'].astype(str).str.strip() + '_' + df_ppi['mutation_prefix']
)

# ----------------------------
# 3. Collect paired values
# ----------------------------
paired_ala = []
paired_cancer_weighted = []
paired_histone = []
paired_keys = []

for key in df_ala['key'].unique():
    sub = df_ppi[
        (df_ppi['key'] == key) &
        df_ppi['PPI_val'].notna() &
        df_ppi['frequency'].notna()
    ]
    if len(sub) == 0:
        continue

    ala_val = df_ala.loc[df_ala['key'] == key, 'PPI_ala'].values
    if len(ala_val) == 0 or pd.isna(ala_val[0]):
        continue

    w = sub['frequency'].values.astype(float)
    v = sub['PPI_val'].values.astype(float)
    weighted = np.sum(w * v) / np.sum(w)

    paired_ala.append(ala_val[0])
    paired_cancer_weighted.append(weighted)
    paired_histone.append(key.split('_')[0])
    paired_keys.append(key)

paired_ala = np.array(paired_ala)
paired_cancer_weighted = np.array(paired_cancer_weighted)
paired_histone = np.array(paired_histone)
paired_keys = np.array(paired_keys)

print(f"Number of paired hotspot sites: {len(paired_ala)}")
for h in ['H2A', 'H2B', 'H3', 'H4']:
    n = int(np.sum(paired_histone == h))
    if n:
        print(f"  {h}: {n} sites")

# ----------------------------
# 4. Pearson correlation (PCC)
# ----------------------------
r_pcc, p_pcc = pearsonr(paired_ala, paired_cancer_weighted)
print(f"\nPCC (Pearson) = {r_pcc:.3f}, p = {p_pcc:.3e}")

# ----------------------------
# 5. Count |ΔΔG_PPI| >= 1.5 on each axis
# ----------------------------
threshold = 1.5
n_ala_ge = int(np.sum(np.abs(paired_ala) >= threshold))
n_cancer_ge = int(np.sum(np.abs(paired_cancer_weighted) >= threshold))

print(f"\n--- |ΔΔG_PPI| >= {threshold} kcal/mol ---")
print(f"Alanine mutations: {n_ala_ge} / {len(paired_ala)} "
      f"({100 * n_ala_ge / len(paired_ala):.1f}%)")
print(f"Cancer mutations (frequency-weighted): {n_cancer_ge} / {len(paired_cancer_weighted)} "
      f"({100 * n_cancer_ge / len(paired_cancer_weighted):.1f}%)")

print("\nSites with |alanine ΔΔG_PPI| >= 1.5 (sorted by value, descending):")
order = np.argsort(paired_ala)[::-1]
for idx in order:
    if abs(paired_ala[idx]) >= threshold:
        print(f"  {paired_keys[idx]:<10s}  "
              f"alanine ΔΔG_PPI = {paired_ala[idx]:+.3f}  |  "
              f"weighted cancer ΔΔG_PPI = {paired_cancer_weighted[idx]:+.3f}")

# ----------------------------
# 6. Scatter plot with |ΔΔG_PPI| = 1.5 threshold lines on both axes
# ----------------------------
fig, ax = plt.subplots(figsize=(8, 7))

POINT_COLOR = '#FF3B1F'
THRESH_COLOR = '#B22222'

ax.scatter(
    paired_ala, paired_cancer_weighted,
    s=120, color=POINT_COLOR, edgecolor='black', linewidth=1.8,
    alpha=0.95, zorder=3
)

# Linear regression line
if len(paired_ala) > 2:
    z = np.polyfit(paired_ala, paired_cancer_weighted, 1)
    xline = np.linspace(paired_ala.min(), paired_ala.max(), 100)
    ax.plot(xline, z[0] * xline + z[1], 'k--', linewidth=2.2, zorder=2)

# y = x reference line
all_vals = np.concatenate([paired_ala, paired_cancer_weighted])
pad = 0.05 * (all_vals.max() - all_vals.min())
lim_lo = all_vals.min() - pad
lim_hi = all_vals.max() + pad
ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi],
        color='gray', linestyle=':', linewidth=1.8, zorder=1)

# Threshold lines at ±1.5 on both axes
for x_thr in [threshold, -threshold]:
    ax.axvline(x_thr, color=THRESH_COLOR, linestyle='--',
               linewidth=1.8, alpha=0.7, zorder=2)

for y_thr in [threshold, -threshold]:
    ax.axhline(y_thr, color=THRESH_COLOR, linestyle='--',
               linewidth=1.8, alpha=0.7, zorder=2)

ax.set_xlim(lim_lo, lim_hi)
ax.set_ylim(lim_lo, lim_hi)

# Labels / title
ax.set_xlabel(
    r'$\mathbf{\Delta\Delta G_{PPI}}$ of alanine mutation (kcal/mol)',
    fontsize=16, fontweight='bold', fontname='Arial', labelpad=10
)
ax.set_ylabel(
    r'$\mathbf{\Delta\Delta G_{PPI}}$ of cancer mutations (frequency-weighted, kcal/mol)',
    fontsize=16, fontweight='bold', fontname='Arial', labelpad=10
)
ax.set_title(
    r'$\mathbf{\Delta\Delta G_{PPI}}$ correlation: alanine vs cancer mutations',
    fontsize=16, fontweight='bold', fontname='Arial', pad=14
)

ax.tick_params(axis='both', which='major', labelsize=14,
               width=2.5, length=8, colors='black')
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontname('Arial')
    label.set_fontweight('normal')

# PCC + p-value annotation
if p_pcc < 1e-4:
    p_str = f"$p$ = {p_pcc:.2e}"
else:
    p_str = f"$p$ = {p_pcc:.4f}"

corr_text = f"PCC = {r_pcc:.2f}\n{p_str}\n$n$ = {len(paired_ala)}"
ax.text(
    0.03, 0.97, corr_text,
    transform=ax.transAxes, ha='left', va='top',
    fontsize=14, fontname='Arial',
    color='black',
    bbox=dict(boxstyle='round,pad=0.4',
              facecolor='white', edgecolor='none', alpha=0.85)
)

for spine in ax.spines.values():
    spine.set_linewidth(2.5)
    spine.set_color('black')

plt.tight_layout()
plt.savefig('ppi_ddg_weighted_pcc_scatter.png', dpi=600, bbox_inches='tight')
plt.show()