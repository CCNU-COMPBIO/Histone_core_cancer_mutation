"""
3D scatter plot of key residues categorized by folding, PPI, and PDI DDG values.
Input: Excel file with four histone sheets (H2A, H2B, H3, H4).
Output: PNG 3D plot with color-coded categories.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

file_path = 'AlanineScan_DDG_Nucleosome_Averages.xlsx'
sheet_names = ['H2A', 'H2B', 'H3', 'H4']

# Read and combine sheets
df_list = []
for sheet in sheet_names:
    df = pd.read_excel(file_path, sheet_name=sheet)
    df['histone_type'] = sheet
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True)

numeric_cols = ['DDG_GeoStab', 'DDG_MutPNI', 'percentage3', 'DDG_PPI_avg']
for col in numeric_cols:
    combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')

combined_df['GeoStab_abs'] = combined_df['DDG_GeoStab'].abs()
combined_df['MutPNI_abs'] = combined_df['DDG_MutPNI'].abs()
combined_df['PPI_abs'] = combined_df['DDG_PPI_avg'].abs()

# Classify residues based on absolute thresholds
condition_geo = combined_df['GeoStab_abs'] >= 2.0
condition_pdi = (combined_df['MutPNI_abs'] >= 1.0) & (combined_df['percentage3'] > 0)
condition_ppi = combined_df['PPI_abs'] >= 1.5

combined_df['category'] = 'None'
combined_df.loc[condition_geo, 'category'] = 'Geo'
combined_df.loc[condition_pdi, 'category'] = 'PDI'
combined_df.loc[condition_ppi, 'category'] = 'PPI'
combined_df.loc[condition_geo & condition_ppi, 'category'] = 'Geo+PPI'

selected_df = combined_df[combined_df['category'] != 'None'].copy()

# Calculate total absolute DDG for ranking
def pdi_adjust_abs(row):
    if (row['MutPNI_abs'] >= 1.0) and (row['percentage3'] > 0):
        return row['MutPNI_abs']
    return 0.0

selected_df['PDI_adjusted_abs'] = selected_df.apply(pdi_adjust_abs, axis=1)
selected_df['total_deltaG'] = (
    selected_df['GeoStab_abs'].fillna(0) +
    selected_df['PPI_abs'].fillna(0) +
    selected_df['PDI_adjusted_abs']
)

selected_df['residue_id'] = selected_df.apply(
    lambda row: f"{row['histone_type']} {row['wild']}{row['site']}", axis=1
)

category_colors = {
    'Geo': '#ff6f91',
    'PDI': '#0081cf',
    'PPI': '#ffc75f',
    'Geo+PPI': '#845ec2',
}

category_labels = {
    'Geo': r'$\Delta\Delta G_{folding} \geq 2$',
    'PDI': r'$\Delta\Delta G_{PDI} \geq 1$',
    'PPI': r'$\Delta\Delta G_{PPI} \geq 1.5$',
    'Geo+PPI': r'$\Delta\Delta G_{folding} \geq 2$ & $\Delta\Delta G_{PPI} \geq 1.5$',
}

# 3D Plot
plt.rcParams.update({
    'font.size': 16,
    'font.family': 'Arial',
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial:bold',
    'mathtext.it': 'Arial:bold:italic',
})

fig = plt.figure(figsize=(10, 7), constrained_layout=True)
ax = fig.add_subplot(111, projection='3d')

existing_categories = selected_df['category'].unique()

for category in existing_categories:
    if category == 'None':
        continue
    category_data = selected_df[selected_df['category'] == category]
    if len(category_data) == 0:
        continue

    x = category_data['DDG_GeoStab'].fillna(0)
    y = category_data['DDG_PPI_avg'].fillna(0)
    z = category_data['DDG_MutPNI'].fillna(0)
    color = category_colors.get(category, 'gray')
    label = category_labels.get(category, category)

    ax.scatter(x, y, z,
               c=color, marker='o', s=100, alpha=0.8,
               edgecolors='black', linewidth=1.0, label=label)

ax.set_xlabel(r'$\Delta\Delta G_{folding}$ (kcal/mol)', fontsize=20, labelpad=10, fontweight='bold')
ax.set_ylabel(r'$\Delta\Delta G_{PPI}$ (kcal/mol)', fontsize=20, labelpad=10, fontweight='bold')
ax.set_zlabel(r'$\Delta\Delta G_{PDI}$ (kcal/mol)', fontsize=20, labelpad=10, fontweight='bold')

LINEWIDTH = 1.8
for spine in ax.spines.values():
    spine.set_linewidth(LINEWIDTH)
    spine.set_edgecolor('black')

ax.xaxis._axinfo['grid']['linewidth'] = LINEWIDTH
ax.yaxis._axinfo['grid']['linewidth'] = LINEWIDTH
ax.zaxis._axinfo['grid']['linewidth'] = LINEWIDTH
ax.tick_params(axis='x', width=LINEWIDTH)
ax.tick_params(axis='y', width=LINEWIDTH)
ax.tick_params(axis='z', width=LINEWIDTH)
ax.tick_params(axis='both', which='major', labelsize=18)

ax.legend(fontsize=18, loc='upper right', frameon=True, edgecolor='none', bbox_to_anchor=(1.1, 0.90))
ax.view_init(elev=20, azim=30)
ax.grid(True, alpha=0.3)

x_min, x_max = selected_df['DDG_GeoStab'].min(), selected_df['DDG_GeoStab'].max()
y_min, y_max = selected_df['DDG_PPI_avg'].min(), selected_df['DDG_PPI_avg'].max()
z_min, z_max = selected_df['DDG_MutPNI'].min(), selected_df['DDG_MutPNI'].max()

x_pad = (x_max - x_min) * 0.1
y_pad = (y_max - y_min) * 0.1
z_pad = (z_max - z_min) * 0.1

ax.set_xlim(x_min - x_pad, x_max + x_pad)
ax.set_ylim(y_min - y_pad, y_max + y_pad)
ax.set_zlim(z_min - z_pad, z_max + z_pad)

output_filename = 'histone_3d_analysis.png'
plt.savefig(output_filename, dpi=600, bbox_inches='tight', pad_inches=0.5)
print(f"\nFigure saved as '{output_filename}'")

print("\n--- Residue Selection Summary ---")
print(f"Total residues in dataset: {len(combined_df)}")
print(f"Residues passing at least one filter: {len(selected_df)}")

print("\n--- Category Distribution ---")
for cat, count in selected_df['category'].value_counts().items():
    print(f"{cat}: {count} residues")

print("\n--- Top 20 Residues (Highest |ΔΔG| Sum) ---")
top_20_df = selected_df.nlargest(20, 'total_deltaG')[
    ['residue_id', 'histone_type', 'wild', 'site',
     'DDG_GeoStab', 'DDG_PPI_avg', 'DDG_MutPNI', 'percentage3', 'category', 'total_deltaG']
].copy()
top_20_df.columns = ['Residue ID', 'Histone', 'Wild', 'Site',
                     'ΔΔG_folding', 'ΔΔG_PPI', 'ΔΔG_PDI', '%3',
                     'Category', 'Total |ΔΔG|']
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
print(top_20_df.to_string(index=False))

plt.show()