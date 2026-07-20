"""
Plot 2D distribution of key residues (folding, PPI, PDI) based on DDG values.
Input: Excel file with four histone sheets (H2A, H2B, H3, H4).
Output: PNG scatter plot with colorbar for PDI.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from adjustText import adjust_text

file_path = 'AlanineScan_DDG_Nucleosome_Averages.xlsx'
sheet_names = ['H2A', 'H2B', 'H3', 'H4']

# Read and combine sheets
df_list = []
for sheet in sheet_names:
    df = pd.read_excel(file_path, sheet_name=sheet)
    df['histone_type'] = sheet
    df_list.append(df)

combined_df = pd.concat(df_list, ignore_index=True)

# Convert numeric columns
numeric_cols = ['DDG_GeoStab', 'DDG_MutPNI', 'percentage3', 'DDG_PPI_avg']
for col in numeric_cols:
    combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')

# Absolute values for filtering
combined_df['GeoStab_abs'] = combined_df['DDG_GeoStab'].abs()
combined_df['MutPNI_abs'] = combined_df['DDG_MutPNI'].abs()
combined_df['PPI_abs'] = combined_df['DDG_PPI_avg'].abs()

# Classify residues based on absolute thresholds
condition_geo = combined_df['GeoStab_abs'] >= 2.0
condition_pdi = (combined_df['MutPNI_abs'] >= 1.0) & (combined_df['percentage3'] > 0)
condition_ppi = combined_df['PPI_abs'] >= 1.5

combined_df['category'] = ''
combined_df.loc[condition_geo, 'category'] += 'Geo;'
combined_df.loc[condition_pdi, 'category'] += 'PDI;'
combined_df.loc[condition_ppi, 'category'] += 'PPI;'

selected_df = combined_df[combined_df['category'] != ''].copy()

# Calculate total absolute DDG for ranking
def pdi_adjust_abs(row):
    if (row['MutPNI_abs'] >= 1.0) and (row['percentage3'] > 0):
        return row['MutPNI_abs']
    else:
        return 0.0

selected_df['PDI_adjusted_abs'] = selected_df.apply(pdi_adjust_abs, axis=1)

selected_df['total_deltaG'] = (
    selected_df['GeoStab_abs'].fillna(0) +
    selected_df['PPI_abs'].fillna(0) +
    selected_df['PDI_adjusted_abs']
)

selected_df['residue_id'] = selected_df.apply(
    lambda row: f"{row['histone_type']} {row['wild']}{row['site']}",
    axis=1
)

# Top 10 based on total absolute DDG
top_10_indices = selected_df.nlargest(10, 'total_deltaG').index
selected_df['is_top_10'] = selected_df.index.isin(top_10_indices)
top_10_df = selected_df.loc[top_10_indices].copy()

# Plotting
plt.rcParams['font.size'] = 18
fig, ax = plt.subplots(figsize=(10, 7))

# Color mapping for PDI values (using original DDG_MutPNI)
cmap = cm.viridis_r
pdi_for_color = selected_df.loc[selected_df['percentage3'] > 0, 'DDG_MutPNI']
if not pdi_for_color.empty:
    norm = plt.Normalize(vmin=pdi_for_color.min(), vmax=pdi_for_color.max())
else:
    norm = plt.Normalize(vmin=0, vmax=1)

plot_config = {
    'Geo': {'label': r'$\Delta\Delta G_{folding} \geq 2$',  'marker': 'o', 's': 150},
    'PDI': {'label': r'$\Delta\Delta G_{PDI} \geq 1$',      'marker': 'o', 's': 150},
    'PPI': {'label': r'$\Delta\Delta G_{PPI} \geq 1.5$',    'marker': 'o', 's': 150},
}

# Scatter points using original DDG_GeoStab (x) and DDG_PPI_avg (y)
for _, row in selected_df.iterrows():
    x = row['DDG_GeoStab']
    y = row['DDG_PPI_avg']
    cats = row['category']

    if pd.isna(x) or pd.isna(y):
        continue

    if 'Geo' in cats:
        cfg = plot_config['Geo']
    elif 'PDI' in cats:
        cfg = plot_config['PDI']
    elif 'PPI' in cats:
        cfg = plot_config['PPI']
    else:
        continue

    if row['percentage3'] > 0 and not pd.isna(row['DDG_MutPNI']):
        face_color = cmap(norm(row['DDG_MutPNI']))
        alpha = 0.9
        size = cfg['s']
    else:
        face_color = 'lightgray'
        alpha = 0.7
        size = cfg['s'] * 0.8

    edge_color = 'black'
    linewidth = 1.0

    ax.scatter(
        x, y,
        c=[face_color],
        alpha=alpha,
        s=size,
        marker='o',
        edgecolors=edge_color,
        linewidths=linewidth,
        zorder=2
    )

# Axis styling
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')

ax.set_xlabel(r'$\mathbf{\Delta\Delta G_{folding}}$ (kcal/mol)', fontsize=22, fontweight='bold')
ax.set_ylabel(r'$\mathbf{\Delta\Delta G_{PPI}}$ (kcal/mol)', fontsize=22, fontweight='bold')
ax.tick_params(axis='both', which='major', labelsize=20)

# Set limits with padding
x_min = selected_df['DDG_GeoStab'].min()
x_max = selected_df['DDG_GeoStab'].max()
y_min = selected_df['DDG_PPI_avg'].min()
y_max = selected_df['DDG_PPI_avg'].max()
pad_x = (x_max - x_min) * 0.08
pad_y = (y_max - y_min) * 0.08
ax.set_xlim(x_min - pad_x, x_max + pad_x)
ax.set_ylim(y_min - pad_y, y_max + pad_y)

ax.grid(True, linestyle='--', alpha=0.35)
ax.axvline(x=2, color='black', linestyle='--', linewidth=1.8, alpha=0.7)
ax.axhline(y=1.5, color='black', linestyle='--', linewidth=1.8, alpha=0.7)

# Colorbar for PDI
if not pdi_for_color.empty:
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.04)
    cbar.set_label(r'$\mathbf{\Delta\Delta G_{PDI}}$ (kcal/mol)', fontsize=22, fontweight='bold')
    cbar.ax.tick_params(labelsize=18)
else:
    ax.text(0.98, 0.02, 'No residues with %3 > 0 to color by PDI',
            transform=ax.transAxes, fontsize=14, style='italic',
            horizontalalignment='right', verticalalignment='bottom')

plt.tight_layout()
output_filename = 'histone_2d_key_residues.png'
plt.savefig(output_filename, dpi=600, bbox_inches='tight')
print(f"\nFigure saved as '{output_filename}'")

print("\n--- Residue Selection Summary ---")
print(f"Total residues in dataset: {len(combined_df)}")
print(f"Residues passing at least one filter: {len(selected_df)}")

print("\n--- Top 10 Residues (Highest |ΔΔG| Sum) ---")
top_10_display = top_10_df[
    ['residue_id', 'histone_type', 'wild', 'site',
     'DDG_GeoStab', 'DDG_PPI_avg', 'DDG_MutPNI', 'percentage3',
     'PDI_adjusted_abs', 'total_deltaG', 'category']
].copy()
top_10_display.columns = [
    'Residue ID', 'Histone', 'Wild', 'Site',
    'ΔΔG_folding(raw)', 'ΔΔG_PPI(raw)', 'ΔΔG_PDI(raw)', '%3',
    '|ΔΔG_PDI|(used)', 'Total |ΔΔG|', 'Category'
]
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
print(top_10_display.to_string(index=False))

plt.show()