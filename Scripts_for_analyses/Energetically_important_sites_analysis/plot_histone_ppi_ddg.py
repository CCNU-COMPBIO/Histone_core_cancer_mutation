"""
Plot mean ΔΔG values from three predictors (RDE-network, SAAMBE-3D, DDMut-PPI)
for each histone, with error bars and highlight histone-histone interface residues.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from adjustText import adjust_text

plt.rcParams['font.family'] = 'Arial'

os.makedirs('histone_plots_final', exist_ok=True)

file_path = 'AlanineScan_DDG_Nucleosome_Averages.xlsx'
sheets = {
    'H2A': pd.read_excel(file_path, sheet_name='H2A'),
    'H2B': pd.read_excel(file_path, sheet_name='H2B'),
    'H3': pd.read_excel(file_path, sheet_name='H3'),
    'H4': pd.read_excel(file_path, sheet_name='H4')
}

for hist_name, df in sheets.items():
    for col in ['percentage1', 'percentage2', 'percentage3', 'DDG_RDE-network', 'DDG_SAAMBE-3D', 'DDG_DDMut-PPI', 'site']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['site_int'] = df['site'].fillna(-1).astype(int)
    df['label'] = df['wild'] + df['site_int'].astype(str)
    df.loc[df['site_int'] == -1, 'site_int'] = np.nan

for hist_name, df in sheets.items():
    df = df.sort_values('site').dropna(subset=['site', 'DDG_RDE-network', 'DDG_SAAMBE-3D', 'DDG_DDMut-PPI'])

    if df.empty:
        print(f"{hist_name} has no valid data")
        continue

    df['mean_ddg'] = df[['DDG_RDE-network', 'DDG_SAAMBE-3D', 'DDG_DDMut-PPI']].mean(axis=1)
    df['std_ddg'] = df[['DDG_RDE-network', 'DDG_SAAMBE-3D', 'DDG_DDMut-PPI']].std(axis=1)
    df['sem_ddg'] = df['std_ddg'] / np.sqrt(3)

    plt.figure(figsize=(18, 5))

    plt.errorbar(df['site'], df['mean_ddg'],
                 yerr=df['sem_ddg'],
                 fmt='o-', color='#888888', linewidth=1.5, markersize=8,
                 markerfacecolor='white', markeredgewidth=1.5,
                 ecolor='#FA7F08', elinewidth=1.5, capsize=5, capthick=1.5,
                 alpha=0.8, label='Mean ΔΔG')

    interface_df = df.query('percentage2 > 0').copy()

    if not interface_df.empty:
        plt.scatter(interface_df['site'], interface_df['mean_ddg'],
                    color='#F24405', s=100, zorder=5,
                    edgecolor='black', linewidth=1.5,
                    label='Histone-Histone Interface')

    texts = []
    for _, row in df.iterrows():
        if abs(row['mean_ddg']) >= 1.5:
            label = f"{row['wild']}{int(row['site'])}"
            texts.append(plt.text(row['site'], row['mean_ddg'] + 0.05,
                         label, ha='center', va='bottom', fontsize=22,
                         color='black', fontweight='bold'))

    adjust_text(texts)

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    plt.title(f'{hist_name} - Mean ΔΔG Values',
              fontsize=24, pad=10, weight='bold')
    plt.xlabel('Residue Position', fontsize=22, weight='bold')
    plt.ylabel('ΔΔG (kcal/mol)', fontsize=22, weight='bold')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    plt.grid(True, linestyle='--', alpha=0.4)

    plt.axhline(y=0, color='black', linestyle='-', alpha=0.5, linewidth=1.5)
    plt.axhline(y=1.5, color='black', linestyle='-', alpha=0.5, linewidth=1.5)

    plt.legend(loc='upper left', frameon=True, framealpha=0.5,
               prop={'size': 20})

    y_min = min(df['mean_ddg'] - df['sem_ddg']) - 0.2
    y_max = max(df['mean_ddg'] + df['sem_ddg']) + 0.5
    plt.ylim(y_min, y_max)

    plt.tight_layout()

    filename = f'histone_plots_final/{hist_name}_RDE_SAAMBE_DDMut_mean.tif'
    plt.savefig(filename, dpi=600, bbox_inches='tight', format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    print(f"Saved plot: {filename}")
    plt.close()

    print(f"{hist_name} - Total residues: {len(df)}")
    if not interface_df.empty:
        print(f"{hist_name} - Histone-histone interface residues: {len(interface_df)}")
    print(f"{hist_name} - Mean ΔΔG: {df['mean_ddg'].mean():.3f} ± {df['mean_ddg'].std():.3f} kcal/mol")
    print("-" * 50)

print("All mean ΔΔG plots saved to 'histone_plots_final' folder!")