"""
Plot GeoStab (DDG_GeoStab) values for buried/internal residues in each histone.
Input: Excel file with four histone sheets (H2A, H2B, H3, H4).
Output: TIFF plots for internal residues with GeoStab labels.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
from adjustText import adjust_text

plt.rcParams['font.family'] = 'Arial'

sns.set(style="whitegrid", palette="pastel", context="paper", font_scale=1.1)

os.makedirs('histone_plots_final', exist_ok=True)

file_path = 'AlanineScan_DDG_Nucleosome_Averages.xlsx'
sheets = {
    'H2A': pd.read_excel(file_path, sheet_name='H2A'),
    'H2B': pd.read_excel(file_path, sheet_name='H2B'),
    'H3': pd.read_excel(file_path, sheet_name='H3'),
    'H4': pd.read_excel(file_path, sheet_name='H4')
}

conditions = {
    'internal': ('percentage1 < 0.2', 'Buried Residues', '#FF0000')
}

for hist_name, df in sheets.items():
    for col in ['percentage1', 'percentage2', 'percentage3', 'DDG_GeoStab']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['site'] = pd.to_numeric(df['site'], errors='coerce')
    df['site_int'] = df['site'].fillna(-1).astype(int)
    df['label'] = df['wild'] + df['site_int'].astype(str)
    df.loc[df['site_int'] == -1, 'site_int'] = np.nan

for hist_name, df in sheets.items():
    df = df.sort_values('site').dropna(subset=['site', 'DDG_GeoStab'])

    for cond_key, (cond_query, cond_title, color) in conditions.items():
        plt.figure(figsize=(18, 5))

        plt.plot(df['site'], df['DDG_GeoStab'],
                 'o-', color='#888888',
                 linewidth=1.5, markersize=8,
                 markerfacecolor='white', markeredgewidth=1.5,
                 alpha=0.7, label='All Residues')

        cond_df = df.query(cond_query).copy()
        if not cond_df.empty:
            plt.scatter(cond_df['site'], cond_df['DDG_GeoStab'],
                        color=color, s=100, zorder=5,
                        edgecolor='black', linewidth=1.2,
                        label=cond_title)

        texts = []
        for _, row in df.iterrows():
            if row['DDG_GeoStab'] >= 2:
                texts.append(plt.text(
                    row['site'], row['DDG_GeoStab'] + 0.05,
                    row['label'],
                    fontsize=22,
                    color='black',
                    fontweight='bold',
                    ha='center',
                    va='bottom'
                ))
        for _, row in df.iterrows():
            if row['DDG_GeoStab'] <= -2:
                texts.append(plt.text(
                    row['site'], row['DDG_GeoStab'] - 0.05,
                    row['label'],
                    fontsize=22,
                    color='black',
                    fontweight='bold',
                    ha='right',
                    va='top'
                ))

        adjust_text(texts)

        ax = plt.gca()
        for spine in ax.spines.values():
            spine.set_linewidth(2)
            spine.set_color('black')

        plt.title(f'{hist_name} {cond_title} - GeoStab', fontsize=24, pad=10, weight='bold')
        plt.xlabel('Residue Position', fontsize=22, weight='bold')
        plt.ylabel('ΔΔG (kcal/mol)', fontsize=22, weight='bold')

        y_min = df['DDG_GeoStab'].min()
        y_max = df['DDG_GeoStab'].max()
        y_pad = (y_max - y_min) * 0.1
        plt.ylim(y_min - y_pad, y_max + y_pad)

        plt.grid(True, linestyle='--', alpha=0.3)
        plt.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1)
        plt.axhline(y=2, color='black', linestyle='-', alpha=0.3, linewidth=1)
        plt.axhline(y=-2, color='black', linestyle='-', alpha=0.3, linewidth=1)

        plt.legend(loc='lower right', frameon=True, framealpha=0.9,
                   prop={'size': 20})
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)

        plt.tight_layout()

        plt.savefig(f'histone_plots_final/{hist_name}_{cond_key}-GeoStab.tif',
                    dpi=600, bbox_inches='tight', format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
        plt.close()

print("All GeoStab plots for buried residues saved to 'histone_plots_final' folder!")