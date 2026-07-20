"""
Plot DDG (MutPNI) values for histone-DNA interface residues.
Input: Excel file with four histone sheets (H2A, H2B, H3, H4).
Output: TIFF bar plots for each histone.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

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

for hist_name, df in sheets.items():
    for col in ['percentage1', 'percentage2', 'percentage3', 'DDG_MutPNI']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['site'] = pd.to_numeric(df['site'], errors='coerce')
    df['site_int'] = df['site'].fillna(-1).astype(int)
    df['label'] = df['wild'] + df['site_int'].astype(str)
    df.loc[df['site_int'] == -1, 'site_int'] = np.nan

for hist_name, df in sheets.items():
    df = df.sort_values('site').dropna(subset=['site', 'DDG_MutPNI'])
    dna_interface_df = df.query('percentage3 > 0').copy()

    if dna_interface_df.empty:
        print(f"{hist_name} has no DNA-interface residues with percentage3 > 0")
        continue

    max_bars = max(
        len(df.query('percentage3 > 0'))
        for df in sheets.values()
    )

    FIXED_FIG_WIDTH = 18
    BASE_BAR_WIDTH = 0.8

    n_bars = len(dna_interface_df)
    x_pos = np.arange(n_bars)
    dynamic_width = BASE_BAR_WIDTH * (n_bars / max_bars)

    plt.figure(figsize=(FIXED_FIG_WIDTH, 5))
    bars = plt.bar(x_pos, dna_interface_df['DDG_MutPNI'],
                   color='#96bbfb',
                   width=max(dynamic_width, 0.1),
                   label='ΔΔG (MutPNI)')

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    plt.title(f'{hist_name} Histone-DNA Interface',
              fontsize=24, pad=20, weight='bold')
    plt.ylabel('ΔΔG (kcal/mol)', fontsize=22, weight='bold')

    plt.xticks(x_pos, dna_interface_df['label'], rotation=90, ha='center',
               fontsize=22, weight='bold')
    plt.yticks(fontsize=20)

    plt.grid(True, linestyle='--', alpha=0.4, axis='y')
    plt.grid(axis='x', visible=False)

    plt.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1)
    plt.axhline(y=1, color='black', linestyle='-', alpha=0.3, linewidth=1)

    y_min = dna_interface_df['DDG_MutPNI'].min() - 0.2
    y_max = dna_interface_df['DDG_MutPNI'].max() + 0.5
    plt.ylim(y_min, y_max)

    plt.tight_layout()

    plt.savefig(f'histone_plots_final/{hist_name}_DNA_interface_ddg.tif',
                dpi=600, bbox_inches='tight', format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

print("All histone-DNA interface DDG bar plots saved to 'histone_plots_final' folder!")