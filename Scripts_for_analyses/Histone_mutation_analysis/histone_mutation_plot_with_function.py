import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
import re
import numpy as np

# ==================== Core region definitions ====================
core_regions = {
    "H2A": (14, 118),
    "H2B": (24, 125),
    "H3": (37, 135),
    "H4": (21, 102)
}

# ==================== Literature-reported key sites ====================
literature_positions = {
    "H2A": [4, 29, 56, 61, 74, 75, 81, 88, 90, 92, 121],
    "H2B": [2, 53, 64, 68, 70, 71, 76, 113],
    "H3": [27, 34, 36, 50, 97, 105],
    "H4": [42, 45, 68, 91, 92]
}

# ==================== Visualization settings ====================
sns.set_context("paper", font_scale=1.8, rc={
    "lines.linewidth": 2.5,
    "axes.labelweight": "bold",
    "font.family": "Arial"
})

sns.set_style("whitegrid", {
    "grid.color": ".8",
    "grid.linestyle": "--",
    "axes.edgecolor": "0.2"
})

# ==================== Helper functions ====================
def extract_position(residue_str):
    """Extract the numeric position from residue_real_site string."""
    if pd.notna(residue_str) and residue_str.strip():
        match = re.search(r'\d+', str(residue_str))
        return int(match.group()) if match else None
    return None

def extract_residue(residue_str):
    """Extract the one-letter residue code from residue_real_site."""
    if pd.notna(residue_str) and residue_str.strip():
        for char in residue_str:
            if char.isalpha():
                return char
    return None

def extract_frequency(freq_str, residue_char):
    """Extract the frequency for a specific residue from the residue_frequency string."""
    if not pd.notna(freq_str) or not residue_char:
        return 0

    if ',' in freq_str:
        parts = freq_str.split(',')
    elif ';' in freq_str:
        parts = freq_str.split(';')
    else:
        parts = [freq_str]

    for part in parts:
        res_freq = part.strip().split(':')
        if len(res_freq) == 2 and res_freq[0].strip() == residue_char:
            try:
                return int(res_freq[1].strip())
            except ValueError:
                return 0
    return 0


# ==================== Data processing function ====================
def process_sheet(df, sheet_name):
    # Standardize column names to lowercase
    df.columns = df.columns.str.strip().str.lower()

    # Extract position numbers
    df['pos'] = df['residue_real_site'].apply(extract_position)

    # Get core region boundaries (for literature site marking)
    start, end = core_regions.get(sheet_name, (0, float('inf')))

    # Filter by conservation (if column exists)
    if 'conservation' in df.columns:
        df['conservation'] = pd.to_numeric(df['conservation'], errors='coerce')
        df = df.dropna(subset=['conservation'])
        df = df[df['conservation'] >= 0.4]

    # Extract residue type and specific frequency
    df['residue_type'] = df['residue_real_site'].apply(extract_residue)
    df['specific_freq'] = df.apply(
        lambda row: extract_frequency(row['residue_frequency'], row['residue_type']),
        axis=1
    )

    # Remove rows with zero frequency
    df = df[df['specific_freq'] > 0].copy()

    # Mark literature sites within core regions
    lit_sites = literature_positions.get(sheet_name, [])
    if sheet_name in core_regions:
        start, end = core_regions[sheet_name]
        df['in_core'] = (df['pos'] >= start) & (df['pos'] <= end)
        df['is_key_site'] = df.apply(
            lambda row: 'literature' if (row['in_core'] and row['pos'] in lit_sites) else 'others',
            axis=1
        )
    else:
        df['is_key_site'] = 'others'

    # Sort by frequency and assign rank for x-axis
    df = df.sort_values('specific_freq', ascending=True)
    df['x_rank'] = range(1, len(df) + 1)

    return df


# ==================== Main program ====================
all_sheets = pd.read_excel("Histone_mutation.xlsx", sheet_name=None)

for sheet_name, df in all_sheets.items():
    if sheet_name not in literature_positions:  # Skip sheets without literature annotations
        continue

    processed_df = process_sheet(df, sheet_name)

    # Extract literature sites within core regions
    key_points = processed_df[processed_df['is_key_site'] == 'literature']

    plt.figure(figsize=(10, 6))

    # Plot non-key sites
    non_key = processed_df[processed_df['is_key_site'] != 'literature']
    key_points_plot = processed_df[processed_df['is_key_site'] == 'literature']

    if not non_key.empty:
        sns.scatterplot(
            data=non_key,
            x='x_rank',
            y='specific_freq',
            color='steelblue',
            s=100,
            edgecolor='white',
            linewidth=0.8,
            alpha=0.8
        )

    if not key_points_plot.empty:
        sns.scatterplot(
            data=key_points_plot,
            x='x_rank',
            y='specific_freq',
            color='firebrick',
            s=100,
            edgecolor='white',
            linewidth=0.8,
            alpha=1.0,
            label='Literature reported'
        )

    plt.legend(loc='upper left', fontsize=20, frameon=False)

    # Optional: annotate literature sites (currently commented out)
    # texts = []
    # for _, row in key_points.iterrows():
    #     texts.append(plt.text(
    #         x=row['x_rank'],
    #         y=row['specific_freq'] + 0.2,
    #         s=f"{row['residue_real_site']}",
    #         color='firebrick',
    #         fontsize=22,
    #         ha='center',
    #         weight='bold',
    #         rotation=90
    #     ))
    # if texts:
    #     adjust_text(texts)

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    plt.title(f"{sheet_name} Mutation Frequency (Sorted)",
              fontsize=24, pad=20, fontweight='bold')
    plt.xlabel("Mutation Position (Sorted by Frequency)", fontsize=22, fontweight='bold')
    plt.ylabel("Mutation Frequency", fontsize=22, fontweight='bold')

    max_freq = processed_df['specific_freq'].max()
    plt.ylim(0, max_freq * 1.1)
    plt.yticks(range(0, int(max_freq) + 5, 5), fontsize=20)
    plt.xticks([])  # Hide x-axis ticks

    plt.tight_layout()
    plt.savefig(f"{sheet_name}_mutation_frequency_sorted-nature.tif", dpi=600, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

print("All plots generated successfully!")