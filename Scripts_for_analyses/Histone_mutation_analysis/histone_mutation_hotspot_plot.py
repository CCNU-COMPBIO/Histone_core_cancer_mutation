import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

# ==================== Core region definitions ====================
core_regions = {
    "H2A": (14, 118),
    "H2B": (24, 125),
    "H3": (37, 135),
    "H4": (21, 102)
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
    """Extract the numeric position from residue_real_site."""
    if isinstance(residue_str, str) and residue_str.strip():
        match = re.search(r'\d+', residue_str)
        if match:
            return int(match.group())
    return None

def extract_residue(residue_str):
    """Extract the one-letter residue code from residue_real_site."""
    if isinstance(residue_str, str) and residue_str.strip():
        for char in residue_str:
            if char.isalpha():
                return char
    return None

def extract_frequency(freq_str, residue_char):
    """Extract the frequency for a specific residue from residue_frequency."""
    if not isinstance(freq_str, str) or not residue_char:
        return 0

    parts = freq_str.split(',')
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
    """Process a single sheet: filter by conservation, extract positions and frequencies."""
    # Filter by conservation >= 0.4
    df = df[df['conservation'] >= 0.4].copy()

    # Extract position from residue_real_site
    df['position'] = df['residue_real_site'].apply(extract_position)
    df = df.dropna(subset=['position'])

    # Filter to core region
    if sheet_name in core_regions:
        start, end = core_regions[sheet_name]
        df = df[(df['position'] >= start) & (df['position'] <= end)].copy()

    # Extract residue type and specific frequency
    df['residue_type'] = df['residue_real_site'].apply(extract_residue)
    df['specific_freq'] = df.apply(
        lambda row: extract_frequency(row['residue_frequency'], row['residue_type']),
        axis=1
    )

    # Keep only rows with frequency > 0
    df = df[df['specific_freq'] > 0].copy()
    df['position'] = df['position'].astype(int)

    return df.sort_values('position').reset_index(drop=True)

# ==================== Main program ====================
file_path = "Histone_mutation.xlsx"
all_sheets = pd.read_excel(file_path, sheet_name=None)

for sheet_name, df in all_sheets.items():
    if sheet_name not in core_regions:
        continue

    processed_df = process_sheet(df, sheet_name)

    # Define high-frequency mutations as those with frequency >= 8
    threshold = 8
    top_points = processed_df[processed_df['specific_freq'] >= threshold]
    non_top_points = processed_df[processed_df['specific_freq'] < threshold]

    plt.figure(figsize=(10, 6))

    # Plot vertical lines for non-top sites (blue)
    for _, row in non_top_points.iterrows():
        plt.vlines(
            x=row['position'],
            ymin=0,
            ymax=row['specific_freq'],
            colors='steelblue',
            linewidth=1.5,
            alpha=0.6
        )

    # Plot vertical lines for top sites (red, highlighted)
    for _, row in top_points.iterrows():
        plt.vlines(
            x=row['position'],
            ymin=0,
            ymax=row['specific_freq'],
            colors='firebrick',
            linewidth=2.0,
            alpha=0.8
        )

    # Scatter for non-top sites (blue)
    if not non_top_points.empty:
        sns.scatterplot(
            data=non_top_points,
            x='position',
            y='specific_freq',
            color='steelblue',
            s=100,
            edgecolor='white',
            linewidth=0.8,
            alpha=0.8,
            legend=False
        )

    # Scatter for top sites (red, highlighted)
    if not top_points.empty:
        sns.scatterplot(
            data=top_points,
            x='position',
            y='specific_freq',
            color='firebrick',
            s=120,
            edgecolor='white',
            linewidth=0.8,
            alpha=1.0,
            legend=False
        )

    # Optional: annotate top sites (currently commented out)
    # texts = []
    # for _, row in top_points.iterrows():
    #     texts.append(plt.text(
    #         x=row['position'],
    #         y=row['specific_freq'] + max_freq * 0.02,
    #         s=row['residue_real_site'],
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

    max_freq = processed_df['specific_freq'].max()
    plt.title(f"{sheet_name} Mutation Frequency (Frequency ≥ {threshold}, {len(top_points)} sites)",
              fontsize=24, pad=20, fontweight='bold')
    plt.xlabel("Mutation Position", fontsize=22, labelpad=15, fontweight='bold')
    plt.ylabel("Mutation Frequency", fontsize=22, labelpad=15, fontweight='bold')

    plt.ylim(0, max_freq * 1.1)
    plt.yticks(range(0, int(max_freq) + 5, 5), fontsize=20)

    # Set x-axis limits to core region
    plt.xlim(core_regions[sheet_name][0] - 1, core_regions[sheet_name][1] + 1)
    plt.xticks(range(core_regions[sheet_name][0], core_regions[sheet_name][1] + 1, 20),
               fontsize=20)

    plt.tight_layout()
    plt.savefig(f"{sheet_name}_mutation_frequency_top10-nature.tif", dpi=600, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

print("All hotspot plots generated successfully!")