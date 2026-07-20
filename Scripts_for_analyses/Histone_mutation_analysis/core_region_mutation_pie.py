import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re

# ===== 1. Load multiple sheets and combine =====
file_path = 'Histone_mutation.xlsx'
sheets = ['H2A', 'H2B', 'H3', 'H4']

dfs = []
for sheet in sheets:
    temp_df = pd.read_excel(file_path, sheet_name=sheet)
    temp_df['histone_type'] = sheet
    dfs.append(temp_df)

df = pd.concat(dfs, ignore_index=True)

# ===== 2. Extract numeric position from residue_real_site =====
df['site_num'] = df['residue_real_site'].astype(str).str.extract(r'(\d+)', expand=False)
df['site_num'] = pd.to_numeric(df['site_num'], errors='coerce')
df = df.dropna(subset=['site_num'])

# ===== 2.5 Reference sequence consistency check =====
def parse_residue_frequency(freq_str, target_residue):
    """Find target residue in 'C:1,G:4' format and return its frequency."""
    if pd.isna(freq_str) or pd.isna(target_residue):
        return False, 0.0
    pairs = re.findall(r'([A-Za-z])\s*:\s*(\d+(?:\.\d+)?)', str(freq_str))
    target = str(target_residue).upper()
    for res, val in pairs:
        if res.upper() == target:
            return True, float(val)
    return False, 0.0

df['real_residue'] = df['residue_real_site'].astype(str).str.extract(r'^([A-Za-z])', expand=False)

match_results = df.apply(
    lambda row: parse_residue_frequency(row['residue_frequency'], row['real_residue']),
    axis=1
)
df['is_ref_consistent'] = match_results.apply(lambda x: x[0])
df['matched_freq'] = match_results.apply(lambda x: x[1])

n_before = len(df)
df = df[df['is_ref_consistent']].copy()
print(f"Reference filter: {n_before} -> {len(df)} ({n_before - len(df)} removed)")

# ===== 3. Define core region ranges per histone =====
CORE_RANGES = {
    "H2A": (14, 118),
    "H2B": (24, None),
    "H3":  (37, None),
    "H4":  (21, None),
}

def in_core_region(row):
    h = row['histone_type']
    pos = row['site_num']
    if h not in CORE_RANGES:
        return False
    start, end = CORE_RANGES[h]
    if end is None:
        return pos >= start
    return start <= pos <= end

df['in_core'] = df.apply(in_core_region, axis=1)
core_df = df[df['in_core']].copy()

# ===== Output mutation counts per histone =====
histone_mutation_counts = core_df.groupby('histone_type')['matched_freq'].sum()
print(histone_mutation_counts)
print(f"Total: {histone_mutation_counts.sum():.0f}")

# ===== 4. Sum matched_freq per histone and plot =====
freq_sums = core_df.groupby('histone_type')['matched_freq'].sum()
total_freq = freq_sums.sum()

if total_freq == 0 or len(freq_sums) == 0:
    print("No valid sites in core regions after filtering.")
else:
    percentages = (freq_sums / total_freq * 100).round(1)
    labels = freq_sums.index.tolist()
    sizes = percentages.values

    colors = sns.color_palette("pastel", n_colors=len(labels))
    sns.set_style("whitegrid")
    plt.figure(figsize=(7, 7))

    wedges, texts, autotexts = plt.pie(
        sizes, labels=labels, colors=colors,
        autopct='%1.1f%%', startangle=140, pctdistance=0.8,
        wedgeprops={'linewidth': 2, 'edgecolor': 'black', 'width': 0.4},
        textprops={'fontsize': 20, 'weight': 'bold'}
    )
    for autotext in autotexts:
        autotext.set_size(20)
        autotext.set_weight('bold')
        autotext.set_color('black')

    plt.title(
        f"Core-region Mutations: {total_freq:.0f}",
        fontsize=24, pad=20, weight='bold'
    )
    plt.tight_layout()
    plt.savefig(
        'core_region_mutation_pie.tif',
        dpi=600, bbox_inches='tight', format='tiff',
        pil_kwargs={'compression': 'tiff_lzw'}
    )
    plt.show()