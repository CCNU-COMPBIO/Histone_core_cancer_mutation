import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ==========================================
# 1. Global style settings
# ==========================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.unicode_minus': False,
    'font.size': 18,
    'axes.labelsize': 20,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18
})

# ==========================================
# 2. Color mapping: directly by mutant abbreviation
# ==========================================
MUTANT_COLOR = {
    "R29Q":  "#D16BA5",   # H2A mutants
    "E56K":  "#D16BA5",
    "K74N":  "#D16BA5",
    "F70L":  "#86A8E7",   # H2B mutant
    "E50K":  "#ff9671",   # H3 mutants
    "E73K":  "#ff9671",
    "E105K": "#ff9671",
}

# ==========================================
# 3. Full mutant names
# ==========================================
MUTANT_FULL_NAME = {
    "R29Q":  "H2AR29Q",
    "E56K":  "H2AE56K",
    "K74N":  "H2AK74N",
    "F70L":  "H2BF70L",
    "E50K":  "H3E50K",
    "E73K":  "H3E73K",
    "E105K": "H3E105K",
}

# ==========================================
# 4. Plotting function (vertical grouped bars, each bar labeled)
# ==========================================
def plot_grouped_bars_vertical(df_plot, filename):
    """
    Required columns in df_plot: Group, Label, Type, Value, SEM, Color, Histone, Chain_Pair
    Only one WT per group (since WT is identical within a group).
    X‑axis tick labels show each bar's Label (WT or mutant full name).
    """
    # Fixed group order by Histone
    histone_order = {"H2A": 0, "H2B": 1, "H3": 2, "H4": 3}
    all_groups = df_plot['Group'].unique()
    def sort_key(g):
        parts = g.split(' ')
        his = parts[0] if len(parts) >= 1 else ''
        return histone_order.get(his, 99)
    groups = sorted(all_groups, key=sort_key)

    # Build final data: per group: 1 WT + N mutants
    final_records = []       # all bar records in plotting order
    group_sizes = []         # number of bars in each group
    group_start_x = []       # starting x coordinate for each group

    bar_width = 0.6
    gap_between_groups = 0.8
    current_x = 0.0

    for g in groups:
        g_data = df_plot[df_plot['Group'] == g].copy()
        wt_data = g_data[g_data['Type'] == 'WT']
        mut_data = g_data[g_data['Type'] == 'Mutant']

        if wt_data.empty:
            continue
        # Take the first WT record (all WT in the same group are identical)
        wt_record = wt_data.iloc[0].to_dict()
        mut_records = mut_data.to_dict('records')
        mut_records.sort(key=lambda x: x['Label'])  # keep consistent ordering

        group_records = [wt_record] + mut_records
        size = len(group_records)
        group_sizes.append(size)
        group_start_x.append(current_x)

        for rec in group_records:
            final_records.append(rec)

        current_x += size * bar_width + gap_between_groups

    # Generate x positions for each bar
    x_positions = []
    for idx, size in enumerate(group_sizes):
        start_x = group_start_x[idx]
        for i in range(size):
            x_positions.append(start_x + i * bar_width)

    # Extract plotting data
    values = [d['Value'] for d in final_records]
    sems = [d['SEM'] for d in final_records]
    colors = [d['Color'] for d in final_records]
    labels = [d['Label'] for d in final_records]   # label for each bar

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(16, 5))
    ax.bar(
        x_positions, values, yerr=sems, width=bar_width,
        color=colors, edgecolor='black', linewidth=1.5, alpha=0.7,
        error_kw={'ecolor': '#333333', 'elinewidth': 1.5, 'capsize': 4, 'capthick': 1.5},
        zorder=2
    )

    # ---- X axis: show each bar's label, rotate 45°, no xlabel ----
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=45, ha='right', va='top', fontsize=18)
    # No x‑axis title

    # ---- Y axis ----
    ax.set_ylabel(r"$\mathbf{S_{RMSF}}$ (Å)", fontweight='bold', fontsize=20)
    y_max = max(values) + max(sems) * 1.2 if values else 1
    ax.set_ylim(0, y_max)

    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_edgecolor('black')

    # ---- Legend: WT + each (Histone, Chain_Pair) combination ----
    pair_groups = df_plot[['Histone', 'Chain_Pair']].drop_duplicates()
    pair_groups['Histone_Order'] = pair_groups['Histone'].map(histone_order)
    pair_groups = pair_groups.sort_values(['Histone_Order', 'Chain_Pair'])

    legend_elements = [Patch(facecolor="#808080", edgecolor='black', label='WT', alpha=0.7)]

    for _, row in pair_groups.iterrows():
        his = row['Histone']
        pair = row['Chain_Pair']
        # Get the mutant color for this histone+pair (take first mutant)
        mask = (df_plot['Histone'] == his) & (df_plot['Chain_Pair'] == pair) & (df_plot['Type'] == 'Mutant')
        if mask.any():
            color = df_plot.loc[mask].iloc[0]['Color']
        else:
            color = "#808080"
        label = f"{his} ({pair})"
        legend_elements.append(
            Patch(facecolor=color, edgecolor='black', label=label, alpha=0.7)
        )

    ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5), frameon=False)

    # Adjust layout to make room for the legend
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ==========================================
# 5. Data loading (single‑pair S_RMSF output)
# ==========================================
csv_path = "SRMSF_SinglePair_Trim5_PerRunDetail.csv"   # confirm the filename
df = pd.read_csv(csv_path)

# Build plotting data: generate one WT record per group
records = []
# Get the unique WT values per (Histone, Chain_Mut, Chain_Copy)
wt_unique = df.groupby(['Histone', 'Chain_Mut', 'Chain_Copy'], as_index=False).first()
for _, row in wt_unique.iterrows():
    histone = row['Histone']
    chain_mut = row['Chain_Mut']
    chain_copy = row['Chain_Copy']
    chain_pair = f"{chain_mut} vs {chain_copy}"
    group = f"{histone} {chain_pair}"

    records.append({
        'Group': group,
        'Label': 'WT',
        'Type': 'WT',
        'Value': row['SRMSF_WT_Mean_A'],
        'SEM': row['SRMSF_WT_SEM_A'],
        'Color': '#808080',
        'Histone': histone,
        'Chain_Pair': chain_pair,
    })

# Add mutant records
for _, row in df.iterrows():
    histone = row['Histone']
    variant = row['Variant']
    chain_mut = row['Chain_Mut']
    chain_copy = row['Chain_Copy']
    chain_pair = f"{chain_mut} vs {chain_copy}"
    group = f"{histone} {chain_pair}"

    mut_name = MUTANT_FULL_NAME.get(variant, variant)
    color = MUTANT_COLOR.get(variant, '#000000')
    records.append({
        'Group': group,
        'Label': mut_name,
        'Type': 'Mutant',
        'Value': row['SRMSF_Mut_Mean_A'],
        'SEM': row['SRMSF_Mut_SEM_A'],
        'Color': color,
        'Histone': histone,
        'Chain_Pair': chain_pair,
    })

df_plot = pd.DataFrame(records)

# Sort groups (H2A, H2B, H3)
histone_order = {"H2A": 0, "H2B": 1, "H3": 2, "H4": 3}
df_plot['Group_Order'] = df_plot['Histone'].map(histone_order)
df_plot = df_plot.sort_values(['Group_Order', 'Chain_Pair', 'Type']).reset_index(drop=True)

# Generate the plot
plot_grouped_bars_vertical(df_plot, "SRMSF_single_pair_vertical.png")
print("\n[DONE] Plotting completed.")