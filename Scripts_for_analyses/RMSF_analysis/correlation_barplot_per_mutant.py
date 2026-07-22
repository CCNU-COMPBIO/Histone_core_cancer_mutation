import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. Global style settings
# ==========================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'axes.unicode_minus': False,
    'font.size': 18,
    'axes.labelsize': 18,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 18
})

palette = sns.color_palette("tab10")[:2]  # two colours: WT and Mutant


# ==========================================
# 2. Plotting function: grouped horizontal bar chart with spacing between groups
# ==========================================
def plot_grouped_bars(df_plot, corr_type, filename):
    """
    df_plot must contain columns: 'Label' (display label), 'Value' (correlation coefficient),
    'SEM', 'Type' ('WT' or 'Mutant'), 'Group' (histone group).
    Groups are plotted with WT first then Mutant inside each group, with a blank row between groups.
    """
    # Sort by Group, then by Type (WT first), and optionally by Label
    df_plot = df_plot.sort_values(['Group', 'Type', 'Label'], ascending=[True, False, True])
    groups = df_plot['Group'].unique()

    y_positions = []
    y_labels = []
    current_y = 0
    group_start_indices = []  # starting index of each group
    for g in groups:
        group_data = df_plot[df_plot['Group'] == g]
        n_items = len(group_data)
        y_positions.extend(list(range(current_y, current_y + n_items)))
        y_labels.extend(group_data['Label'].tolist())
        group_start_indices.append(current_y)
        current_y += n_items + 1  # one blank row between groups

    values = df_plot['Value'].tolist()
    sems = df_plot['SEM'].tolist()
    types = df_plot['Type'].tolist()

    color_map = {'WT': palette[0], 'Mutant': palette[1]}

    fig, ax = plt.subplots(figsize=(7, 7))

    for i, (y, val, sem, typ) in enumerate(zip(y_positions, values, sems, types)):
        ax.barh(y, val, xerr=sem if not np.isnan(sem) else None,
                color=color_map[typ], edgecolor='black', linewidth=1.0,
                capsize=4, error_kw={'ecolor': '#333333', 'linewidth': 1.5, 'capthick': 1.5})

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels)

    # Draw a grey horizontal line between groups to separate them
    for start_idx in group_start_indices[1:]:
        ax.axhline(y=start_idx - 1, color='gray', linestyle='-', linewidth=1, alpha=1)

    ax.set_xlabel(f"{corr_type} Correlation Coefficient", fontweight='bold')
    ax.set_xlim(0, 1.05)

    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_edgecolor('black')

    # Legend (commented out; can be added if needed)
    # from matplotlib.patches import Patch
    # legend_elements = [Patch(facecolor=palette[0], edgecolor='black', label='WT'),
    #                    Patch(facecolor=palette[1], edgecolor='black', label='Mutant')]
    # ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)

    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close(fig)


# ==========================================
# 3. Main workflow
# ==========================================
csv_path = "RMSF_correlation_symmetry_pairs.csv"
df_raw = pd.read_csv(csv_path)

# Ensure there is a WT row per histone/chain pair
has_wt = df_raw['Mutant'].astype(str).str.upper().eq('WT').any()
if not has_wt:
    print("[INFO] No WT row detected in CSV; constructing from existing data...")
    wt_key_cols = ['Histone', 'Chain_Pair']
    df_wt_dedup = df_raw.drop_duplicates(subset=wt_key_cols, keep='first').copy()

    wt_rows = []
    for _, row in df_wt_dedup.iterrows():
        wt_row = {'Mutant': 'WT'}
        for kc in wt_key_cols:
            wt_row[kc] = row[kc]

        for c in df_raw.columns:
            if c in wt_key_cols or c == 'Mutant':
                continue
            if c.startswith('Symmetry_WT_'):
                wt_row[c] = row[c]
            else:
                wt_row[c] = np.nan

        wt_rows.append(wt_row)

    if wt_rows:
        df_wt = pd.DataFrame(wt_rows)
        df_raw = pd.concat([df_wt, df_raw], ignore_index=True)
        print(f"[OK] Added {len(wt_rows)} WT row(s) based on each histone/chain pair.")
    else:
        print("[WARN] Could not construct WT rows.")

is_wt = df_raw['Mutant'].astype(str).str.upper() == 'WT'

# Build plotting data
records = []
groups = df_raw.loc[~is_wt, 'Histone'].unique()
for g in groups:
    wt_mask = is_wt & (df_raw['Histone'] == g)
    if wt_mask.sum() == 0:
        print(f"[WARN] No WT data for histone {g}, skipping")
        continue
    wt_row = df_raw.loc[wt_mask].iloc[0]
    wt_pearson = wt_row['Symmetry_WT_Pearson']
    wt_pearson_sem = wt_row['Symmetry_WT_Pearson_SEM']
    wt_spearman = wt_row['Symmetry_WT_Spearman']
    wt_spearman_sem = wt_row['Symmetry_WT_Spearman_SEM']

    # WT entry
    records.append({
        'Label': f'{g} WT',
        'Group': g,
        'Type': 'WT',
        'Pearson_Value': wt_pearson,
        'Pearson_SEM': wt_pearson_sem,
        'Spearman_Value': wt_spearman,
        'Spearman_SEM': wt_spearman_sem,
    })

    # Mutant entries
    mut_mask = (~is_wt) & (df_raw['Histone'] == g)
    for _, mut_row in df_raw.loc[mut_mask].iterrows():
        mutant = mut_row['Mutant']
        records.append({
            'Label': f'{g}{mutant}',
            'Group': g,
            'Type': 'Mutant',
            'Pearson_Value': mut_row['Symmetry_Mutant_Pearson'],
            'Pearson_SEM': mut_row['Symmetry_Mutant_Pearson_SEM'],
            'Spearman_Value': mut_row['Symmetry_Mutant_Spearman'],
            'Spearman_SEM': mut_row['Symmetry_Mutant_Spearman_SEM'],
        })

df_plot_all = pd.DataFrame(records)
if df_plot_all.empty:
    raise ValueError("No data available for plotting")

# --- Pearson plot ---
df_pearson = df_plot_all[['Label', 'Group', 'Type', 'Pearson_Value', 'Pearson_SEM']].rename(
    columns={'Pearson_Value': 'Value', 'Pearson_SEM': 'SEM'}
)
plot_grouped_bars(df_pearson, 'Pearson',
                  filename='Pearson_correlations_grouped.png')

# --- Spearman plot ---
df_spearman = df_plot_all[['Label', 'Group', 'Type', 'Spearman_Value', 'Spearman_SEM']].rename(
    columns={'Spearman_Value': 'Value', 'Spearman_SEM': 'SEM'}
)
plot_grouped_bars(df_spearman, 'Spearman',
                  filename='Spearman_correlations_grouped.png')

print("\n[DONE] Generated two grouped bar plots with group spacing.")