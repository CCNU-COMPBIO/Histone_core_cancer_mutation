import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

palette = sns.color_palette("tab10")[:2]

# ==========================================
# 2. Single-mutant plotting function (no partner chain)
# ==========================================
def plot_single_mutant(df_plot, y_order, title, filename):
    df_long = df_plot.melt(
        id_vars=['Label'],
        value_vars=['Pearson', 'Spearman'],
        var_name='Correlation Type',
        value_name='Coefficient'
    )

    fig, ax = plt.subplots(figsize=(8, 4))

    sns.barplot(
        data=df_long,
        y="Label", x="Coefficient",
        hue="Correlation Type",
        order=y_order,
        palette=palette,
        orient="h",
        ax=ax,
        gap=0,
        width=0.6,
        edgecolor="black", linewidth=1.0
    )

    # SEM error bars
    n_rows = len(y_order)
    sem_map = {
        'Pearson': df_plot.set_index('Label')['Pearson_SEM'].to_dict(),
        'Spearman': df_plot.set_index('Label')['Spearman_SEM'].to_dict()
    }
    hue_order = ['Pearson', 'Spearman']
    for i, patch in enumerate(ax.patches):
        hue_idx = i % len(hue_order)
        row_idx = i // len(hue_order)
        if row_idx >= n_rows:
            continue
        label = y_order[row_idx]
        corr_type = hue_order[hue_idx]
        sem_val = sem_map[corr_type].get(label, 0)
        if not np.isnan(sem_val) and sem_val > 0:
            width = patch.get_width()
            y_center = patch.get_y() + patch.get_height() / 2
            ax.hlines(y=y_center, xmin=width - sem_val, xmax=width + sem_val,
                      colors='#333333', linewidth=1.5, zorder=5)
            cap_half = 0.012
            for xv in [width - sem_val, width + sem_val]:
                ax.vlines(x=xv, ymin=y_center - cap_half, ymax=y_center + cap_half,
                          colors='#333333', linewidth=1.5, zorder=5)

    # Style
    for spine in ax.spines.values():
        spine.set_linewidth(2.0)
        spine.set_edgecolor('black')

    ax.set_ylabel("")
    ax.set_xlabel("Correlation Coefficient", fontweight='bold')
    ax.set_xlim(0, 1.05)
    ax.margins(y=0.15)

    ax.set_title(title, fontsize=20, fontweight='bold', pad=12)

    legend = ax.legend(title="", bbox_to_anchor=(1.02, 1), loc="upper left",
                       frameon=False, borderaxespad=0)
    for p in legend.get_patches():
        p.set_edgecolor('black')
        p.set_linewidth(1.0)

    plt.tight_layout()
    plt.savefig(filename, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close(fig)


# ==========================================
# 3. Main workflow
# ==========================================
csv_path = "RMSF_symmetry_correlation.csv"
df_raw = pd.read_csv(csv_path)

# --- Add WT row per system if missing (using new column names) ---
has_wt = df_raw['Mutant'].astype(str).str.upper().eq('WT').any()
if not has_wt:
    print("[INFO] No WT row detected in CSV, constructing per system...")
    wt_key_cols = ['Histone', 'Chain_Pair']
    df_wt_dedup = df_raw.drop_duplicates(subset=wt_key_cols, keep='first').copy()

    wt_rows = []
    for _, row in df_wt_dedup.iterrows():
        wt_row = {'Mutant': 'WT'}
        for kc in wt_key_cols:
            wt_row[kc] = row[kc]

        # Copy only columns starting with 'Symmetry_WT_' (originally Cross_WT)
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
        print(f"[OK] Added {len(wt_rows)} WT records per system")
        print("\n[DEBUG] WT Pearson values per system:")
        for _, r in df_wt.iterrows():
            cross_val = r.get('Symmetry_WT_Pearson', np.nan)
            print(f"  {r['Histone']} {r['Chain_Pair']}: Cross={cross_val:.4f}")
    else:
        print("[WARN] Could not build WT rows")

is_wt = df_raw['Mutant'].astype(str).str.upper() == 'WT'

# --- Get all mutants ---
mutants = df_raw.loc[~is_wt, 'Mutant'].unique().tolist()
print(f"\n[INFO] Found {len(mutants)} mutants: {mutants}\n")

# --- Plot per mutant ---
for mutant in mutants:
    mut_mask = (~is_wt) & (df_raw['Mutant'] == mutant)
    if mut_mask.sum() == 0:
        print(f"[SKIP] {mutant}: no mutant data")
        continue

    histone = df_raw.loc[mut_mask, 'Histone'].iloc[0]
    pair_raw = df_raw.loc[mut_mask, 'Chain_Pair'].iloc[0]
    pair = str(pair_raw).replace('↔', ' vs ')

    # Exact WT match
    precise_wt_mask = is_wt & (df_raw['Histone'] == histone) & (df_raw['Chain_Pair'] == pair_raw)
    if precise_wt_mask.sum() != 1:
        print(f"[WARN] {mutant}: exact WT match count={precise_wt_mask.sum()}, skipping!")
        continue

    # Build plotting data (Label simplified, no parentheses)
    pair_display = str(pair_raw).replace('↔', ' vs ')

    records = []
    # WT
    records.append({
        'Label': 'WT',
        'Pearson': df_raw.loc[precise_wt_mask, 'Symmetry_WT_Pearson'].values[0],
        'Spearman': df_raw.loc[precise_wt_mask, 'Symmetry_WT_Spearman'].values[0],
        'Pearson_SEM': df_raw.loc[precise_wt_mask, 'Symmetry_WT_Pearson_SEM'].values[0],
        'Spearman_SEM': df_raw.loc[precise_wt_mask, 'Symmetry_WT_Spearman_SEM'].values[0],
    })
    # Mutant
    records.append({
        'Label': f"{histone}{mutant}",
        'Pearson': df_raw.loc[mut_mask, 'Symmetry_Mutant_Pearson'].values[0],
        'Spearman': df_raw.loc[mut_mask, 'Symmetry_Mutant_Spearman'].values[0],
        'Pearson_SEM': df_raw.loc[mut_mask, 'Symmetry_Mutant_Pearson_SEM'].values[0],
        'Spearman_SEM': df_raw.loc[mut_mask, 'Symmetry_Mutant_Spearman_SEM'].values[0],
    })

    df_plot = pd.DataFrame(records)

    # New title format
    title = f"{histone}{mutant} | {pair_display}"
    filename = f"rmsf_correlation_{mutant}.png"

    y_order = df_plot['Label'].tolist()[::-1]

    print(f">>> Plotting {title} ...")
    plot_single_mutant(df_plot, y_order, title, filename)
    print(f"    [OK] Saved: {filename}")

print("\n[DONE] All per-mutant plots generated.")
