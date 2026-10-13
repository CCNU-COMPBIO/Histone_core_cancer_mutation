import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# ===== Seaborn style and Arial font =====
sns.set_style("white")
sns.set_context("notebook", font_scale=1.2)
plt.rcParams['font.family'] = 'Arial'

# ===== Configuration =====
base_dir = "."
runs = [1, 2, 3]

mutants_info = {
    "H2A_R29Q":  ("C29", "H2AR29Q"),
    "H2A_E56K":  ("C56", "H2AE56K"),
    "H2A_K74N":  ("C74", "H2AK74N"),
    "H2B_F70L":  ("D70", "H2BF70L"),
    "H2B_E105K": ("D105", "H2BE105K"),
    "H3_E50K":   ("A50", "H3E50K"),
    "H3_E73K":   ("A73", "H3E73K"),
    "H3_E105K":  ("A105", "H3E105K"),
}

# Colors consistent with Figure a and Figure c
color_wt  = '#C2C2C2'   # WT: light gray
color_mut = '#6A95D6'   # Mutant: blue

# ===== Read a single file =====
def read_hbonds(filepath):
    if not os.path.isfile(filepath):
        print(f"   [Warning] File does not exist: {filepath}")
        return None, None
    data = np.loadtxt(filepath)
    if data.ndim == 1:
        frames = np.arange(len(data))
        hbonds = data
    else:
        frames = data[:, 0]
        hbonds = data[:, 1] if data.shape[1] >= 2 else data[:, 0]
    print(f"   Successfully read: {os.path.basename(filepath)} frames={len(hbonds)}")
    return frames, hbonds

# ===== Load three runs =====
def load_system(folder_name, file_prefix, resid):
    all_hbonds = []
    for run in runs:
        fname = f"{file_prefix}_run{run}_resid_{resid}-hbonds.dat"
        path = os.path.join(base_dir, folder_name, f"run{run}", fname)
        print(f"Trying to read: {path}")
        frames, hbonds = read_hbonds(path)
        if frames is not None:
            all_hbonds.append(hbonds)
            print(f"   Using all {len(hbonds)} frames")

    if not all_hbonds:
        print(f"   [Error] No valid data found for system {folder_name}")
        return None, None, None

    n_runs = len(all_hbonds)
    run_means = [np.mean(hb) for hb in all_hbonds]
    overall_mean = np.mean(run_means)
    sem = np.std(run_means, ddof=1) / np.sqrt(n_runs) if n_runs > 1 else 0.0
    print(f"   [Info] {folder_name}: {n_runs} runs kept, mean={overall_mean:.2f}, SEM={sem:.2f}")
    return overall_mean, sem, run_means

# ===== Main pipeline =====
group_data = []

for mutant, (resid, full_name) in mutants_info.items():
    print(f"\n===== Processing {mutant} (residue {resid}) -> {full_name} =====")
    mut_folder = f"{mutant}_rw"
    mut_file_prefix = mutant.split('_', 1)[1] + "_rw"

    wt_folder = "WT_rw"
    wt_file_prefix = "WT_rw"

    wt_mean, wt_sem, wt_runs = load_system(wt_folder, wt_file_prefix, resid)
    mut_mean, mut_sem, mut_runs = load_system(mut_folder, mut_file_prefix, resid)

    if wt_mean is None or mut_mean is None:
        print(f"   [Skip] Missing data for {mutant}")
        continue

    group_data.append({
        'group': mutant,
        'full_name': full_name,
        'wt_mean': wt_mean,
        'wt_sem': wt_sem,
        'wt_runs': wt_runs,
        'mut_mean': mut_mean,
        'mut_sem': mut_sem,
        'mut_runs': mut_runs,
    })

if not group_data:
    print("No valid data, exiting")
    exit()

# ===== Plotting =====
fig, ax = plt.subplots(figsize=(16, 5))

n_groups = len(group_data)
x = np.arange(n_groups)
width = 0.35

wt_means = [d['wt_mean'] for d in group_data]
wt_sems  = [d['wt_sem'] for d in group_data]
ax.bar(x - width/2, wt_means, yerr=wt_sems, width=width,
       color=color_wt, alpha=0.7, edgecolor='black',
       capsize=6, error_kw=dict(ecolor='black', lw=2, capsize=6, capthick=2),
       label='WT')

mut_means = [d['mut_mean'] for d in group_data]
mut_sems  = [d['mut_sem'] for d in group_data]
ax.bar(x + width/2, mut_means, yerr=mut_sems, width=width,
       color=color_mut, alpha=0.7, edgecolor='black',
       capsize=6, error_kw=dict(ecolor='black', lw=2, capsize=6, capthick=2),
       label='Mutant')

max_offset = width * 0.5 * 0.6

for i, d in enumerate(group_data):
    wt_vals = d['wt_runs']
    n_wt = len(wt_vals)
    if n_wt > 0:
        offsets_wt = np.linspace(-max_offset, max_offset, n_wt)
        ax.scatter(x[i] - width/2 + offsets_wt, wt_vals,
                   s=45, facecolors='white', edgecolors='black', linewidth=1.5,
                   zorder=10, label='_nolegend_')
    mut_vals = d['mut_runs']
    n_mut = len(mut_vals)
    if n_mut > 0:
        offsets_mut = np.linspace(-max_offset, max_offset, n_mut)
        ax.scatter(x[i] + width/2 + offsets_mut, mut_vals,
                   s=45, facecolors='white', edgecolors='black', linewidth=1.5,
                   zorder=10, label='_nolegend_')

ax.set_title("Hydrogen bonds formed by the mutated residue in WT and mutants",
             fontsize=22, weight='bold', pad=12)
ax.set_xticks(x)
ax.set_xticklabels([d['full_name'] for d in group_data], rotation=45, ha='right')
ax.set_ylabel('Number of H-bonds', fontweight='bold', fontsize=20)
ax.tick_params(axis='both', labelsize=18)

for spine in ax.spines.values():
    spine.set_linewidth(2)
    spine.set_color('black')

legend_elements = [
    Patch(facecolor=color_wt,  edgecolor='black', alpha=0.7, label='WT'),
    Patch(facecolor=color_mut, edgecolor='black', alpha=0.7, label='Mutant'),
]
ax.legend(handles=legend_elements, fontsize=18, loc='upper right', frameon=True)

ax.grid(False)

plt.tight_layout()
outfile = "WT_vs_all_mutants_hbonds.png"
plt.savefig(outfile, dpi=600, bbox_inches='tight')
print(f"\n[Saved] {outfile}")
plt.show()

# ===== Summary statistics =====
print("\nSummary statistics:")
print(f"{'Mutant':<10}{'WT_mean':>10}{'WT_SEM':>10}{'Mut_mean':>10}{'Mut_SEM':>10}")
for d in group_data:
    print(f"{d['group']:<10}{d['wt_mean']:>10.2f}{d['wt_sem']:>10.2f}"
          f"{d['mut_mean']:>10.2f}{d['mut_sem']:>10.2f}")