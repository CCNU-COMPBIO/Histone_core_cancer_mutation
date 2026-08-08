#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Heavy-atom contact counts between H3(A) and H4(B): WT vs mutants
Run: cd /media/lenovo/MyDrive/rm_wat && python contacts_AB.py
Time window: 0–500 ns (trajectory starts at 6 ns, shifted by -6 ns);
only runs 1/2/3 used; short runs discarded
"""

import os, glob, warnings
import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

warnings.filterwarnings('ignore', category=UserWarning)

# ================= Configuration =================
H3_SEL = 'resid 1:97 or resid 377:475 and not name H*'
H4_SEL = 'resid 98:175 or resid 476:557 and not name H*'
CUT = 4.5
FRAME_STEP = 50
T_START_NS = 100.0   # shifted time (original 106 ns)
T_END_NS = 500.0     # shifted time (original 506 ns)
MIN_FRAMES = 350     # drop run if effective frames < this
RUN_NAMES = ['run1', 'run2', 'run3']
OUT_DIR = 'contacts_out'

SYSTEMS = [
    ('WT', 'WT', '#000000'),
    ('H3E50K', 'E50K', '#e41a1c'),
    ('H3E73K', 'E73K', '#377eb8'),
    ('H3E105K', 'E105K', '#4daf4a'),
]

# ================= File locating =================
def find_system_dir(keyword):
    dirs = [d for d in sorted(os.listdir('.')) if os.path.isdir(d)]
    for cand in (f'{keyword}_rw', keyword):
        if cand in dirs:
            return cand
    hits = [d for d in dirs if keyword.lower() in d.lower()]
    if len(hits) > 1:
        print(f'  [warn] {keyword} matches multiple {hits}, using {hits[0]}')
    return hits[0] if hits else None


def find_topology(sys_dir, keyword):
    exact = os.path.join(sys_dir, f'{keyword}_rw.prmtop')
    if os.path.isfile(exact):
        return exact
    prmtops = glob.glob(os.path.join(sys_dir, '*.prmtop'))
    if len(prmtops) == 1:
        return prmtops[0]
    for p in prmtops:
        if keyword.lower() in os.path.basename(p).lower():
            return p
    return prmtops[0] if prmtops else None


# ================= Single run: exact shifted time window =================
def compute_one_run(top, traj):
    u = mda.Universe(top, traj)
    h3 = u.select_atoms(H3_SEL)
    h4 = u.select_atoms(H4_SEL)
    if len(h3) == 0 or len(h4) == 0:
        raise ValueError(f'Selection empty: H3={len(h3)} H4={len(h4)}')

    times, counts = [], []
    if len(u.trajectory) == 0:
        raise ValueError('Trajectory is empty')

    first_ts = u.trajectory[0].time / 1000.0
    last_ts = u.trajectory[-1].time / 1000.0
    print(f'  [debug] trajectory time range: {first_ts:.1f}–{last_ts:.1f} ns')

    for ts in u.trajectory[::FRAME_STEP]:
        t_ns = ts.time / 1000.0
        t = t_ns - 6.0   # shift: original time - 6 ns
        if t < T_START_NS:
            continue
        if t > T_END_NS:
            break
        dm = mda.lib.distances.distance_array(h3.positions, h4.positions)
        times.append(t)
        counts.append(float((dm < CUT).sum()))

    if not times:
        raise ValueError(f'No frames in {T_START_NS}–{T_END_NS} ns')
    return np.asarray(times), np.asarray(counts)


# ================= Single system analysis =================
def analyze_system(label, sys_dir, keyword):
    top = find_topology(sys_dir, keyword)
    if top is None:
        print(f'  [error] no topology')
        return None
    print(f'  topo: {os.path.basename(top)}')

    valid = []
    for rn in RUN_NAMES:
        rd = os.path.join(sys_dir, rn)
        if not os.path.isdir(rd):
            print(f'  [skip] {rn} does not exist')
            continue
        trajs = sorted(glob.glob(os.path.join(rd, '*.nc')) +
                       glob.glob(os.path.join(rd, '*.ncdf')))
        if not trajs:
            print(f'  [skip] {rn} no trajectory')
            continue
        try:
            t, y = compute_one_run(top, trajs[0])
            if len(y) < MIN_FRAMES:
                print(f'  [drop] {rn}: {len(y)} frames (<{MIN_FRAMES}), discarding')
                continue
            valid.append((rn, t, y))
            print(f'  {rn}: {len(y)} frames, '
                  f'{t[0]:.1f}–{t[-1]:.1f} ns, mean={y.mean():.0f}')
        except Exception as e:
            print(f'  [warn] {rn} failed: {e}')

    if not valid:
        return None
    if len(valid) < 2:
        print(f'  [warn] only {len(valid)} valid runs, SEM=0')

    t_lengths = [len(t) for _, t, _ in valid]
    if len(set(t_lengths)) > 1:
        print(f'  [warn] time axis lengths differ: {t_lengths}')
        min_len = min(t_lengths)
        print(f'  [info] truncating to shortest length {min_len}')
        mat = np.vstack([y[:min_len] for _, _, y in valid])
        t_ref = valid[0][1][:min_len]
    else:
        L = min(len(y) for _, _, y in valid)
        mat = np.vstack([y[:L] for _, _, y in valid])
        t_ref = valid[0][1][:L]

    mean = mat.mean(axis=0)
    sem = (mat.std(axis=0, ddof=1) / np.sqrt(mat.shape[0])
           if mat.shape[0] > 1 else np.zeros_like(mean))
    run_means = mat.mean(axis=1)
    return {'t': t_ref, 'mean': mean, 'sem': sem,
            'n_runs': mat.shape[0],
            'overall_mean': run_means.mean(),
            'overall_sem': (run_means.std(ddof=1) / np.sqrt(len(run_means))
                            if len(run_means) > 1 else 0.0)}


# ================= Main workflow =================
def main():
    # ========== Font and size settings ==========
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.labelsize'] = 16
    plt.rcParams['axes.titlesize'] = 18
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 14
    # ============================================

    os.makedirs(OUT_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 5))
    summary = []

    for label, kw, color in SYSTEMS:
        sys_dir = find_system_dir(kw)
        if sys_dir is None:
            print(f'[error] System {label} ({kw}) not found')
            continue
        print(f'=== {label} -> {sys_dir} ===')
        res = analyze_system(label, sys_dir, kw)
        if res is None:
            continue

        np.savetxt(os.path.join(OUT_DIR, f'{label}_mean_sem.csv'),
                   np.column_stack([res['t'], res['mean'], res['sem']]),
                   header='time_ns  mean_contacts  sem', fmt='%.3f')

        sm = uniform_filter1d(res['mean'], size=10)
        sm_s = uniform_filter1d(res['sem'], size=10)
        ax.fill_between(res['t'], sm - sm_s, sm + sm_s,
                        color=color, alpha=0.18, linewidth=0)
        ax.plot(res['t'], sm, lw=1.6, color=color, label=label)
        summary.append((label, res['n_runs'],
                        res['overall_mean'], res['overall_sem']))

    if not summary:
        print('[error] No system succeeded')
        return

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    ax.set_xlabel('Time (ns)', fontweight='bold')
    ax.set_ylabel(f'Heavy-atom contacts (< {CUT} Å)', fontweight='bold')
    ax.legend(ncol=2, loc='lower right')
    plt.tight_layout()
    out_png = os.path.join(OUT_DIR, 'contacts_AB.png')
    plt.savefig(out_png, dpi=600)
    print(f'\nSaved: {out_png}')

    print(f"\n{'System':<10}{'n_runs':>7}{'Mean±SEM':>15}")
    for label, n, m, s in summary:
        print(f'{label:<10}{n:>7}{m:>10.0f} ± {s:4.0f}')


if __name__ == '__main__':
    main()