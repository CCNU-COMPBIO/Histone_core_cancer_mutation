#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Terminal DNA unwrapping -- contact-based criterion (step 1: per-frame analysis)

Writes to <dir>/unwrapping_contacts_out/<run>/:
  <end>_contacts.dat    frame time_ns contacts min_dist n_nt_lost n_bp_lost
  <end>_bp_matrix.dat   frame time_ns bp1..bp10   (1 = bp has no contact < CUT)
  <end>_bp_mindist.dat  frame time_ns bp1..bp10   (min heavy-atom distance per bp)

bp1 = outermost terminal bp.
Usage: python calc_unwrapping_contacts.py [--only H2AE56K:run2,H2AK74N:run1]
"""
import os
import glob
import argparse
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

warnings.filterwarnings('ignore')

# ---------------- Configuration ----------------
CUT           = 4.5          # heavy-atom contact cutoff (A)
FRAMES_PER_NS = 50           # 20 ps/frame
START_FRAME   = 100 * FRAMES_PER_NS + 1
END_FRAME     = 500 * FRAMES_PER_NS
OUT_DIRNAME   = 'unwrapping_contacts_out'
USE_PBC       = True

PROT_SEL = 'resid 1:753'
HEAVY = ' and not name H* and not name 1H* and not name 2H* and not name 3H*'

# Two 10-bp ends: s1 and s2 are reverse-paired
END_DEF = {
    'end1': dict(sel='resid 754:763 or resid 1036:1045',
                 s1=range(754, 764), s2=range(1045, 1035, -1)),
    'end2': dict(sel='resid 890:899 or resid 900:909',
                 s1=range(890, 900), s2=range(909, 899, -1)),
}
# bp1 = outermost terminal bp (end2 reversed)
BP_PAIRS = {}
for k, v in END_DEF.items():
    pairs = list(zip(v['s1'], v['s2']))
    BP_PAIRS[k] = pairs if k == 'end1' else pairs[::-1]

SYSTEMS = [
    ('WT',       'WT_rw',        'WT'),
    ('H2AR29Q',  'H2A_R29Q_rw',  'R29Q'),
    ('H2AE56K',  'H2A_E56K_rw',  'E56K'),
    ('H2AK74N',  'H2A_K74N_rw',  'K74N'),
    ('H2BE105K', 'H2B_E105K_rw', 'H2BE105K'),
    ('H2BF70L',  'H2B_F70L_rw',  'F70L'),
    ('H3E50K',   'H3_E50K_rw',   'E50K'),
    ('H3E73K',   'H3_E73K_rw',   'E73K'),
    ('H3E105K',  'H3_E105K_rw',  'E105K'),
]
RUNS = ['run1', 'run2', 'run3']


def find_traj(sys_dir, kw, run):
    fname = '%s_rw_%s_500ns.nc' % (kw, run)
    cands = [os.path.join(sys_dir, run, fname),
             os.path.join(sys_dir, '%s_traj' % kw, fname)]
    cands += sorted(glob.glob(os.path.join(sys_dir, '**', fname), recursive=True))
    for c in cands:
        if os.path.isfile(c):
            return c
    return None


def analyse(top, traj, end_name):
    u = mda.Universe(top, traj)
    prot = u.select_atoms(PROT_SEL + HEAVY)
    end = u.select_atoms(END_DEF[end_name]['sel'] + HEAVY)
    if len(prot) == 0 or len(end) == 0:
        raise ValueError('empty selection prot=%d end=%d' % (len(prot), len(end)))

    resids = end.atoms.resids
    res_pos = {r: np.where(resids == r)[0] for r in np.unique(resids)}
    pairs = BP_PAIRS[end_name]

    rows, bpmat, bpdist = [], [], []
    for i, ts in enumerate(u.trajectory[START_FRAME - 1:END_FRAME:FRAMES_PER_NS], start=1):
        box = u.dimensions if (USE_PBC and u.dimensions is not None) else None
        dm = distance_array(end.positions, prot.positions, box=box)
        atom_min = dm.min(axis=1)
        contacts = int((dm < CUT).sum())

        res_min = {r: float(atom_min[idx].min()) for r, idx in res_pos.items()}
        nt_lost = sum(1 for r in res_min if res_min[r] >= CUT)
        bp_mind = [min(res_min[r1], res_min[r2]) for r1, r2 in pairs]
        bp_lost_vec = [1 if x >= CUT else 0 for x in bp_mind]

        t = 100.02 + (i - 1) * 1.0
        rows.append([i, t, contacts, float(atom_min.min()), nt_lost, int(sum(bp_lost_vec))])
        bpmat.append([i, t] + bp_lost_vec)
        bpdist.append([i, t] + bp_mind)

        if i % 100 == 0:
            print('      ... %d/400 frames' % i, flush=True)

    return np.array(rows), np.array(bpmat), np.array(bpdist), len(prot), len(end)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default='')
    args = ap.parse_args()
    only = set(tuple(x.split(':')) for x in args.only.split(',') if x.strip())

    hdr_bp = ' '.join('bp%d' % (k + 1) for k in range(10))

    for label, sys_dir, kw in SYSTEMS:
        if not os.path.isdir(sys_dir):
            print('[skip] %s' % sys_dir)
            continue
        top = os.path.join(sys_dir, '%s_rw.prmtop' % kw)
        if not os.path.isfile(top):
            print('[skip] %s missing %s_rw.prmtop' % (sys_dir, kw))
            continue
        print('=== %s -> %s ===' % (label, sys_dir))

        for run in RUNS:
            if only and (label, run) not in only:
                continue
            traj = find_traj(sys_dir, kw, run)
            if traj is None:
                print('   [skip] %s %s trajectory not found' % (sys_dir, run))
                continue
            out_dir = os.path.join(sys_dir, OUT_DIRNAME, run)
            os.makedirs(out_dir, exist_ok=True)

            for end_name in END_DEF:
                try:
                    rows, bpmat, bpdist, nprot, nend = analyse(top, traj, end_name)
                except Exception as e:
                    print('   [warn] %s %s failed: %s' % (run, end_name, e))
                    continue

                np.savetxt(os.path.join(out_dir, '%s_contacts.dat' % end_name), rows,
                           header='frame time_ns contacts min_dist n_nt_lost n_bp_lost',
                           fmt='%d %.2f %d %.3f %d %d')
                np.savetxt(os.path.join(out_dir, '%s_bp_matrix.dat' % end_name), bpmat,
                           header='frame time_ns ' + hdr_bp,
                           fmt='%d %.2f ' + ' '.join(['%d'] * 10))
                np.savetxt(os.path.join(out_dir, '%s_bp_mindist.dat' % end_name), bpdist,
                           header='frame time_ns ' + hdr_bp,
                           fmt='%d %.2f ' + ' '.join(['%.3f'] * 10))

                nb = rows[:, 5]
                k = int(np.argmax(nb))
                print('   %s %s: %d frames (prot=%d, end=%d) | max bp_lost=%d @ %.1f ns '
                      '(minDist=%.1f A)'
                      % (run, end_name, len(rows), nprot, nend,
                         int(nb[k]), rows[k, 1], rows[k, 3]))

if __name__ == '__main__':
    main()