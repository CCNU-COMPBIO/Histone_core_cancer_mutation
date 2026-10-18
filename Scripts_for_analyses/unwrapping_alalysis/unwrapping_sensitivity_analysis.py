#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sensitivity of the excluded-segment set to the histone-DNA contact cut-off.

The contact-free base-pair threshold (8) and the persistence threshold (15 ns)
are held fixed at the values used throughout this work, so that only one
parameter is varied at a time.

Uses the precomputed per-base-pair minimum distances, so no trajectory is
re-read. Output: sensitivity_contact_cutoff.csv plus a table on stdout.
"""
import csv
import numpy as np
from pathlib import Path

BASE = Path('.')
SYSTEMS = [('WT', 'WT_rw'), ('H2AR29Q', 'H2A_R29Q_rw'), ('H2AE56K', 'H2A_E56K_rw'),
           ('H2AK74N', 'H2A_K74N_rw'), ('H2BE105K', 'H2B_E105K_rw'), ('H2BF70L', 'H2B_F70L_rw'),
           ('H3E50K', 'H3_E50K_rw'), ('H3E73K', 'H3_E73K_rw'), ('H3E105K', 'H3_E105K_rw')]
RUNS, ENDS = ['run1', 'run2', 'run3'], ['end1', 'end2']

BP_K       = 8                    # fixed
MIN_CONSEC = 15                   # fixed (ns)
CUTS       = [4.0, 4.5, 5.0]      # the only parameter varied
OUT_CSV    = BASE / 'sensitivity_contact_cutoff.csv'

data = {}
for label, sd in SYSTEMS:
    for run in RUNS:
        for end in ENDS:
            m = BASE / sd / 'unwrapping_contacts_out' / run / ('%s_bp_mindist.dat' % end)
            c = BASE / sd / 'unwrapping_contacts_out' / run / ('%s_contacts.dat' % end)
            if m.exists() and c.exists():
                D = np.loadtxt(m, comments='#')
                C = np.loadtxt(c, comments='#')
                data[(label, run, end)] = dict(t=D[:, 1], bpmin=D[:, 2:12], mind=C[:, 3])
print('read %d series' % len(data))


def segments(rec, cut):
    flag = (rec['bpmin'] >= cut).sum(axis=1) >= BP_K
    out, n, st = [], 0, None
    for i, a in enumerate(flag):
        if a:
            if n == 0:
                st = i
            n += 1
        else:
            if n >= MIN_CONSEC:
                out.append((st, i - 1))
            n = 0
    if n >= MIN_CONSEC:
        out.append((st, len(flag) - 1))
    return out


res = {}
for cut in CUTS:
    for key, rec in data.items():
        segs = segments(rec, cut)
        if segs:
            res[(cut, key)] = [(round(float(rec['t'][a]), 1), round(float(rec['t'][b]), 1),
                                b - a + 1) for a, b in segs]

keys = sorted({k for (_, k) in res})
rows = []
for cut in CUTS:
    for key in keys:
        segs = res.get((cut, key), [])
        rows.append(dict(cutoff_A=cut, system=key[0], run=key[1], end=key[2],
                         n_segments=len(segs),
                         segments='; '.join('%g-%g (%d)' % s for s in segs)))
with open(OUT_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['cutoff_A', 'system', 'run', 'end',
                                      'n_segments', 'segments'])
    w.writeheader()
    w.writerows(rows)
print('wrote %s' % OUT_CSV.name)

# ---------------- table ----------------
colw = 30
hdr = '%-11s|' % 'cut-off (A)'
for key in keys:
    hdr += (' %s %s' % (key[0], key[1])).ljust(colw) + '|'
print('\n' + hdr)
print('-' * len(hdr))
for cut in CUTS:
    line = '%-11g|' % cut
    for key in keys:
        segs = res.get((cut, key), [])
        cell = '; '.join('%g-%g' % (s[0], s[1]) for s in segs) if segs else 'not identified'
        line += ' ' + cell.ljust(colw - 1) + '|'
    print(line)