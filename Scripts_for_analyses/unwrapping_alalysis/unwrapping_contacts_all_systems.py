#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Terminal-DNA unwrapping: identify the trajectory segments excluded from the
analysis, using the histone-DNA contact criterion.

Criterion: at least BP_K of the ten terminal base pairs of a DNA end have no
histone heavy-atom contact within CUT A, for at least MIN_CONSEC consecutive
frames (1 frame = 1 ns). bp1 = outermost base pair.

If <end>_bp_mindist.dat is present, CUT can be varied freely (sensitivity
analysis, no trajectory re-reading); otherwise the precomputed n_bp_lost
obtained with a 4.5 A cut-off is used.

Output: excluded_segments_contacts.csv, one row per segment.
"""
import csv
import numpy as np
from pathlib import Path

BASE = Path('.')

CUT        = 4.5     # contact distance cut-off (A)
BP_K       = 8       # minimum number of contact-free terminal base pairs
MIN_CONSEC = 15      # minimum number of consecutive frames (ns)
OUT_CSV    = BASE / 'excluded_segments_contacts.csv'

SYSTEMS = [('WT', 'WT_rw'), ('H2AR29Q', 'H2A_R29Q_rw'), ('H2AE56K', 'H2A_E56K_rw'),
           ('H2AK74N', 'H2A_K74N_rw'), ('H2BE105K', 'H2B_E105K_rw'), ('H2BF70L', 'H2B_F70L_rw'),
           ('H3E50K', 'H3_E50K_rw'), ('H3E73K', 'H3_E73K_rw'), ('H3E105K', 'H3_E105K_rw')]
RUNS, ENDS = ['run1', 'run2', 'run3'], ['end1', 'end2']

series = {}
for label, sd in SYSTEMS:
    for run in RUNS:
        for end in ENDS:
            outc = BASE / sd / 'unwrapping_contacts_out' / run / ('%s_contacts.dat' % end)
            if not outc.exists():
                continue
            d = np.loadtxt(outc, comments='#')
            t, mind, nbp = d[:, 1], d[:, 3], d[:, 5]

            outm = BASE / sd / 'unwrapping_contacts_out' / run / ('%s_bp_mindist.dat' % end)
            if outm.exists():
                bpmin = np.loadtxt(outm, comments='#')[:, 2:12]
                flag = (bpmin >= CUT).sum(axis=1) >= BP_K
            else:
                flag = nbp >= BP_K
                if abs(CUT - 4.5) > 1e-9:
                    print('[warn] %s %s %s: bp_mindist.dat missing, CUT ignored'
                          % (label, run, end))

            series[(label, run, end)] = dict(t=t, mind=mind, flag=flag)

print('read %d series' % len(series))


def runs(flag, min_consec=MIN_CONSEC):
    segs, n, st = [], 0, None
    for i, a in enumerate(flag):
        if a:
            if n == 0:
                st = i
            n += 1
        else:
            if n >= min_consec:
                segs.append((st, i - 1, n))
            n = 0
    if n >= min_consec:
        segs.append((st, len(flag) - 1, n))
    return segs


rows = []
for key in sorted(series):
    s = series[key]
    for st, en, n in runs(s['flag']):
        rows.append(dict(system=key[0], run=key[1], end=key[2],
                         start_ns=round(float(s['t'][st]), 1),
                         end_ns=round(float(s['t'][en]), 1),
                         n_frames=n,
                         min_dist_min=round(float(s['mind'][st:en + 1].min()), 2),
                         min_dist_max=round(float(s['mind'][st:en + 1].max()), 2)))

fields = ['system', 'run', 'end', 'start_ns', 'end_ns', 'n_frames',
          'min_dist_min', 'min_dist_max']
with open(OUT_CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print('criterion: bp_open>=%d, contact < %.1f A, >= %d consecutive frames'
      % (BP_K, CUT, MIN_CONSEC))
print('segments: %d' % len(rows))
print('%-9s %-5s %-5s | %-40s | %6s | %s'
      % ('system', 'run', 'end', 'segment (ns)', 'frames', 'min distance (A)'))
print('-' * 100)
for key in sorted({(r['system'], r['run'], r['end']) for r in rows}):
    g = [r for r in rows if (r['system'], r['run'], r['end']) == key]
    iv = ', '.join('%g-%g' % (r['start_ns'], r['end_ns']) for r in g)
    tot = sum(r['n_frames'] for r in g)
    dmin = min(r['min_dist_min'] for r in g)
    dmax = max(r['min_dist_max'] for r in g)
    print('%-9s %-5s %-5s | %-40s | %6d | %.1f-%.1f'
          % (key[0], key[1], key[2], iv, tot, dmin, dmax))