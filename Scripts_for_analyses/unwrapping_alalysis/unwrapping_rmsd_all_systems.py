#!/usr/bin/env python3
import csv
import numpy as np
from pathlib import Path

# ========== Configuration ==========
BASE_DIR = Path(".").resolve()

# (display label, directory name)
SYSTEMS = [
    ("H2AR29Q",  "R29Q_rw"),
    ("H2AE56K",  "E56K_rw"),
    ("H2AK74N",  "K74N_rw"),
    ("H2BE105K",  "H2BE105K_rw"),
    ("H2BF70L",  "F70L_rw"),
    ("H3E50K",   "E50K_rw"),
    ("H3E73K",   "E73K_rw"),
    ("H3E105K",  "E105K_rw"),
    ("WT",       "WT_rw"),
]

RUNS = ["run1", "run2", "run3"]
ENDS = ["end1", "end2"]
THRESHOLDS = [10.0, 12.0]
MIN_CONSEC = 15          # consecutive frames; each frame is 1 ns

OUT_CSV = BASE_DIR / "unwrapping_all_systems.csv"

# ========== Utility functions ==========
def frame_to_ns(frame):
    """Output frame 1 -> 100.02 ns; each subsequent frame adds 1 ns."""
    return 100.02 + (frame - 1) * 1.0

def read_rmsd(path):
    data = np.loadtxt(path, comments="#")
    frames = data[:, 0]
    rmsd = data[:, 1]
    times = np.array([frame_to_ns(f) for f in frames])
    return times, rmsd

def find_persistent_segments(times, rmsd, threshold, min_consec=MIN_CONSEC):
    """Segments where RMSD > threshold for >= min_consec consecutive frames: [(start_ns, end_ns, n_frames), ...]"""
    above = rmsd > threshold
    segments = []
    count = 0
    start_idx = None
    for i, a in enumerate(above):
        if a:
            if count == 0:
                start_idx = i
            count += 1
        else:
            if count >= min_consec:
                segments.append((times[start_idx], times[i - 1], count))
            count = 0
            start_idx = None
    if count >= min_consec:
        segments.append((times[start_idx], times[-1], count))
    return segments

# ========== Main workflow ==========
rows = []
missing = []

for label, sys_dirname in SYSTEMS:
    for run in RUNS:
        for end in ENDS:
            dat_path = BASE_DIR / sys_dirname / "unwrapping_out" / run / f"{end}_rmsd.dat"
            if not dat_path.exists():
                missing.append(str(dat_path))
                continue

            times, rmsd = read_rmsd(dat_path)
            total_frames = len(rmsd)

            for thr in THRESHOLDS:
                segs = find_persistent_segments(times, rmsd, thr)
                total_unwrapped = sum(s[2] for s in segs)
                frac = total_unwrapped / total_frames if total_frames else 0.0

                if not segs:
                    rows.append({
                        "system": label,
                        "run": run,
                        "end": end,
                        "threshold_A": f"{thr:.1f}",
                        "n_segments": 0,
                        "segment_start_ns": "",
                        "segment_end_ns": "",
                        "segment_duration_ns": "",
                        "segment_n_frames": "",
                        "total_unwrapped_frames": 0,
                        "total_unwrapped_fraction": f"{frac:.3f}",
                    })
                else:
                    for s_start, s_end, n in segs:
                        rows.append({
                            "system": label,
                            "run": run,
                            "end": end,
                            "threshold_A": f"{thr:.1f}",
                            "n_segments": len(segs),
                            "segment_start_ns": f"{s_start:.1f}",
                            "segment_end_ns": f"{s_end:.1f}",
                            "segment_duration_ns": f"{s_end - s_start:.1f}",
                            "segment_n_frames": n,
                            "total_unwrapped_frames": total_unwrapped,
                            "total_unwrapped_fraction": f"{frac:.3f}",
                        })

# ========== Write CSV ==========
fieldnames = [
    "system", "run", "end", "threshold_A", "n_segments",
    "segment_start_ns", "segment_end_ns", "segment_duration_ns", "segment_n_frames",
    "total_unwrapped_frames", "total_unwrapped_fraction",
]

with open(OUT_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Saved: {OUT_CSV}")
print(f"Total rows: {len(rows)}")
if missing:
    print(f"\nMissing {len(missing)} dat files:")
    for m in missing:
        print(f"  {m}")