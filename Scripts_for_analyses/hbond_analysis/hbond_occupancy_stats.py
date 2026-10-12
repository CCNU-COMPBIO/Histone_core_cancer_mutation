#!/usr/bin/env python3
"""
Automatically discover all chain pairs, count hydrogen bonds with occupancy >= 50%,
and compute mean and SEM across 3 runs.

Directory structure expected:
    E105K_rw/
        run1/
            E105K_rw_run1_chainA-chainB-detail.dat
            E105K_rw_run1_chainA-chainIJ-detail.dat
            ...
        run2/ ...
        run3/ ...
"""

import os
import re
import math
from collections import defaultdict

# ========== User parameters ==========
BASE_DIR = "E105K_rw"               # top-level directory
MUTATION = "E105K"                  # mutation name (used in file prefix)
RUNS = [1, 2, 3]
OCCUPANCY_THRESHOLD = 50.0          # threshold (%); count bonds with occupancy >= this
# ====================================

def parse_file(filepath, threshold):
    """Return the number of hydrogen bonds with occupancy >= threshold."""
    count = 0
    with open(filepath, 'r') as f:
        lines = f.readlines()

    data_start = False
    for line in lines:
        if line.strip().startswith("donor"):
            data_start = True
            continue
        if not data_start:
            continue
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        occupancy_str = parts[2]
        try:
            occupancy_val = float(occupancy_str.rstrip('%'))
        except ValueError:
            continue
        if occupancy_val >= threshold:
            count += 1
    return count

def chain_pair_to_short(chain_pair):
    """
    Convert 'chainA-chainB' -> 'AB', 'chainA-chainIJ' -> 'AIJ'.
    """
    parts = chain_pair.split('-')
    if len(parts) != 2:
        return chain_pair
    left = parts[0].replace('chain', '')
    right = parts[1].replace('chain', '')
    return left + right

def main():
    # Store data: { short_name: { run_number: count } }
    all_data = defaultdict(lambda: defaultdict(int))

    # Walk through each run directory and collect matching .dat files
    for run in RUNS:
        run_dir = os.path.join(BASE_DIR, f"run{run}")
        if not os.path.isdir(run_dir):
            print(f"Warning: directory not found, skipping: {run_dir}")
            continue

        for filename in os.listdir(run_dir):
            # Match pattern: {MUTATION}_rw_run{run}_chain??-chain??-detail.dat
            pattern = rf"^{re.escape(MUTATION)}_rw_run{run}_(chain\S+-chain\S+)-detail\.dat$"
            match = re.match(pattern, filename)
            if not match:
                continue

            chain_pair = match.group(1)
            short_name = chain_pair_to_short(chain_pair)

            filepath = os.path.join(run_dir, filename)
            count = parse_file(filepath, OCCUPANCY_THRESHOLD)
            all_data[short_name][run] = count

    # ---------- Output statistics ----------
    print(f"\n===== Hydrogen bond statistics (occupancy >= {OCCUPANCY_THRESHOLD}%) =====")
    if not all_data:
        print("No data files found.")
        return

    for short_name in sorted(all_data.keys()):
        run_counts = all_data[short_name]
        print(f"\nChain pair: {short_name}")

        run_values = []
        for run in RUNS:
            val = run_counts.get(run)
            if val is not None:
                print(f"  Run{run}: {val}")
                run_values.append(val)
            else:
                print(f"  Run{run}: N/A")

        if not run_values:
            print("  No valid data, skipping statistics.")
            continue

        n = len(run_values)
        mean_val = sum(run_values) / n
        if n >= 2:
            std_dev = math.sqrt(sum((x - mean_val) ** 2 for x in run_values) / (n - 1))
            sem = std_dev / math.sqrt(n)
        else:
            sem = 0.0

        print(f"  Mean: {mean_val:.2f}")
        print(f"  SEM:  {sem:.2f}")

if __name__ == "__main__":
    main()