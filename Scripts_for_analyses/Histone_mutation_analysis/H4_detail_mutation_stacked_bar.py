#!/usr/bin/env python3
# encoding: utf-8
"""
Generate stacked bar plot for H4 histone mutations (frequency >= 5).
Reads Histone_mutation.xlsx, sheet "H4".
Only keeps core-region sites with total mutation count >= 5.
Matches mutations by original amino acid only (ignores site number in mutation token).
Output: plots_mutations_H4/residue_detail_mutations_H4_freq5.tif
"""

import re
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# ================== CONFIG ==================
INPUT_FILE = "Histone_mutation.xlsx"
SHEET_NAME = "H4"
OUT_DIR = "plots_mutations"
os.makedirs(OUT_DIR, exist_ok=True)

# Core region for H4: start=21, end=None (no upper bound)
CORE_RANGES = {"H4": (21, None)}
RESIDUE_COL = "residue_real_site"
MUT_COL_KEYWORDS = ["replication-dependent", "replication-independent"]

# Charge sets
POS = {"K", "R"}
NEG = {"D", "E"}

# Colors for destination residues
POS_COLOR = "#1f77b4"   # blue
NEG_COLOR = "#d62728"   # red
NEUTRAL_COLOR = "#7f7f7f"  # gray

# ================== UTILITY FUNCTIONS ==================
def in_core(histone_sheet: str, site: int) -> bool:
    """Check if site is within the core region."""
    start, end = CORE_RANGES[histone_sheet]
    if site < start:
        return False
    if end is not None and site > end:
        return False
    return True

def parse_residue_real_site(x):
    """Parse 'E56' into (orig_aa, site, normalized_str) or (None,None,None)."""
    if pd.isna(x):
        return None, None, None
    s = str(x).strip()
    m = re.match(r"^([A-Za-z])(\d+)$", s)
    if not m:
        return None, None, None
    orig = m.group(1).upper()
    site = int(m.group(2))
    norm = f"{orig}{site}"
    return orig, site, norm

def find_mutation_columns(df):
    """Find columns containing mutation data based on keywords."""
    cols = []
    for c in df.columns:
        cl = str(c).lower()
        for kw in MUT_COL_KEYWORDS:
            if kw.lower() in cl:
                cols.append(c)
                break
    seen = set()
    out = []
    for c in cols:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out

def split_items(s):
    """Split semicolon-separated string into list of (token, count)."""
    if pd.isna(s) or str(s).strip() == "":
        return []
    parts = [p.strip() for p in str(s).split(";") if p.strip()]
    out = []
    for p in parts:
        m = re.match(r"^(.*?)(?:\s*:\s*(\d+))?$", p)
        if not m:
            continue
        token = m.group(1).strip()
        cnt = int(m.group(2)) if m.group(2) else 1
        out.append((token, cnt))
    return out

def parse_token(token: str):
    """Parse a mutation token like 'H2A,Q93077,E56K' into dict or None."""
    fields = [f.strip() for f in token.split(",")]
    if len(fields) < 3:
        return None
    change = fields[2].strip()
    m = re.match(r"^([A-Za-z])(\d+)([A-Za-z])$", change)
    if not m:
        return None
    orig = m.group(1).upper()
    site = int(m.group(2))
    mut = m.group(3).upper()
    return {
        "orig": orig,
        "mut": mut,
        "orig_site_norm": f"{orig}{site}",
    }

def get_dest_color(res):
    """Return color based on charge of destination residue."""
    if res in POS:
        return POS_COLOR
    elif res in NEG:
        return NEG_COLOR
    else:
        return NEUTRAL_COLOR

# ================== MAIN PROCESSING ==================
def main():
    # Data structure: data[site] = {destination_residue: total_count}
    data = defaultdict(lambda: defaultdict(int))
    site_order = []

    df = pd.read_excel(INPUT_FILE, sheet_name=SHEET_NAME, engine="openpyxl")

    if RESIDUE_COL not in df.columns:
        print(f"Missing column '{RESIDUE_COL}' in sheet {SHEET_NAME}.")
        return

    mut_cols = find_mutation_columns(df)
    if not mut_cols:
        print(f"No mutation columns found in sheet {SHEET_NAME}.")
        return

    for idx, row in df.iterrows():
        rr = row.get(RESIDUE_COL)
        rr_orig, rr_site, rr_norm = parse_residue_real_site(rr)
        if rr_site is None:
            continue
        # Core region check
        if not in_core(SHEET_NAME, rr_site):
            continue

        # Collect all mutation tokens and their counts from all mutation columns
        all_tokens = []
        total_count = 0
        for col in mut_cols:
            for token, cnt in split_items(row.get(col)):
                parsed = parse_token(token)
                if not parsed:
                    continue
                # Match by original amino acid only (ignore site number)
                if parsed["orig"] != rr_orig:
                    continue
                all_tokens.append((parsed["mut"], cnt))
                total_count += cnt

        # Filter: keep site only if total_count >= 5 (H4-specific threshold)
        if total_count < 5:
            continue

        # Store data
        site_key = rr_norm
        if site_key not in site_order:
            site_order.append(site_key)
        for dest_res, cnt in all_tokens:
            data[site_key][dest_res] += cnt

    # ================== PLOTTING ==================
    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = "Arial"

    if not site_order:
        print("No sites with total frequency >= 5 found for H4.")
        return

    # Calculate total frequency per site for sorting descending
    site_totals = {s: sum(data[s].values()) for s in site_order}
    sites_sorted = sorted(site_order, key=lambda s: site_totals[s], reverse=True)

    n_items = len(sites_sorted)
    bar_height = min(0.8, 0.8 * n_items / 10)

    fig, ax = plt.subplots(figsize=(10, 7))
    lefts = [0] * len(sites_sorted)

    for site_idx, site in enumerate(sites_sorted):
        dest_items = sorted(data[site].items(), key=lambda x: x[1], reverse=True)
        for dest_res, cnt in dest_items:
            color = get_dest_color(dest_res)
            ax.barh(site_idx, cnt, left=lefts[site_idx],
                    height=bar_height, color=color, edgecolor='white')
            if cnt > 0:
                x_text = lefts[site_idx] + cnt / 2
                ax.text(x_text, site_idx, dest_res,
                        va='center', ha='center',
                        color='white', fontsize=18, fontweight='bold')
            lefts[site_idx] += cnt

    # Y-axis labels
    ax.set_yticks(range(len(sites_sorted)))
    ax.set_yticklabels(sites_sorted)
    # Color tick labels by original residue charge
    for label, site in zip(ax.get_yticklabels(), sites_sorted):
        orig = site[0]
        if orig in POS:
            label.set_color(POS_COLOR)
        elif orig in NEG:
            label.set_color(NEG_COLOR)
        else:
            label.set_color("black")

    if len(sites_sorted) == 1:
        ax.set_ylim(-0.5, 0.5)

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    ax.set_xlabel("Frequency", fontsize=20, fontweight='bold')
    ax.set_title(f"Detailed mutation profile of residues in H4 (frequency ≥ 5)", fontsize=24, pad=15, fontweight='bold')
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

    ax.grid(axis='x', linestyle=':', linewidth=0.6, alpha=0.7)
    ax.grid(axis='y', visible=False)

    # Legend
    legend_handles = [
        Patch(facecolor=POS_COLOR, label='Positive charge (K/R)'),
        Patch(facecolor=NEG_COLOR, label='Negative charge (D/E)'),
        Patch(facecolor=NEUTRAL_COLOR, label='Neutral charge'),
    ]
    ax.legend(handles=legend_handles, loc='upper right', prop={'size': 20})

    max_total = max(site_totals.values()) if sites_sorted else 1
    ax.set_xlim(0, max_total * 1.05 if max_total > 0 else 1)

    plt.tight_layout()

    save_path = os.path.join(OUT_DIR, "residue_detail_mutations_H4_freq5.tif")
    fig.savefig(save_path, dpi=600, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)

    print(f"Plot saved to: {save_path}")

if __name__ == "__main__":
    main()