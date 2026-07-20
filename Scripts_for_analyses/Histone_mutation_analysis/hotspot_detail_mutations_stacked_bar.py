#!/usr/bin/env python3
# encoding: utf-8
"""
Generate stacked bar plots for each histone showing mutation destinations per residue site.
Reads Histone_mutation.xlsx (sheets: H2A, H2B, H3, H4).
Only keeps core-region sites with total mutation count >= 8.
Matches mutations by original amino acid only (ignores site number in mutation token).
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
SHEETS = ["H2A", "H2B", "H3", "H4"]
OUT_DIR = "plots_mutations"
os.makedirs(OUT_DIR, exist_ok=True)

# Core region ranges (start, end), end=None means no upper bound
CORE_RANGES = {
    "H2A": (14, 118),
    "H2B": (24, None),
    "H3":  (37, None),
    "H4":  (21, None),
}
RESIDUE_COL = "residue_real_site"
MUT_COL_KEYWORDS = ["replication-dependent", "replication-independent"]

# Charge sets
POS = {"K", "R"}
NEG = {"D", "E"}

# Colors for destination residues
POS_COLOR = "#1f77b4"   # blue for positive
NEG_COLOR = "#d62728"   # red for negative
NEUTRAL_COLOR = "#7f7f7f"  # gray for neutral

# ================== UTILITY FUNCTIONS ==================
def in_core(histone_sheet: str, site: int) -> bool:
    """Check if a residue position is within the core region of the given histone."""
    start, end = CORE_RANGES[histone_sheet]
    if site < start:
        return False
    if end is not None and site > end:
        return False
    return True

def parse_residue_real_site(x):
    """
    Parse a string like 'E56' into (original_aa, site_number, normalized_name).
    Returns (None, None, None) on failure.
    """
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
    """Return column names in df that contain any of the MUT_COL_KEYWORDS (case-insensitive)."""
    cols = []
    for c in df.columns:
        cl = str(c).lower()
        for kw in MUT_COL_KEYWORDS:
            if kw.lower() in cl:
                cols.append(c)
                break
    # Deduplicate preserving order
    seen = set()
    out = []
    for c in cols:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out

def split_items(s):
    """Split a semicolon-separated string into list of (token, count)."""
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
    """
    Parse a mutation token like 'H2A,Q93077,E56K' into dict with orig, mut, orig_site_norm.
    Returns None if parsing fails.
    """
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

def get_dest_charge(res):
    """Return 'positive', 'negative', or 'neutral' based on the residue."""
    if res in POS:
        return "positive"
    elif res in NEG:
        return "negative"
    else:
        return "neutral"

def get_dest_color(res):
    """Return the color for a destination residue based on its charge."""
    if res in POS:
        return POS_COLOR
    elif res in NEG:
        return NEG_COLOR
    else:
        return NEUTRAL_COLOR

# ================== MAIN PROCESSING ==================
def main():
    # Data structure: data[histone][site][destination_residue] = total_count
    data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    # Keep order of sites per histone (for plotting, sorted by total frequency later)
    site_order = defaultdict(list)

    for sheet in SHEETS:
        df = pd.read_excel(INPUT_FILE, sheet_name=sheet, engine="openpyxl")

        if RESIDUE_COL not in df.columns:
            print(f"[{sheet}] Missing column '{RESIDUE_COL}'. Skipping sheet.")
            continue

        mut_cols = find_mutation_columns(df)
        if not mut_cols:
            print(f"[{sheet}] No mutation columns found. Skipping sheet.")
            continue

        for idx, row in df.iterrows():
            rr = row.get(RESIDUE_COL)
            rr_orig, rr_site, rr_norm = parse_residue_real_site(rr)
            if rr_site is None:
                continue
            # Core region check
            if not in_core(sheet, rr_site):
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

            # Filter: keep site only if total_count >= 8
            if total_count < 8:
                continue

            # Store data
            site_key = rr_norm  # use the row's site label
            if site_key not in site_order[sheet]:
                site_order[sheet].append(site_key)
            for dest_res, cnt in all_tokens:
                data[sheet][site_key][dest_res] += cnt

    # ================== PLOTTING ==================
    sns.set_style("whitegrid")
    plt.rcParams["font.family"] = "Arial"

    def plot_histone(hist, save_path=None, show=False):
        sites = site_order[hist]
        if not sites:
            print(f"No sites to plot for {hist}")
            return

        # Calculate total frequency per site for sorting descending
        site_totals = {s: sum(data[hist][s].values()) for s in sites}
        sites_sorted = sorted(sites, key=lambda s: site_totals[s], reverse=True)

        n_items = len(sites_sorted)
        # Adjust bar height based on number of sites
        bar_height = min(0.8, 0.8 * n_items / 10)

        fig, ax = plt.subplots(figsize=(10, 7))
        lefts = [0] * len(sites_sorted)

        for site_idx, site in enumerate(sites_sorted):
            dest_items = sorted(data[hist][site].items(), key=lambda x: x[1], reverse=True)
            for dest_res, cnt in dest_items:
                color = get_dest_color(dest_res)
                ax.barh(site_idx, cnt, left=lefts[site_idx],
                        height=bar_height, color=color, edgecolor='white')
                if cnt > 0:
                    # Label with destination residue letter
                    x_text = lefts[site_idx] + cnt / 2
                    ax.text(x_text, site_idx, dest_res,
                            va='center', ha='center',
                            color='white', fontsize=18, fontweight='bold')
                lefts[site_idx] += cnt

        # Y-axis labels: site names
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

        # Adjust y-axis for single site
        if len(sites_sorted) == 1:
            ax.set_ylim(-0.5, 0.5)

        for spine in ax.spines.values():
            spine.set_linewidth(2)
            spine.set_color('black')
        ax.tick_params(color='black', which='both')

        ax.set_xlabel("Frequency", fontsize=20, fontweight='bold')
        ax.set_title(f"Detailed mutation profile of residues in {hist}", fontsize=24, pad=15, fontweight='bold')
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)

        # Grid: only vertical
        ax.grid(axis='x', linestyle=':', linewidth=0.6, alpha=0.7)
        ax.grid(axis='y', visible=False)

        # Legend
        legend_handles = [
            Patch(facecolor=POS_COLOR, label='Positive charge (K/R)'),
            Patch(facecolor=NEG_COLOR, label='Negative charge (D/E)'),
            Patch(facecolor=NEUTRAL_COLOR, label='Neutral charge'),
        ]
        ax.legend(handles=legend_handles, loc='upper right', prop={'size': 20})

        # Set x-axis limit with a small margin
        max_total = max(site_totals.values()) if sites_sorted else 1
        ax.set_xlim(0, max_total * 1.05 if max_total > 15 else 15)

        plt.tight_layout()
        if save_path:
            fig.savefig(save_path, dpi=600, bbox_inches='tight',
                        format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
        if show:
            plt.show()
        plt.close(fig)
        return save_path

    # Generate plots for each histone sheet
    saved_files = []
    for hist in sorted(data.keys()):
        fname = f"residue_detail_mutations_{hist}.tif".replace("/", "_")
        path = os.path.join(OUT_DIR, fname)
        plot_histone(hist, save_path=path, show=False)
        saved_files.append(path)

    print("Saved files:")
    for p in saved_files:
        print(p)

if __name__ == "__main__":
    main()