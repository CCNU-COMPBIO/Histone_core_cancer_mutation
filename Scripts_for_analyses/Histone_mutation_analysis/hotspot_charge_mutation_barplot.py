#!/usr/bin/env python3
# encoding: utf-8
"""
Compute counts of charge-reversal / change / no-change per histone
from Histone_mutation.xlsx, only for core-region sites where the total
mutation count at that site (using row's residue_real_site) >= 8.

Matching rule: mutations with the same original amino acid (regardless of site number)
are counted toward the site defined in the row.
Example: row has E56, mutations E55K and E100R are both counted as E56.

Input: Histone_mutation.xlsx (sheets: H2A, H2B, H3, H4)
Output: charge_mutation_by_histone.tif
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ---- configuration ----
INPUT_FILE = "Histone_mutation.xlsx"          # Updated filename
OUTPUT_PNG = "charge_mutation_by_histone.tif"

SHEETS = ["H2A", "H2B", "H3", "H4"]
HISTONES = SHEETS[:]
ALL_LABEL = "All Histones"

CORE_RANGES = {
    "H2A": (14, 118),
    "H2B": (24, None),
    "H3":  (37, None),
    "H4":  (21, None),
}

RESIDUE_COL = "residue_real_site"
MUT_COL_KEYWORDS = [
    "replication-dependent",
    "replication-independent",
]

POSITIVE = {"K", "R"}
NEGATIVE = {"D", "E"}

sns.set_style("whitegrid")
plt.rcParams["font.family"] = "Arial"


def aa_charge_str(aa: str) -> str:
    """Return charge category of a single amino acid: positive, negative, or neutral."""
    if not aa or len(aa) != 1:
        return "neutral"
    aa = aa.upper()
    if aa in POSITIVE:
        return "positive"
    if aa in NEGATIVE:
        return "negative"
    return "neutral"


def in_core(histone_sheet: str, site: int) -> bool:
    """Check if a residue position falls within the defined core region."""
    start, end = CORE_RANGES[histone_sheet]
    if site < start:
        return False
    if end is not None and site > end:
        return False
    return True


def parse_residue_real_site(x):
    """Parse residue_real_site string into original amino acid, site number, and normalized name."""
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


def find_mutation_columns(df: pd.DataFrame):
    """Find columns containing mutation data based on keywords."""
    cols = []
    for c in df.columns:
        cl = str(c).lower()
        for kw in MUT_COL_KEYWORDS:
            if kw.lower() in cl:
                cols.append(c)
                break
    # Deduplicate while preserving order
    seen = set()
    out = []
    for c in cols:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def split_items(s):
    """Split a semicolon-separated string into tokens with counts."""
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
    """Parse a mutation token into components: histone, Uniprot, change, etc."""
    fields = [f.strip() for f in token.split(",")]
    if len(fields) < 3:
        return None
    hist = fields[0].upper()
    uniprot = fields[1]
    change = fields[2].strip()
    m = re.match(r"^([A-Za-z])(\d+)([A-Za-z])$", change)
    if not m:
        return None
    orig = m.group(1).upper()
    site = int(m.group(2))
    mut = m.group(3).upper()
    return {
        "histone_in_token": hist,
        "uniprot": uniprot,
        "change": change,
        "orig": orig,
        "site": site,
        "mut": mut,
        "orig_site_norm": f"{orig}{site}",
    }


def categorize_charge(orig: str, mut: str) -> str:
    """Classify a mutation as charge reversal, change, or no change."""
    oc = aa_charge_str(orig)
    mc = aa_charge_str(mut)
    if oc == mc:
        return "no_change"
    if (oc == "positive" and mc == "negative") or (oc == "negative" and mc == "positive"):
        return "reversal"
    return "change"


def main():
    xlsx = Path(INPUT_FILE)
    if not xlsx.exists():
        print(f"File not found: {xlsx.resolve()}", file=sys.stderr)
        sys.exit(1)

    categories = ["reversal", "change", "no_change"]
    counts = {h: {c: 0 for c in categories} for h in HISTONES}

    for sheet in SHEETS:
        df = pd.read_excel(xlsx, sheet_name=sheet, engine="openpyxl")

        if RESIDUE_COL not in df.columns:
            print(f"[{sheet}] Missing required column: {RESIDUE_COL}", file=sys.stderr)
            sys.exit(1)

        mut_cols = find_mutation_columns(df)
        if not mut_cols:
            print(f"[{sheet}] Could not find mutation columns.", file=sys.stderr)
            sys.exit(1)

        for idx, row in df.iterrows():
            rr = row.get(RESIDUE_COL)
            rr_orig, rr_site, rr_norm = parse_residue_real_site(rr)
            if rr_site is None:
                continue
            # Core region filter: use the site from the row
            if not in_core(sheet, rr_site):
                continue

            # ---- Collect all mutation entries and their counts for this site ----
            all_tokens = []
            total_site_count = 0
            for col in mut_cols:
                for token, cnt in split_items(row.get(col)):
                    all_tokens.append((token, cnt))
                    total_site_count += cnt

            # ---- Site frequency filter (total count >= 8) ----
            if total_site_count < 8:
                continue

            # ---- Categorize each mutation entry for this retained site ----
            for token, cnt in all_tokens:
                parsed = parse_token(token)
                if not parsed:
                    continue
                # Match only by original amino acid (ignore site number)
                if parsed["orig"] != rr_orig:
                    continue
                cat = categorize_charge(parsed["orig"], parsed["mut"])
                counts[sheet][cat] += cnt

    # ---- Aggregate across all histones ----
    total_all = {c: 0 for c in categories}
    for h in HISTONES:
        for c in categories:
            total_all[c] += counts[h][c]

    # ---- Prepare data for plotting ----
    result_counts = {}
    for h in HISTONES:
        result_counts[h] = counts[h]
    result_counts[ALL_LABEL] = total_all

    print("Counts of mutation charge categories per histone (site total >= 8; matched by original amino acid only):")
    print(pd.DataFrame(result_counts).T.round(0))

    # ---- Plot ----
    plot_df = pd.DataFrame(result_counts).T.reset_index().rename(columns={"index": "histone"})
    desired_order = HISTONES + [ALL_LABEL]
    plot_df['histone'] = pd.Categorical(plot_df['histone'], categories=desired_order, ordered=True)
    plot_df = plot_df.sort_values('histone')

    long = plot_df.melt(id_vars="histone", value_vars=categories,
                        var_name="category", value_name="count")

    color_map = {
        "reversal": "#e53935",   # Bright red
        "change": "#ff9100",     # Bright orange
        "no_change": "#38b000",  # Bright green
    }
    display_names = {
        "reversal": "ΔCharge = ±2",
        "change": "ΔCharge = ±1",
        "no_change": "ΔCharge = 0"
    }

    plt.figure(figsize=(10, 7))
    ax = sns.barplot(data=long, x="histone", y="count", hue="category",
                     palette=color_map, hue_order=categories)

    # Styling
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')
    ax.tick_params(color='black', which='both')

    ax.set_ylabel("Count", fontsize=22, fontweight='bold')
    ax.set_xlabel("", fontsize=10, fontweight='bold')
    ax.set_title("Charge changes at mutation hotspots in histones", fontsize=24, pad=15, fontweight='bold')
    ax.tick_params(axis='x', labelsize=22)
    ax.tick_params(axis='y', labelsize=22)

    max_count = long['count'].max() if not long.empty else 0.0
    ax.set_ylim(0, max_count + max_count * 0.1 + 1)

    handles, labels = ax.get_legend_handles_labels()
    new_labels = [display_names.get(l, l) for l in labels]
    legend = ax.legend(handles=handles, labels=new_labels,
                       loc='upper left', borderaxespad=0.5,
                       fontsize=22)

    # Optional: add count annotations (currently commented out)
    # for p in ax.patches:
    #     height = p.get_height()
    #     if height > 0.5:
    #         ax.annotate(f"{int(height)}", (p.get_x() + p.get_width() / 2., height),
    #                     ha='center', va='bottom', fontsize=20, fontweight='bold',
    #                     color='black', xytext=(0, 3), textcoords='offset points')

    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    print(f"Plot saved to: {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()