#!/usr/bin/env python3
# encoding: utf-8
"""
Compute counts of charge-reversal / change / no-change for H4 histone
from Histone_mutation.xlsx, sheet H4, only for core-region sites
where total mutation count at that site >= 5.

Matching rule: mutations with same original amino acid (regardless of site number)
are counted toward the site defined in the row.

Input: Histone_mutation.xlsx (sheet H4)
Output: h4_charge_mutation_pie.tif (pie chart)
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ---- configuration ----
INPUT_FILE = "Histone_mutation.xlsx"
OUTPUT_PNG = "h4_charge_mutation_pie.tif"

SHEETS = ["H4"]
HISTONES = SHEETS[:]

CORE_RANGES = {
    "H4": (21, None),
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
    """Return charge category: positive, negative, or neutral."""
    if not aa or len(aa) != 1:
        return "neutral"
    aa = aa.upper()
    if aa in POSITIVE:
        return "positive"
    if aa in NEGATIVE:
        return "negative"
    return "neutral"


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


def find_mutation_columns(df: pd.DataFrame):
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
    """Categorize mutation as reversal, change, or no_change based on charge."""
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
            if not in_core(sheet, rr_site):
                continue

            all_tokens = []
            total_site_count = 0
            for col in mut_cols:
                for token, cnt in split_items(row.get(col)):
                    all_tokens.append((token, cnt))
                    total_site_count += cnt

            if total_site_count < 5:
                continue

            for token, cnt in all_tokens:
                parsed = parse_token(token)
                if not parsed:
                    continue
                # Match by original amino acid only (ignore site number)
                if parsed["orig"] != rr_orig:
                    continue
                cat = categorize_charge(parsed["orig"], parsed["mut"])
                counts[sheet][cat] += cnt

    result_counts = counts["H4"]

    print("Counts of mutation charge categories for H4 (site total >= 5; matched by original amino acid only):")
    print(pd.DataFrame([result_counts], index=["H4"]).round(0))

    # ==================== Plot: Pie chart ====================
    labels_map = {
        "reversal": "ΔCharge = ±2",
        "change": "ΔCharge = ±1",
        "no_change": "ΔCharge = 0"
    }
    colors = {
        "reversal": "#ff9999",  # light red
        "change": "#ffcc99",    # light orange
        "no_change": "#99dd99"  # light green
    }

    ordered_cats = ["reversal", "change", "no_change"]
    values = [result_counts[cat] for cat in ordered_cats]
    labels = [labels_map[cat] for cat in ordered_cats]
    pie_colors = [colors[cat] for cat in ordered_cats]

    if sum(values) == 0:
        print("No data to plot for H4.")
        return

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        values,
        labels=labels,
        autopct=lambda p: f'{p:.1f}%\n({int(p * sum(values) / 100)})',
        colors=pie_colors,
        startangle=90,
        textprops={'fontsize': 20, 'fontweight': 'bold'},
        wedgeprops={'edgecolor': 'black', 'linewidth': 2}
    )
    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontsize(20)
        autotext.set_fontweight('bold')

    ax.set_title("Charge changes in H4 (frequency ≥ 5)", fontsize=24, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight',
                format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    print(f"Plot saved to: {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()