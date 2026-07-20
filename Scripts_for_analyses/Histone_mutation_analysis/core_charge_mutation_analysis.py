#!/usr/bin/env python3
# encoding: utf-8

"""
Stat charge-change categories (ΔCharge = 0/1/2) for histone core-region mutations
from Histone_mutation.xlsx (sheets: H2A, H2B, H3, H4).

Rules:
- Only keep rows whose residue_real_site position is in the histone core region.
- Parse mutation tokens like: "H2A,Q93077,E56K: 1; H2A,Q8IUE6,E56Q: 2; ..."
- Keep a token if its original residue (e.g., E in E56K) matches the row's original residue,
  regardless of the site number.
- Histidine (H) is treated as neutral (0).
- Charge per residue: K/R=+1, D/E=-1, others=0.
- Category uses absolute delta: 0,1,2.

Outputs:
- core_charge_mutation_long.tsv
- core_charge_mutation_summary.tsv
- core_charge_mutation_by_histone.tif
"""

import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ---- configuration ----
INPUT_FILE = "Histone_mutation.xlsx"          # Updated filename
OUTPUT_LONG = "core_charge_mutation_long.tsv"
OUTPUT_SUMMARY = "core_charge_mutation_summary.tsv"
OUTPUT_PNG = "core_charge_mutation_by_histone.tif"

SHEETS = ["H2A", "H2B", "H3", "H4"]

# Core region ranges: (start, end). end=None means no upper bound
CORE_RANGES = {
    "H2A": (14, 118),
    "H2B": (24, None),
    "H3":  (37, None),
    "H4":  (21, None),
}

RESIDUE_COL = "residue_real_site"

# Mutation column keywords (case-insensitive fuzzy search)
MUT_COL_KEYWORDS = [
    "replication-dependent",
    "replication-independent",
]

# Plot style
sns.set_style("whitegrid")
plt.rcParams["font.family"] = "Arial"

POSITIVE = {"K", "R"}   # +1
NEGATIVE = {"D", "E"}   # -1
# Histidine is neutral per requirement, so NOT in POSITIVE


def aa_charge(aa: str) -> int:
    """Return charge: +1 / -1 / 0 for one-letter amino acid (H is 0 here)."""
    if not aa or len(aa) != 1:
        return 0
    aa = aa.upper()
    if aa in POSITIVE:
        return 1
    if aa in NEGATIVE:
        return -1
    return 0


def in_core(histone_sheet: str, site: int) -> bool:
    """Check if site falls within the defined core region."""
    start, end = CORE_RANGES[histone_sheet]
    if site < start:
        return False
    if end is not None and site > end:
        return False
    return True


def parse_residue_real_site(x):
    """
    Parse 'E56' -> ('E', 56, 'E56')
    Return (orig_aa, site, normalized_str) or (None, None, None)
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


def find_mutation_columns(df: pd.DataFrame):
    """
    Find columns whose names contain any keyword in MUT_COL_KEYWORDS (case-insensitive).
    Return list of column names in df.
    """
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
    """
    Split '...; ...; ...' into list of (token, count)
    token example: 'H2A,Q93077,E56K'
    """
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
    token like: 'H2A,Q93077,E56K'
    Return dict or None.
    """
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


def main():
    xlsx = Path(INPUT_FILE)
    if not xlsx.exists():
        print(f"File not found: {xlsx.resolve()}", file=sys.stderr)
        sys.exit(1)

    long_rows = []

    for sheet in SHEETS:
        df = pd.read_excel(xlsx, sheet_name=sheet, engine="openpyxl")

        if RESIDUE_COL not in df.columns:
            print(f"[{sheet}] Missing required column: {RESIDUE_COL}", file=sys.stderr)
            print(f"[{sheet}] Columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)

        mut_cols = find_mutation_columns(df)
        if not mut_cols:
            print(f"[{sheet}] Could not find mutation columns by keywords {MUT_COL_KEYWORDS}", file=sys.stderr)
            print(f"[{sheet}] Columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)

        for idx, row in df.iterrows():
            rr = row.get(RESIDUE_COL)
            rr_orig, rr_site, rr_norm = parse_residue_real_site(rr)
            if rr_site is None:
                continue
            if not in_core(sheet, rr_site):
                continue

            for col in mut_cols:
                for token, cnt in split_items(row.get(col)):
                    parsed = parse_token(token)
                    if not parsed:
                        continue

                    # Match only by original amino acid (ignore site number)
                    if parsed["orig"] != rr_orig:
                        continue

                    # Compute |Δcharge|
                    d = abs(aa_charge(parsed["mut"]) - aa_charge(parsed["orig"]))
                    if d not in (0, 1, 2):
                        continue

                    long_rows.append({
                        "sheet": sheet,
                        "row_index": int(idx),
                        "source_col": str(col),
                        "residue_real_site": rr_norm,

                        "histone_in_token": parsed["histone_in_token"],
                        "uniprot": parsed["uniprot"],
                        "change": parsed["change"],
                        "orig": parsed["orig"],
                        "site": parsed["site"],
                        "mut": parsed["mut"],
                        "count": int(cnt),

                        "delta_charge_abs": int(d),
                    })

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(OUTPUT_LONG, sep="\t", index=False)

    # ---- Summary: per histone (sheet) and all histones ----
    categories = [0, 1, 2]  # |ΔCharge|
    histones = SHEETS[:]
    all_label = "All Histones"

    if long_df.empty:
        summary = pd.DataFrame({
            "histone": histones + [all_label],
            "delta0_count": 0,
            "delta1_count": 0,
            "delta2_count": 0,
            "total": 0,
            "delta0_pct": 0.0,
            "delta1_pct": 0.0,
            "delta2_pct": 0.0,
        })
        summary.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)
        print("No valid mutations found after filtering; outputs written.")
        return

    # Counts per sheet
    count_tbl = (
        long_df
        .groupby(["sheet", "delta_charge_abs"])["count"]
        .sum()
        .unstack(fill_value=0)
    )
    for c in categories:
        if c not in count_tbl.columns:
            count_tbl[c] = 0
    count_tbl = count_tbl[categories]

    # Overall counts across all sheets
    all_counts = long_df.groupby("delta_charge_abs")["count"].sum()
    all_row = pd.Series({c: int(all_counts.get(c, 0)) for c in categories}, name=all_label)

    count_tbl2 = pd.concat([count_tbl, all_row.to_frame().T], axis=0)

    # Build summary with percentages
    rows = []
    for h in histones + [all_label]:
        c0 = int(count_tbl2.loc[h, 0]) if h in count_tbl2.index else 0
        c1 = int(count_tbl2.loc[h, 1]) if h in count_tbl2.index else 0
        c2 = int(count_tbl2.loc[h, 2]) if h in count_tbl2.index else 0
        total = c0 + c1 + c2
        if total == 0:
            p0 = p1 = p2 = 0.0
        else:
            p0 = c0 / total * 100.0
            p1 = c1 / total * 100.0
            p2 = c2 / total * 100.0
        rows.append({
            "histone": h,
            "delta0_count": c0,
            "delta1_count": c1,
            "delta2_count": c2,
            "total": total,
            "delta0_pct": p0,
            "delta1_pct": p1,
            "delta2_pct": p2,
        })

    summary = pd.DataFrame(rows)
    summary.to_csv(OUTPUT_SUMMARY, sep="\t", index=False)

    print("Summary (counts & %):")
    print(summary.round(2).to_string(index=False))

    # ============= Plotting =============
    plot_df = summary.copy()
    plot_df["histone"] = pd.Categorical(plot_df["histone"], categories=histones + [all_label], ordered=True)
    plot_df = plot_df.sort_values("histone")

    long_plot = plot_df.melt(
        id_vars="histone",
        value_vars=["delta0_count", "delta1_count", "delta2_count"],
        var_name="category",
        value_name="count"
    )
    cat_map = {
        "delta0_count": "ΔCharge = 0",
        "delta1_count": "ΔCharge = ±1",
        "delta2_count": "ΔCharge = ±2",
    }
    long_plot["category"] = long_plot["category"].map(cat_map)

    color_map = {
        "ΔCharge = ±2": "#e53935",   # Bright red
        "ΔCharge = ±1": "#ff9100",   # Bright orange
        "ΔCharge = 0": "#38b000",    # Bright green
    }
    hue_order = ["ΔCharge = ±2", "ΔCharge = ±1", "ΔCharge = 0"]

    plt.figure(figsize=(10, 7))
    ax = sns.barplot(
        data=long_plot, x="histone", y="count",
        hue="category", hue_order=hue_order,
        palette=color_map
    )

    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')
    ax.tick_params(color='black', which='both')

    ax.set_ylabel("Count", fontsize=22, fontweight="bold")
    ax.set_xlabel("")
    ax.set_title("Charge change of core-region mutations in histones", fontsize=24, pad=15, fontweight="bold")
    ax.tick_params(axis="x", labelsize=22)
    ax.tick_params(axis="y", labelsize=22)

    max_count = long_plot["count"].max() if not long_plot.empty else 0.0
    ax.set_ylim(0, max_count + max_count * 0.1 + 1)

    # Optionally add value annotations (commented out)
    # for p in ax.patches:
    #     height = p.get_height()
    #     if height > 0.5:
    #         ax.annotate(f"{int(height)}",
    #                     (p.get_x() + p.get_width() / 2., height),
    #                     ha='center', va='bottom',
    #                     fontsize=20, fontweight='bold',
    #                     color='black', xytext=(0, 3),
    #                     textcoords='offset points')

    legend = ax.legend(loc="upper left", fontsize=22)
    try:
        legend.get_title().set_fontsize(22)
    except Exception:
        pass

    plt.tight_layout()
    plt.savefig(
        OUTPUT_PNG, dpi=600, bbox_inches="tight",
        format="tiff", pil_kwargs={"compression": "tiff_lzw"}
    )
    print(f"Files written:\n- {OUTPUT_LONG}\n- {OUTPUT_SUMMARY}\n- {OUTPUT_PNG}")
    plt.show()


if __name__ == "__main__":
    main()