#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RMSF_wt_vs_mut - windows.py
Per-residue RMSF comparison between WT and mutants + difference statistics.

Two variants, generated in one run:
  Variant A "excluded": all systems x run1-run3, reads rmsf-<H>_chain_<C>.dat, output name without suffix
  Variant B "no-exclusion": only processes runs that actually have _noexcl files (K74N run1, F70L run3),
                 reads rmsf-<H>_chain_<C>_noexcl.dat, output name with _noexcl
                 -- systems with no excluded segments and not recomputed are skipped in Variant B,
                    no duplicate computation, no duplicate plots.
WT has no excluded segments, so both variants use the regular files.

Folder naming is unified as: <histone>_<mutation>_rw, e.g. H2B_E105K_rw.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_theme(style="whitegrid", context="paper")
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# ========== Paths and run parameters ==========
BASE_DIR = Path(__file__).resolve().parent
RUNS = ["run1", "run2", "run3"]
SEQUENCE_CSV = "sequence.csv"
WT_SYSTEM = "WT_rw"
NOEXCL = "_noexcl"

# ★ Variant definitions: mode="regular" reads regular files; mode="noexcl" reads only _noexcl files ★
VARIANTS = [
    {"mode": "regular", "out_suffix": "",     "label": "excluded segments"},
    {"mode": "noexcl",  "out_suffix": NOEXCL, "label": "no excluded segments"},
]

CHAIN_START = {
    "A": 38, "E": 37,
    "B": 25, "F": 21,
    "C": 14, "G": 15,
    "D": 30, "H": 33,
}

# ★ Folder names are already unified as H2A_R29Q_rw format ★
#   key is the folder name without "_rw"
MUTATION_MAP = {
    "H2A_R29Q":  {"histone": "H2A", "mutation": "R29Q",  "chain": "C", "x_resnum": 29},
    "H2A_E56K":  {"histone": "H2A", "mutation": "E56K",  "chain": "C", "x_resnum": 56},
    "H2B_E105K": {"histone": "H2B", "mutation": "E105K", "chain": "D", "x_resnum": 105},
    "H2A_K74N":  {"histone": "H2A", "mutation": "K74N",  "chain": "C", "x_resnum": 74},
    "H2B_F70L":  {"histone": "H2B", "mutation": "F70L",  "chain": "D", "x_resnum": 70},
    "H3_E50K":   {"histone": "H3",  "mutation": "E50K",  "chain": "A", "x_resnum": 50},
    "H3_E73K":   {"histone": "H3",  "mutation": "E73K",  "chain": "A", "x_resnum": 73},
    "H3_E105K":  {"histone": "H3",  "mutation": "E105K", "chain": "A", "x_resnum": 105},
}

HISTONE_CHAIN_PAIRS = {
    "H3":  ("A", "B"),
    "H2A": ("C", "D"),
    "H2B": ("C", "D"),
}

CHAIN_TO_HISTONE = {
    "A": "H3", "B": "H4",
    "C": "H2A", "D": "H2B",
    "E": "H3", "F": "H4",
    "G": "H2A", "H": "H2B",
}

PLOT_SEM_BAND = True
DPI = 600
DIFF_THRESHOLD = 0.5
LEGEND_SHOW_SUFFIX = False


# ========== Utility functions ==========
def read_sequence_csv(path: Path) -> pd.DataFrame:
    encodings = ["utf-8", "utf-8-sig", "gbk", "cp936", "latin1"]
    last = None
    for enc in encodings:
        try:
            df = pd.read_csv(path, encoding=enc)
            break
        except UnicodeDecodeError as e:
            last = e
            df = None
    if df is None:
        raise UnicodeDecodeError("unknown", b"", 0, 1,
                                 f"Failed to decode {path}; last={last}")

    rename = {}
    if "#Res" in df.columns: rename["#Res"] = "Res"
    if "chain" in df.columns: rename["chain"] = "Chain"
    if "residue" in df.columns: rename["residue"] = "ResidueName"
    if "#Orig" in df.columns: rename["#Orig"] = "Orig"
    df = df.rename(columns=rename)

    if "Res" not in df.columns or "Chain" not in df.columns:
        raise ValueError(f"{path} must contain #Res/Res and chain/Chain. "
                         f"Got: {list(df.columns)}")

    df["Res"] = pd.to_numeric(df["Res"], errors="coerce")
    df = df.dropna(subset=["Res"]).copy()
    df["Res"] = df["Res"].astype(int)
    df["Chain"] = df["Chain"].astype(str).str.strip()
    return df


def build_chain_res_to_pos(seq_df: pd.DataFrame, chain_id: str) -> pd.DataFrame:
    sub = seq_df[seq_df["Chain"] == chain_id].copy()
    if sub.empty:
        chains = sorted(seq_df["Chain"].unique().tolist())
        raise ValueError(f"Chain '{chain_id}' not found in sequence.csv. "
                         f"Available: {chains}")
    sub = sub.sort_values("Res").reset_index(drop=True)
    sub["Pos"] = np.arange(1, len(sub) + 1)
    return sub[["Res", "Pos"]]


def read_cpptraj_rmsf(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", delim_whitespace=True, header=None,
                     names=["Res", "RMSF"], dtype={"Res": float, "RMSF": float})
    if df.empty:
        raise ValueError(f"{path} parsed as empty dataframe.")
    df["Res"] = df["Res"].round().astype(int)
    return df


def runs_with_noexcl(system) -> list:
    """Return which runs actually have _noexcl files for this system (probe using H3 chain A)."""
    found = []
    for run in RUNS:
        probe = BASE_DIR / system / run / f"rmsf-H3_chain_A{NOEXCL}.dat"
        if probe.exists():
            found.append(run)
    return found


def rmsf_path_for(system, run, file_histone, chain, mode):
    d = BASE_DIR / system / run
    suffix = NOEXCL if mode == "noexcl" else ""
    return d / f"rmsf-{file_histone}_chain_{chain}{suffix}.dat"


def load_system_chain_runs(system, file_histone, chain, res_to_pos,
                           start_index, mode, runs):
    per_run = []
    for run in runs:
        f = rmsf_path_for(system, run, file_histone, chain, mode)
        if not f.exists():
            raise FileNotFoundError(f"Missing file: {f}")
        df = read_cpptraj_rmsf(f)
        merged = df.merge(res_to_pos, on="Res", how="inner")
        if merged.empty:
            raise ValueError(f"No residues matched between {f} and "
                             f"sequence.csv for chain={chain}.")
        merged["X"] = start_index + merged["Pos"] - 1
        merged = merged.sort_values("X")
        per_run.append((merged["X"].to_numpy(dtype=int),
                        merged["RMSF"].to_numpy(dtype=float)))
    return per_run


def align_and_stack(series_list):
    sets = [set(x.tolist()) for x, _ in series_list]
    common = sorted(set.intersection(*sets))
    if not common:
        raise ValueError("No common X positions across runs.")
    idx = {v: i for i, v in enumerate(common)}
    Y = np.full((len(series_list), len(common)), np.nan, dtype=float)
    for i, (x, y) in enumerate(series_list):
        for xv, yv in zip(x, y):
            j = idx.get(int(xv))
            if j is not None:
                Y[i, j] = float(yv)
    valid = ~np.isnan(Y).any(axis=0)
    return np.array(common)[valid], Y[:, valid]


def summarize_runs(per_run):
    common_x, Y = align_and_stack(per_run)
    n_runs = Y.shape[0]
    mean = Y.mean(axis=0)
    std = Y.std(axis=0, ddof=1) if n_runs > 1 else np.zeros_like(mean)
    return common_x, mean, std / np.sqrt(n_runs)


def compute_system_mean_for_chain(system, chain_id, res_to_pos, mode, runs):
    file_histone = CHAIN_TO_HISTONE[chain_id]
    per_run = load_system_chain_runs(
        system, file_histone, chain_id, res_to_pos, CHAIN_START[chain_id], mode, runs)
    return summarize_runs(per_run)


def annotate_single_mutation(ax, mut_label, mut_info, chain_id, res_to_pos, y_map):
    if mut_info["chain"] != chain_id:
        return
    x_mut = int(mut_info["x_resnum"])
    pos_mut = int(x_mut - CHAIN_START[chain_id] + 1)
    if pos_mut < 1 or not (res_to_pos["Pos"] == pos_mut).any():
        print(f"[WARN] Mutation {mut_label} (X={x_mut}) not in chain {chain_id}")
        return
    if x_mut not in y_map:
        print(f"[WARN] Mutation {mut_label} X={x_mut} not found on mean curve "
              f"for chain {chain_id}")
        return
    ax.plot([x_mut], [y_map[x_mut]],
            marker="o", markersize=10,
            markerfacecolor="crimson", markeredgecolor="white",
            markeredgewidth=1.2, linestyle="None", zorder=10)


def diff_table_for_chain(chain_id, res_to_pos, sys_a, sys_b,
                         mode_a, mode_b, runs_a, runs_b):
    x_a, mean_a, _ = compute_system_mean_for_chain(sys_a, chain_id, res_to_pos,
                                                   mode_a, runs_a)
    x_b, mean_b, _ = compute_system_mean_for_chain(sys_b, chain_id, res_to_pos,
                                                   mode_b, runs_b)
    common = np.intersect1d(x_a, x_b)
    if common.size == 0:
        raise ValueError(f"No common X between {sys_a} and {sys_b} "
                         f"for chain {chain_id}")

    idx_a = {int(x): i for i, x in enumerate(x_a)}
    idx_b = {int(x): i for i, x in enumerate(x_b)}
    rows = []
    for x in common:
        wt = float(mean_a[idx_a[int(x)]])
        mut = float(mean_b[idx_b[int(x)]])
        rows.append({
            "chain": chain_id, "resid": int(x),
            f"{sys_a}_mean": wt, f"{sys_b}_mean": mut,
            "diff": mut - wt, "abs_diff": abs(mut - wt),
        })
    return pd.DataFrame(rows).sort_values("resid").reset_index(drop=True)


def print_diff_hits(df, threshold, label):
    hits = df[df["RMSF_delta_abs"] >= float(threshold)].copy()
    if hits.empty:
        print(f"[INFO] [{label}] No sites with |RMSF_delta| >= {threshold}")
        return
    hits = hits.sort_values("RMSF_delta_abs", ascending=False)
    print(f"\n[DIFF] [{label}] Sites with |RMSF_delta| >= {threshold}")
    for _, r in hits.iterrows():
        print(f"  Chain {r['chain']}  resid={int(r['resid'])}  "
              f"RMSF_delta={r['RMSF_delta']:.3f}  (abs={r['RMSF_delta_abs']:.3f})")


# ========== Main workflow ==========
def main():
    seq_path = BASE_DIR / SEQUENCE_CSV
    if not seq_path.exists():
        raise FileNotFoundError(f"Cannot find: {seq_path}")
    seq_df = read_sequence_csv(seq_path)
    pal = sns.color_palette("deep", 2)
    color_wt, color_mut = pal[0], pal[1]

    print("=" * 68)
    print("[Variant B scan] Which runs have _noexcl files (only these runs will use Variant B)")
    noexcl_map = {}
    for mut_id in MUTATION_MAP:
        sys_name = f"{mut_id}_rw"
        r = runs_with_noexcl(sys_name)
        noexcl_map[sys_name] = r
        print(f"  {sys_name:15s} -> {r if r else '(none, Variant B skipped)'}")
    print("=" * 68)

    for mut_id, mut_info in MUTATION_MAP.items():
        histone = mut_info["histone"]
        mutation = mut_info["mutation"]

        mut_system = f"{mut_id}_rw"              # e.g. H2B_E105K_rw
        display_label = f"{histone}{mutation}"   # e.g. H2BE105K
        file_stem = mut_id                       # e.g. H2B_E105K (use key directly)

        chain1, chain2 = HISTONE_CHAIN_PAIRS[histone]

        chain_to_res_to_pos = {}
        for cid in (chain1, chain2):
            try:
                chain_to_res_to_pos[cid] = build_chain_res_to_pos(seq_df, cid)
            except ValueError as e:
                print(f"[ERROR] {e}")

        for variant in VARIANTS:
            mode = variant["mode"]
            out_sfx = variant["out_suffix"]
            vlabel = variant["label"]

            # ---------- Variant B only processes runs with _noexcl ----------
            if mode == "noexcl":
                mut_runs = noexcl_map.get(mut_system, [])
                if not mut_runs:
                    print(f"\n[SKIP-B] {display_label}: no _noexcl files, "
                          f"Variant B will not recompute (same as Variant A)")
                    continue
                wt_runs = list(mut_runs)
                wt_mode = "regular"
            else:
                mut_runs = list(RUNS)
                wt_runs = list(RUNS)
                wt_mode = "regular"

            print(f"\n{'='*68}")
            print(f"{display_label} ({mut_system}) | chains {chain1},{chain2} | {vlabel}")
            print(f"  mutant runs = {mut_runs} | WT runs = {wt_runs} | WT mode = {wt_mode}")
            print(f"{'='*68}")

            for chain_id in (chain1, chain2):
                if chain_id not in chain_to_res_to_pos:
                    print(f"[SKIP] Chain {chain_id} not available")
                    continue

                res_to_pos = chain_to_res_to_pos[chain_id]
                file_histone = CHAIN_TO_HISTONE[chain_id]

                # ---------- Difference statistics ----------
                try:
                    df_diff = diff_table_for_chain(
                        chain_id, res_to_pos,
                        sys_a=WT_SYSTEM, sys_b=mut_system,
                        mode_a=wt_mode, mode_b=mode,
                        runs_a=wt_runs, runs_b=mut_runs)

                    df_diff.rename(columns={
                        f"{WT_SYSTEM}_mean": "RMSF_WT",
                        f"{mut_system}_mean": f"RMSF_{display_label}",
                        "diff": "RMSF_delta",
                        "abs_diff": "RMSF_delta_abs",
                    }, inplace=True)

                    print_diff_hits(df_diff, DIFF_THRESHOLD,
                                    f"{display_label} [{vlabel}]")

                    out_csv = (BASE_DIR /
                               f"RMSF_delta_{file_stem}_chain{chain_id}{out_sfx}.csv")
                    df_diff.to_csv(out_csv, index=False, float_format='%.4f')
                    print(f"  [OK] {out_csv.name}")
                except Exception as e:
                    print(f"[ERROR] Diff failed for {display_label} chain {chain_id}: {e}")
                    continue

                # ---------- Plotting ----------
                fig, ax = plt.subplots(figsize=(8, 4))
                mean_by_system = {}

                ls_extra = NOEXCL if (mode == "noexcl" and LEGEND_SHOW_SUFFIX) else ""
                for system, color, ls, legend_label, m, rr in [
                    (WT_SYSTEM, color_wt, "-", "WT", wt_mode, wt_runs),
                    (mut_system, color_mut, "--", display_label + ls_extra, mode, mut_runs),
                ]:
                    try:
                        x, mean, sem = compute_system_mean_for_chain(
                            system, chain_id, res_to_pos, m, rr)
                    except Exception as e:
                        print(f"[ERROR] Cannot load {system} chain {chain_id} "
                              f"(mode={m}, runs={rr}): {e}")
                        continue

                    mean_by_system[system] = (x, mean)
                    ax.plot(x, mean, color=color, linestyle=ls,
                            linewidth=2.0, label=legend_label)
                    if PLOT_SEM_BAND and len(rr) > 1:
                        ax.fill_between(x, mean - sem, mean + sem,
                                        color=color, alpha=0.20, linewidth=0)

                if mut_system in mean_by_system:
                    x_curve, y_curve = mean_by_system[mut_system]
                    y_map = {int(x): float(y) for x, y in zip(x_curve, y_curve)}
                    annotate_single_mutation(ax, display_label, mut_info,
                                             chain_id, res_to_pos, y_map)

                ax = plt.gca()
                for spine in ax.spines.values():
                    spine.set_linewidth(2)
                    spine.set_color('black')

                ax.set_xlabel("Residue index", fontsize=20, weight='bold')
                ax.set_ylabel("RMSF (Å)", fontsize=20, weight='bold')
                ax.set_title(f"WT vs {display_label} — {CHAIN_TO_HISTONE[chain_id]} "
                             f"(chain {chain_id})", fontsize=22, pad=10, weight='bold')
                ax.grid(True, alpha=0.25)
                ax.legend(frameon=False, ncol=2, fontsize=20, loc="upper center")
                for lab in ax.get_xticklabels() + ax.get_yticklabels():
                    lab.set_fontsize(18)

                fig.tight_layout()
                out = BASE_DIR / (f"RMSF_WT_vs_{file_stem}_chain_{chain_id}"
                                  f"{out_sfx}_seaborn.png")
                fig.savefig(out, bbox_inches='tight', dpi=DPI)
                print(f"  [OK] {out.name}")
                plt.show()


if __name__ == "__main__":
    main()