#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rmsf_convergence_replica_only.py

RMSF convergence: compare only across the three runs; no split-half of a single trajectory.
Trim 5 residues from each end of every chain (same residue set as S_RMSF); set TRIM to 0 if not needed.

Directory names are auto-detected: <histone>_<mut>_rw (e.g. H2A_R29Q_rw, H2B_E105K_rw).
WT supports WT_rw / H2A_WT_rw / WT_H2A_rw / H3_WT_rw.

Outputs (one row per system-chain):
  conv_summary.csv          excluded version: mean PCC (average over the 3 run pairs) + mean SEM (average of per-residue SEM)
  conv_sem.csv              per-residue SEM details, used to show the "range of SEM"
  conv_summary_noexcl.csv   no-exclusion version (only for systems that have _noexcl)
  conv_replica.csv          pairwise run PCC details
"""

from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

# ========== Configuration ==========
BASE_DIR = Path(__file__).resolve().parent

RUNS = ["run1", "run2", "run3"]
NOEXCL = "_noexcl"

TRIM_N_TERM = 5
TRIM_C_TERM = 5

# Key must match the folder name (without "_rw")
#   Folder lookup: f"{mut}_rw" -> H2A_R29Q_rw / H2B_E105K_rw / H3_E50K_rw ...
#   Also fixes the original bug where "E105K" was defined twice and the H2B version was overwritten.
MUTATION_MAP = {
    "H2A_R29Q":  {"histone": "H2A", "chains": ("C", "D")},
    "H2A_E56K":  {"histone": "H2A", "chains": ("C", "D")},
    "H2A_K74N":  {"histone": "H2A", "chains": ("C", "D")},
    "H2B_E105K": {"histone": "H2B", "chains": ("C", "D")},
    "H2B_F70L":  {"histone": "H2B", "chains": ("C", "D")},
    "H3_E50K":   {"histone": "H3",  "chains": ("A", "B")},
    "H3_E73K":   {"histone": "H3",  "chains": ("A", "B")},
    "H3_E105K":  {"histone": "H3",  "chains": ("A", "B")},
}

CHAIN_TO_HISTONE = {
    "A": "H3", "B": "H4", "C": "H2A", "D": "H2B",
    "E": "H3", "F": "H4", "G": "H2A", "H": "H2B",
}


# ========== Helpers ==========
def read_rmsf(path, trim_n=TRIM_N_TERM, trim_c=TRIM_C_TERM):
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None,
                     names=["Res", "RMSF"])
    df["Res"] = df["Res"].round().astype(int)
    df = df.sort_values("Res").reset_index(drop=True)
    n = len(df)
    if (trim_n > 0 or trim_c > 0) and n > (trim_n + trim_c):
        df = df.iloc[trim_n: n - trim_c].reset_index(drop=True)
    return df


def load_aligned(files):
    """files: {run: Path}; returns (res, {run: ndarray}) aligned to common residues."""
    data = {r: read_rmsf(p) for r, p in files.items()}
    common = None
    for df in data.values():
        s = set(df["Res"].tolist())
        common = s if common is None else (common & s)
    if not common:
        return None, None
    common = sorted(common)
    idx = {r: i for i, r in enumerate(common)}
    out = {}
    for run, df in data.items():
        y = np.full(len(common), np.nan)
        for r, v in zip(df["Res"].tolist(), df["RMSF"].tolist()):
            j = idx.get(r)
            if j is not None:
                y[j] = v
        out[run] = y
    return np.array(common), out


def mean_pairwise_pcc(Y):
    vals = []
    for r1, r2 in combinations(sorted(Y), 2):
        x, y = Y[r1], Y[r2]
        if len(x) > 2 and np.std(x) > 0 and np.std(y) > 0:
            vals.append(float(stats.pearsonr(x, y)[0]))
    return (float(np.mean(vals)) if vals else np.nan), vals


# ========== Directory-name resolution ==========
def resolve_systems():
    systems = {}
    for mut, info in MUTATION_MAP.items():
        cand = f"{mut}_rw"
        if (BASE_DIR / cand).is_dir():
            systems[cand] = {"histone": info["histone"],
                             "chains": info["chains"], "mut": mut}
        else:
            print(f"[WARN] Directory not found for {mut} (expected {cand})")
    for cand in ("WT_rw", "H2A_WT_rw", "WT_H2A_rw", "H3_WT_rw"):
        if (BASE_DIR / cand).is_dir():
            systems[cand] = {"histone": None, "chains": ("A", "B", "C", "D"),
                             "mut": "WT"}
            break
    else:
        print("[WARN] WT directory not found")
    return systems


def runs_with_noexcl(system):
    out = []
    for run in RUNS:
        d = BASE_DIR / system / run
        if any((d / f"rmsf-{h}_chain_{c}{NOEXCL}.dat").exists()
               for h, c in (("H3", "A"), ("H2A", "C"))):
            out.append(run)
    return out


# ========== Computation ==========
def run_variant(systems, suffix, skip_map, wt_runs, label):
    rows, sem_rows, detail = [], [], []
    print(f"\n### VARIANT: {label}  (suffix='{suffix or 'none'}', WT runs={wt_runs})")

    for system, info in systems.items():
        skip = set(skip_map.get(system, []))
        use_runs = [r for r in (wt_runs if info["mut"] == "WT" else RUNS)
                    if r not in skip]
        if len(use_runs) < 2:
            print(f"[SKIP] {system}: fewer than 2 runs")
            continue

        for chain in info["chains"]:
            fh = CHAIN_TO_HISTONE[chain]
            files = {}
            for run in use_runs:
                p = BASE_DIR / system / run / f"rmsf-{fh}_chain_{chain}{suffix}.dat"
                if p.exists():
                    files[run] = p
            if len(files) < 2:
                print(f"[SKIP] {system} chain {chain}: fewer than 2 valid files")
                continue

            res, Y = load_aligned(files)
            if Y is None:
                print(f"[SKIP] {system} chain {chain}: no common residues")
                continue

            runs_sorted = sorted(Y)
            mean_pcc, vals = mean_pairwise_pcc(Y)
            for (r1, r2), v in zip(combinations(runs_sorted, 2), vals):
                detail.append({"system": system, "chain": chain, "histone": fh,
                               "pair": f"{r1}_vs_{r2}", "PCC": v})

            # Per-residue SEM across replicates: sd(ddof=1)/sqrt(n)
            M = np.vstack([Y[r] for r in runs_sorted])      # (n_run, n_res)
            n_run = M.shape[0]
            sd = M.std(axis=0, ddof=1)
            sem = sd / np.sqrt(n_run)
            mean_rmsf = float(M.mean())

            for r, s_ in zip(res, sem):
                sem_rows.append({"system": system, "chain": chain, "histone": fh,
                                 "resid": int(r), "SEM": float(s_)})

            rows.append({
                "system": system, "chain": chain, "histone": fh,
                "n_runs": n_run,
                "mean_PCC": mean_pcc,
                "mean_SEM": float(np.mean(sem)),
                "min_SEM": float(np.min(sem)),
                "max_SEM": float(np.max(sem)),
                "mean_RMSF": mean_rmsf,
                "SEM_over_mean_RMSF": float(np.mean(sem) / mean_rmsf) if mean_rmsf else np.nan,
            })
            print(f"  [OK] {system:16s} chain {chain} ({fh})  runs={runs_sorted}  "
                  f"mean PCC = {mean_pcc:.3f}  mean SEM = {np.mean(sem):.3f} Å")

    cols = ["system", "chain", "histone", "n_runs", "mean_PCC",
            "mean_SEM", "min_SEM", "max_SEM", "mean_RMSF", "SEM_over_mean_RMSF"]
    df_sum = pd.DataFrame(rows, columns=cols)
    df_sem = pd.DataFrame(sem_rows, columns=["system", "chain", "histone",
                                             "resid", "SEM"])
    df_rep = pd.DataFrame(detail, columns=["system", "chain", "histone",
                                           "pair", "PCC"])

    df_sum.to_csv(BASE_DIR / f"conv_summary{suffix}.csv", index=False, float_format="%.4f")
    df_sem.to_csv(BASE_DIR / f"conv_sem{suffix}.csv", index=False, float_format="%.4f")
    df_rep.to_csv(BASE_DIR / f"conv_replica{suffix}.csv", index=False, float_format="%.4f")

    print(f"[OK] conv_summary{suffix}.csv  rows={len(df_sum)}")
    print(f"[OK] conv_sem{suffix}.csv      rows={len(df_sem)}")
    return df_sum


# ========== Main ==========
def main():
    systems = resolve_systems()
    if not systems:
        return

    noexcl = {s: runs_with_noexcl(s) for s, i in systems.items() if i["mut"] != "WT"}
    noexcl = {s: r for s, r in noexcl.items() if r}
    wt_runs_b = sorted({r for rs in noexcl.values() for r in rs})

    print("Resolved systems:")
    for s, i in systems.items():
        print(f"  {s:18s} mut={i['mut']:10s} chains={i['chains']}")
    print("Systems with _noexcl: " + (", ".join(f"{s}{noexcl[s]}" for s in noexcl) or "none"))

    sum_a = run_variant(systems, "", {}, list(RUNS), "excluded segments")
    if noexcl:
        skip_map = {s: [r for r in RUNS if r not in rs] for s, rs in noexcl.items()}
        sum_b = run_variant(systems, NOEXCL, skip_map, wt_runs_b,
                            "no excluded segments")
        if len(sum_a) and len(sum_b):
            m = (sum_a[sum_a["system"].isin(noexcl)]
                 .merge(sum_b[sum_b["system"].isin(noexcl)],
                        on=["system", "chain"], suffixes=("_excl", "_noexcl")))
            m.to_csv(BASE_DIR / "conv_compare_excl_vs_noexcl.csv",
                     index=False, float_format="%.4f")
            print("\n[Comparison] excluded vs no-exclusion")
            print(m[["system", "chain", "mean_PCC_excl", "mean_PCC_noexcl",
                     "mean_SEM_excl", "mean_SEM_noexcl"]]
                  .to_string(index=False, float_format="%.3f"))


if __name__ == "__main__":
    main()