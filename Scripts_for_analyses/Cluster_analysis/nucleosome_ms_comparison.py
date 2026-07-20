#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nucleosome_ms_comparison.py
Calculate and compare Ms values (mean(1/r_ij)) between two residue sets,
with significance testing using permutation.
Uses residue centroids (geometric center of all atoms) rather than C-alpha atoms.
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
import random

# Parameters
EXCEL_PATH = "AlanineScan_DDG_Nucleosome_Averages.xlsx"
PDB_LOCAL = "3afa.pdb"
PDB_ID = "3afa"
CHAIN_MAPPING = {'C': 'H2A', 'D': 'H2B', 'A': 'H3', 'B': 'H4'}

# Statistical parameters
BOOTSTRAP_ITERATIONS = 10000
SEED = 42

# Residue name mapping (one-letter to three-letter)
ONE2THREE = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}


def get_residue_centroid(residue):
    """Compute the centroid (geometric center of all atoms) for a residue."""
    atom_coords = [atom.get_coord() for atom in residue]
    if atom_coords:
        return np.mean(atom_coords, axis=0)
    return None


def read_functional_excel(path):
    """Read functional data from the four histone sheets."""
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        print("Biopython is required. Please install: pip install biopython")
        exit(1)

    sheets = ['H2A', 'H2B', 'H3', 'H4']
    frames = []
    for s in sheets:
        df = pd.read_excel(path, sheet_name=s)
        df.columns = df.columns.str.strip()
        df['histone_type'] = s
        frames.append(df)
    func_df = pd.concat(frames, ignore_index=True)
    func_df['site'] = pd.to_numeric(func_df['site'], errors='coerce').astype('Int64')
    func_df['wild'] = func_df['wild'].astype(str).str.strip().str.upper()
    return func_df


def get_structure_coords(pdb_file):
    """Parse PDB structure and extract residue centroids."""
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, pdb_file)
    rows = []

    for model in structure:
        for chain in model:
            chain_id = chain.id
            histone = CHAIN_MAPPING.get(chain_id)
            if histone is None:
                continue
            for residue in chain:
                if residue.get_id()[0] != ' ':
                    continue
                resnum = residue.get_id()[1]
                resname = residue.get_resname().strip().upper()
                centroid = get_residue_centroid(residue)
                if centroid is not None:
                    rows.append({
                        'histone_type': histone,
                        'wild': resname,
                        'site': int(resnum),
                        'x': float(centroid[0]),
                        'y': float(centroid[1]),
                        'z': float(centroid[2])
                    })

    struct_df = pd.DataFrame(rows)
    print(f"Parsed structure: {len(struct_df)} residues with centroid coordinates")
    return struct_df


def unify_wild_names(df):
    """Convert single-letter residue codes to three-letter codes."""
    df = df.copy()
    df['wild'] = df['wild'].astype(str).str.strip().str.upper()
    df['wild'] = df['wild'].apply(lambda w: ONE2THREE.get(w, w) if len(w) == 1 else w)
    return df


def calc_Ms(centroid_coords):
    """Calculate Ms = mean(1/r_ij) for a set of centroid coordinates."""
    if len(centroid_coords) < 2:
        return np.nan
    distances = pdist(centroid_coords, metric='euclidean')
    distances = distances[distances > 1e-8]
    if len(distances) == 0:
        return np.nan
    return float(np.mean(1.0 / distances))


def calculate_p_value(centroid_coords1, centroid_coords2, n_iterations=10000, seed=42):
    """
    Perform permutation test to compute p-value for the difference in Ms between two sets.
    """
    random.seed(seed)
    np.random.seed(seed)

    Ms1 = calc_Ms(centroid_coords1)
    Ms2 = calc_Ms(centroid_coords2)
    observed_diff = Ms1 - Ms2

    if np.isnan(observed_diff):
        return np.nan, observed_diff, Ms1, Ms2

    all_coords = np.vstack([centroid_coords1, centroid_coords2])
    n1, n2 = len(centroid_coords1), len(centroid_coords2)
    count_extreme = 0

    for _ in range(n_iterations):
        indices = np.random.permutation(n1 + n2)
        pseudo_diff = calc_Ms(all_coords[indices[:n1]]) - calc_Ms(all_coords[indices[n1:]])
        if pseudo_diff >= observed_diff:
            count_extreme += 1

    p_value = (count_extreme + 1) / (n_iterations + 1)
    return p_value, observed_diff, Ms1, Ms2


def main():
    # 1. Load functional data
    func_df = read_functional_excel(EXCEL_PATH)

    if not os.path.exists(PDB_LOCAL):
        print(f"PDB file not found: {PDB_LOCAL}")
        return
    struct_df = get_structure_coords(PDB_LOCAL)

    # 2. Merge structural and functional data
    s = unify_wild_names(struct_df)
    f = unify_wild_names(func_df)
    merged = pd.merge(s, f, on=['histone_type', 'wild', 'site'], how='inner')

    if merged.empty:
        print("No matching residues found between structure and functional data")
        return

    # 3. Define two residue sets
    # Set 1: frequency > 0
    condition1 = merged['frequency'] > 0
    # Set 2: frequency >= 8
    condition2 = (merged['frequency'] >= 8)

    subset1 = merged[condition1]
    subset2 = merged[condition2]

    if len(subset1) < 2 or len(subset2) < 2:
        print("At least one set has fewer than 2 residues; cannot compute Ms")
        return

    coords1 = subset1[['x', 'y', 'z']].values
    coords2 = subset2[['x', 'y', 'z']].values

    # 4. Compute Ms and p-value
    p_value, diff, Ms1, Ms2 = calculate_p_value(coords1, coords2, BOOTSTRAP_ITERATIONS, SEED)

    # 5. Output results
    print("=" * 50)
    print("Ms comparison results (using residue centroids)")
    print("=" * 50)
    print(f"Set 1 (frequency > 0): {len(subset1)} residues, Ms = {Ms1:.6e}")
    print(f"Set 2 (frequency ≥ 8): {len(subset2)} residues, Ms = {Ms2:.6e}")
    print(f"Ms difference: {diff:.6e}")
    print(f"Permutation test p-value: {p_value:.6e}")

    if p_value < 0.01:
        print(f"Conclusion: Significant at 1% level (P = {p_value:.4f})")
    elif p_value < 0.05:
        print(f"Conclusion: Significant at 5% level (P = {p_value:.4f})")
    else:
        print(f"Conclusion: Not significant (P = {p_value:.4f})")


if __name__ == "__main__":
    main()