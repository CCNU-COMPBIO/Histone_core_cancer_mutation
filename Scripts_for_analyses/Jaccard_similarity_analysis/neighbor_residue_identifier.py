#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neighbor_residue_identifier.py

Analyze two types of neighboring residues for each residue in 3afa.pdb:
1. Spatial neighbors within specified cutoff (3/5/7 Å) based on minimum atom distance (not centroid).
2. Sequence neighbors within ±1 position along the same histone chain.
"""

import os
import numpy as np
import pandas as pd
from scipy.spatial import KDTree
from collections import defaultdict

# ==================== GLOBAL PARAMETERS ====================
EXCEL_PATH = "3afa.xlsx"
PDB_LOCAL = "3afa.pdb"
PDB_ID = "3afa"
CHAIN_MAPPING = {'C': 'H2A', 'D': 'H2B', 'A': 'H3', 'B': 'H4'}
CUTOFFS = [3, 5, 7]          # Supported distance cutoffs (Å)
SEQUENCE_WINDOW = 1            # Sequence window size (±1)

# Residue name mapping (one-letter to three-letter)
ONE2THREE = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}
THREE2ONE = {v: k for k, v in ONE2THREE.items()}

# ==================== HELPER FUNCTIONS ====================
def get_structure_coords(pdb_file):
    """
    Parse PDB file using Biopython.
    Returns:
        struct_df : DataFrame with residue info, indexed by residue ID.
        atom_coords : numpy array of all atom coordinates.
        atom_residue_indices : list mapping each atom to its residue index in struct_df.
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, pdb_file)
    rows = []
    all_atom_coords = []
    atom_residue_indices = []
    residue_index = 0

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
                atoms = list(residue.get_atoms())
                if not atoms:
                    continue

                coords = [atom.get_coord() for atom in atoms]
                centroid = np.mean(coords, axis=0)
                one_letter = THREE2ONE.get(resname, resname[0] if resname else 'X')

                rows.append({
                    'histone_type': histone,
                    'chain': chain_id,
                    'wild': resname,
                    'wild_one': one_letter,
                    'site': int(resnum),
                    'x': centroid[0],
                    'y': centroid[1],
                    'z': centroid[2]
                })

                for atom in atoms:
                    all_atom_coords.append(atom.get_coord())
                    atom_residue_indices.append(residue_index)

                residue_index += 1

    struct_df = pd.DataFrame(rows)
    print(f"Parsed PDB: {len(struct_df)} residues, {len(all_atom_coords)} atoms (for spatial neighborhood)")
    return struct_df, np.array(all_atom_coords), atom_residue_indices

def read_functional_excel(path):
    """Read functional annotation Excel file (sheets H2A, H2B, H3, H4)."""
    sheets = ['H2A', 'H2B', 'H3', 'H4']
    frames = []
    for s in sheets:
        df = pd.read_excel(path, sheet_name=s)
        df.columns = df.columns.str.strip()
        df['histone_type'] = s
        frames.append(df)
    func_df = pd.concat(frames, ignore_index=True)
    numeric_cols = ['GeoStab', 'MutPNI', 'percentage3', 'PPI', 'frequency']
    for col in numeric_cols:
        func_df[col] = pd.to_numeric(func_df[col], errors='coerce')
    func_df['site'] = pd.to_numeric(func_df['site'], errors='coerce').astype('Int64')
    func_df['wild'] = func_df['wild'].astype(str).str.strip().str.upper()
    return func_df

def unify_wild_names(df):
    """Standardize residue names (one-letter ↔ three-letter)."""
    df = df.copy()
    df['wild'] = df['wild'].astype(str).str.strip().str.upper()
    df['wild'] = df['wild'].apply(lambda w: ONE2THREE.get(w, w) if len(w)==1 else w)
    return df

def find_neighbors_within_cutoff(atom_coords, atom_residue_indices, cutoff):
    """
    Find spatial neighbors based on atom distances.
    Returns: dict {residue_index: list_of_neighbor_residue_indices} (unique).
    """
    kdtree = KDTree(atom_coords)
    neighbors_dict = defaultdict(set)

    for i, coord in enumerate(atom_coords):
        neighbors_idx = kdtree.query_ball_point(coord, cutoff)
        res_i = atom_residue_indices[i]
        for j in neighbors_idx:
            res_j = atom_residue_indices[j]
            if res_i != res_j:
                neighbors_dict[res_i].add(res_j)

    return {k: list(v) for k, v in neighbors_dict.items()}

def find_sequence_neighbors(struct_df, window):
    """
    Find sequence neighbors (same histone, site difference ≤ window).
    Returns: dict {residue_index: list_of_neighbor_residue_indices}.
    """
    histone_site_to_idx = defaultdict(list)
    for idx, row in struct_df.iterrows():
        key = (row['histone_type'], row['site'])
        histone_site_to_idx[key].append(idx)

    sequence_neighbors = defaultdict(list)
    for idx, row in struct_df.iterrows():
        histone = row['histone_type']
        site = row['site']
        for delta in range(-window, window + 1):
            if delta == 0:
                continue
            neighbor_site = site + delta
            key = (histone, neighbor_site)
            if key in histone_site_to_idx:
                for nb_idx in histone_site_to_idx[key]:
                    sequence_neighbors[idx].append(nb_idx)
    return sequence_neighbors

def format_residue_name(histone, wild_one, site):
    """Format residue name as 'H2A R29'."""
    return f"{histone} {wild_one}{site}"

def analyze_neighbors_for_condition(condition_df, struct_df, spatial_dict, seq_dict, cutoff):
    """
    For residues in condition_df, gather spatial and sequence neighbors.
    Returns DataFrame with neighbor lists and counts.
    """
    results = []
    for idx, row in condition_df.iterrows():
        histone = row['histone_type']
        wild_one = row['wild_one']
        site = row['site']

        spatial_neighbors = spatial_dict.get(idx, [])
        seq_neighbors = seq_dict.get(idx, [])

        spatial_list = []
        for nb_idx in spatial_neighbors:
            nb_row = struct_df.loc[nb_idx]
            spatial_list.append(format_residue_name(
                nb_row['histone_type'], nb_row['wild_one'], nb_row['site']))

        seq_list = []
        for nb_idx in seq_neighbors:
            nb_row = struct_df.loc[nb_idx]
            seq_list.append(format_residue_name(
                nb_row['histone_type'], nb_row['wild_one'], nb_row['site']))

        all_set = set(spatial_list) | set(seq_list)

        results.append({
            'Center_Residue': format_residue_name(histone, wild_one, site),
            f'Spatial_Neighbors_{cutoff}A': ", ".join(sorted(spatial_list)),
            'Sequence_Neighbors_±1': ", ".join(sorted(seq_list)),
            'All_Neighbors_Combined': ", ".join(sorted(all_set)),
            'Spatial_Neighbor_Count': len(spatial_list),
            'Sequence_Neighbor_Count': len(seq_list),
            'Total_Neighbor_Count': len(all_set)
        })
    return pd.DataFrame(results)

def run_analysis_for_cutoff(cutoff, struct_df, atom_coords, atom_residue_indices,
                            seq_dict, merged, output_file):
    """Run full analysis for a single cutoff and save results to Excel."""
    print(f"\n===== Analyzing cutoff {cutoff}Å (atom‑based minimum distance) =====")
    print("Finding spatial neighbors...")
    spatial_dict = find_neighbors_within_cutoff(atom_coords, atom_residue_indices, cutoff)

    # Define conditions based on functional metrics (GeoStab excluded)
    conditions = [
        ('PDI', (merged['MutPNI'] >= 1.0) & (merged['percentage3'] > 0)),
        ('PPI', merged['PPI'] >= 1.5),
        ('frequency', merged['frequency'] >= 8)
    ]

    all_results = []
    condition_names = []

    for cond_name, mask in conditions:
        cond_df = merged[mask].copy()
        if cond_df.empty:
            print(f"  Condition {cond_name}: no residues")
            continue
        print(f"  Condition {cond_name}: {len(cond_df)} residues")
        neighbors_df = analyze_neighbors_for_condition(
            cond_df, struct_df, spatial_dict, seq_dict, cutoff)
        if not neighbors_df.empty:
            all_results.append(neighbors_df)
            condition_names.append(cond_name)

    if not all_results:
        print(f"Cutoff {cutoff}Å produced no results")
        return

    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for cond_name, res_df in zip(condition_names, all_results):
            res_df = res_df.sort_values('Center_Residue')
            sheet_name = cond_name[:31]
            res_df.to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"  Results saved to: {output_file}")

# ==================== MAIN ====================
def main():
    try:
        from Bio.PDB import PDBParser
        from scipy.spatial import KDTree
        import openpyxl
    except ImportError as e:
        print(f"Missing dependency: {e}")
        return

    print("Reading Excel data...")
    func_df = read_functional_excel(EXCEL_PATH)

    if not os.path.exists(PDB_LOCAL):
        print(f"PDB file not found: {PDB_LOCAL}")
        return

    print("Parsing PDB structure...")
    struct_df, atom_coords, atom_residue_indices = get_structure_coords(PDB_LOCAL)
    if struct_df.empty:
        print("No valid residues in PDB")
        return

    struct_df = unify_wild_names(struct_df)
    func_df = unify_wild_names(func_df)
    merged = pd.merge(struct_df, func_df, on=['histone_type', 'wild', 'site'], how='inner')
    if merged.empty:
        print("No matching residues between structure and functional table")
        return
    print(f"Matched {len(merged)} residues")

    print(f"Finding sequence neighbors (window ±{SEQUENCE_WINDOW})...")
    seq_dict = find_sequence_neighbors(struct_df, SEQUENCE_WINDOW)

    for cutoff in CUTOFFS:
        out_file = f"Nuc_NeighborSets_AtomBased_Jaccard_{cutoff}A.xlsx"
        run_analysis_for_cutoff(cutoff, struct_df, atom_coords, atom_residue_indices,
                                seq_dict, merged, out_file)

    print("\nAll cutoff analyses completed.")

if __name__ == "__main__":
    main()