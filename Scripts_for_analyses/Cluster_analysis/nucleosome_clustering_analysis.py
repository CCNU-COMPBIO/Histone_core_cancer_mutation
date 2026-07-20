#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nucleosome_cluster_ms.py
Based on Ms (mean(1/r_ij)) + randomization test and average-linkage subclustering.
Applies to 3afa structure and functional Excel with four sheets (H2A,H2B,H3,H4).
Output: clusters.csv, ms_random_test.txt, clusters_3d.png
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt

# Optional: biopython for PDB parsing and download
try:
    from Bio.PDB import PDBParser, PDBList
except Exception:
    PDBParser = None
    PDBList = None

# ------------- Parameters (modifiable) -------------
EXCEL_PATH = "AlanineScan_DDG_Nucleosome_Averages.xlsx"  # Your Excel file
PDB_LOCAL = "3afa.pdb"  # Local PDB file path if available
PDB_ID = "3afa"  # PDB ID for download if local not found
CHAIN_MAPPING = {'C': 'H2A', 'D': 'H2B', 'A': 'H3', 'B': 'H4'}
N_RANDOM = 1000
MIN_CLUSTER_SIZE = 2
OUT_PREFIX = "nucleosome_clusters"
# ------------------------------------------------

# Some residue name mappings (one-letter -> three-letter)
ONE2THREE = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}


def get_structure_coords(pdb_file):
    """
    Parse PDB with Biopython, extract centroid coordinates (geometric center of all atoms) for each residue.
    Modified: uses residue centroids instead of C-alpha coordinates.
    """
    if PDBParser is None:
        raise ImportError("Biopython not installed, please pip install biopython")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(PDB_ID, pdb_file)
    rows = []
    for model in structure:
        for chain in model:
            chain_id = chain.id
            histone = CHAIN_MAPPING.get(chain_id, None)
            if histone is None:
                continue
            for residue in chain:
                if residue.get_id()[0] != ' ':  # skip heteroatoms
                    continue
                resnum = residue.get_id()[1]
                resname = residue.get_resname().strip().upper()

                # Collect all atom coordinates in this residue
                atom_coords = []
                for atom in residue:
                    coord = atom.get_coord()
                    atom_coords.append(coord)

                if atom_coords:
                    atom_coords_array = np.array(atom_coords)
                    centroid = np.mean(atom_coords_array, axis=0)

                    rows.append({
                        'histone_type': histone,
                        'chain': chain_id,
                        'wild': resname,
                        'site': int(resnum),
                        'x': float(centroid[0]),
                        'y': float(centroid[1]),
                        'z': float(centroid[2]),
                        'all_atom_coords': atom_coords_array
                    })
    df = pd.DataFrame(rows)
    print(f"Parsed structure: {len(df)} residues with centroid coordinates from chains {sorted(df['chain'].unique())}")
    return df


def calculate_centroid_distance_matrix(centroid_coords):
    """
    Calculate Euclidean distance matrix based on residue centroid coordinates.
    """
    if centroid_coords.shape[0] == 0:
        return np.array([])
    dist_vector = pdist(centroid_coords, metric='euclidean')
    dist_matrix = squareform(dist_vector)
    return dist_matrix


def calc_Ms_centroid_based(centroid_coords):
    """
    Calculate Ms = mean(1/r_ij) based on centroid Euclidean distances.
    """
    if centroid_coords.shape[0] < 2:
        return np.nan
    dist_matrix = calculate_centroid_distance_matrix(centroid_coords)
    n = dist_matrix.shape[0]
    dist_list = []
    for i in range(n):
        for j in range(i + 1, n):
            dist = dist_matrix[i, j]
            if dist > 1e-8:
                dist_list.append(dist)
    if not dist_list:
        return np.nan
    inv = 1.0 / np.array(dist_list)
    return float(np.mean(inv))


def random_Ms_distribution_centroid_based(all_centroid_coords, k, n_random=N_RANDOM):
    """
    Generate random Ms distribution by sampling k residues from all centroid coordinates.
    """
    n_all = all_centroid_coords.shape[0]
    if k < 2 or n_all < k:
        return np.array([])
    Ms_rand = []
    for _ in range(n_random):
        idx = np.random.choice(n_all, k, replace=False)
        sampled_coords = all_centroid_coords[idx]
        Ms_rand.append(calc_Ms_centroid_based(sampled_coords))
    return np.array(Ms_rand)


def hierarchical_subclusters_centroid_based(centroid_coords, threshold=None, min_size=MIN_CLUSTER_SIZE):
    """
    Hierarchical clustering (average-linkage) based on centroid distances.
    """
    if centroid_coords.shape[0] == 0:
        return np.array([], dtype=int)
    if centroid_coords.shape[0] == 1:
        return np.array([0], dtype=int)

    if threshold is None:
        threshold = 15.0
        print(f"Warning: using default threshold {threshold} Å for centroid distance clustering")

    dist_matrix = calculate_centroid_distance_matrix(centroid_coords)
    Z = linkage(squareform(dist_matrix), method='average', metric='precomputed')
    labels = fcluster(Z, t=threshold, criterion='distance')

    # Handle small clusters
    lab_series = pd.Series(labels)
    counts = lab_series.value_counts()
    small = counts[counts < min_size].index.tolist()
    lab_series[lab_series.isin(small)] = 0

    # Renumber non-zero clusters
    nonzero = sorted(lab_series[lab_series > 0].unique())
    mapping = {old: new for new, old in enumerate(nonzero, start=1)}
    lab_series = lab_series.replace(mapping).astype(int)

    return lab_series.values


def calculate_dynamic_threshold(atom_coords_list):
    """
    Compute threshold based on the mean of all pair-wise maximum interatomic distances.
    """
    if len(atom_coords_list) < 2:
        print(f"Warning: insufficient residues, using default 15.0 Å")
        return 15.0

    n_residues = len(atom_coords_list)
    all_pair_distances = []

    for i in range(n_residues):
        for j in range(i + 1, n_residues):
            dist = max_interatomic_distance(atom_coords_list[i], atom_coords_list[j])
            all_pair_distances.append(dist)

    if len(all_pair_distances) == 0:
        print(f"Warning: cannot compute distances, using default 15.0 Å")
        return 15.0

    mean_distance = np.mean(all_pair_distances)
    dynamic_threshold = mean_distance / 2.0

    print(f"Dynamic threshold (based on all-atom distances): residues={n_residues}, pairs={len(all_pair_distances)}, mean_dist={mean_distance:.2f} Å, threshold={dynamic_threshold:.2f} Å")
    return dynamic_threshold


import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage as scipy_linkage
from scipy.spatial.distance import squareform as scipy_squareform
import matplotlib.pyplot as plt

THREE2ONE = {v: k for k, v in ONE2THREE.items()}


def plot_dendrogram_with_residue_labels(subset_df, centroid_coords, threshold, out_png):
    """
    Draw dendrogram with residue labels (e.g., H3 E97).
    """
    sns.set_style("white")
    plt.rcParams["font.family"] = "Arial"

    labels = []
    for _, row in subset_df.iterrows():
        histone = row["histone_type"]
        wild3 = row["wild"]
        wild1 = THREE2ONE.get(wild3, wild3)
        site = int(row["site"])
        labels.append(f"{histone} {wild1}{site}")

    dist_matrix = calculate_centroid_distance_matrix(centroid_coords)
    condensed_dist = scipy_squareform(dist_matrix)
    Z = scipy_linkage(condensed_dist, method="average")

    plt.figure(figsize=(20, 8))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=14)

    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)
        spine.set_color('black')

    plt.axhline(
        y=threshold,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Threshold = {threshold:.2f} Å"
    )

    plt.title(r'Hierarchical clustering dendrogram of residues with mutation hotspots',
              fontsize=24, pad=20, fontweight='bold')
    plt.ylabel("Centroid distance (Å)", fontsize=22, fontweight="bold")
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=22)

    plt.tight_layout()
    plt.savefig(out_png, dpi=600, format='tiff', pil_kwargs={'compression': 'tiff_lzw'})
    plt.close()

    return Z


def read_functional_excel(path):
    """Read the four histone sheets and unify column names (use 'DDG_PPI_avg' for energetic column)."""
    sheets = ['H2A', 'H2B', 'H3', 'H4']
    frames = []
    for s in sheets:
        try:
            df = pd.read_excel(path, sheet_name=s)
            df = df.copy()
            df.columns = df.columns.str.strip()
            df['histone_type'] = s
            # Keep DDG_PPI_avg if present
            if 'DDG_PPI_avg' in df.columns:
                df['DDG_PPI_avg'] = df['DDG_PPI_avg']
            frames.append(df)
            print(f"Loaded sheet {s}: {len(df)} rows; cols: {df.columns.tolist()}")
        except Exception as e:
            print(f"Warning: cannot read sheet {s} from {path}: {e}")
    if not frames:
        raise FileNotFoundError("No H2A/H2B/H3/H4 sheets found in Excel.")
    func_df = pd.concat(frames, ignore_index=True)
    if 'wild' not in func_df.columns or 'site' not in func_df.columns:
        raise KeyError("Excel sheet must contain 'wild' and 'site' columns.")
    func_df['site'] = pd.to_numeric(func_df['site'], errors='coerce').astype('Int64')
    func_df['wild'] = func_df['wild'].astype(str).str.strip().str.upper()
    return func_df


def max_interatomic_distance(res1_coords, res2_coords):
    """
    Compute the maximum distance between any atom in res1 and any atom in res2.
    """
    diff = res1_coords[:, np.newaxis, :] - res2_coords[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff**2, axis=2))
    return np.max(distances)


def unify_wild_names(df):
    """Convert single-letter amino acid codes to three-letter codes."""
    def conv(w):
        if len(w) == 1:
            return ONE2THREE.get(w, w)
        return w
    df = df.copy()
    df['wild'] = df['wild'].astype(str).str.strip().str.upper()
    df['wild'] = df['wild'].apply(conv)
    return df


def merge_struct_func(struct_df, func_df):
    """Merge structural and functional data on histone_type, wild, site."""
    s = struct_df.copy()
    f = func_df.copy()
    s = unify_wild_names(s)
    f = unify_wild_names(f)
    merged = pd.merge(s, f, on=['histone_type', 'wild', 'site'], how='inner', suffixes=('', '_f'))
    print(f"Merged structural+functional: {len(merged)} rows")
    return merged


def empirical_pvalue(ms_obs, ms_rand):
    """Empirical p-value: P(Ms_rand >= Ms_obs)."""
    if len(ms_rand) == 0 or np.isnan(ms_obs):
        return np.nan
    return (np.sum(ms_rand >= ms_obs) + 1) / (len(ms_rand) + 1)


def main():
    # 1. Read functional table
    func_df = read_functional_excel(EXCEL_PATH)

    # 2. Read structure
    pdb_file = PDB_LOCAL
    if not os.path.exists(pdb_file):
        if PDBList is not None:
            print(f"{pdb_file} not found; attempting download of {PDB_ID}...")
            pdbl = PDBList()
            pdb_file = pdbl.retrieve_pdb_file(PDB_ID, pdir='.', file_format='pdb')
            print("Downloaded PDB:", pdb_file)
        else:
            raise FileNotFoundError(f"{pdb_file} not found and Biopython PDBList not available.")

    struct_df = get_structure_coords(pdb_file)

    # 3. Merge
    merged = merge_struct_func(struct_df, func_df)
    if merged.empty:
        print("Warning: merged is empty — no matching residues found.")
        return

    # 4. Filter by frequency (keep frequency >= 8)
    merged['frequency'] = pd.to_numeric(merged['frequency'], errors='coerce')
    condition = (merged['frequency'] >= 8)
    subset = merged[condition].copy()

    print(f"Total merged residues: {len(merged)}; selected {len(subset)} with frequency >= 8")
    if subset.empty:
        print("No high-frequency residues found; exiting.")
        return

    # Extract coordinates
    atom_coords_list = subset['all_atom_coords'].tolist()
    centroid_coords = subset[['x', 'y', 'z']].values

    # 5. Ms calculation
    Ms_obs = calc_Ms_centroid_based(centroid_coords)

    # Randomization test using all centroid coordinates
    all_centroid_coords = struct_df[['x', 'y', 'z']].values
    ms_rand = random_Ms_distribution_centroid_based(all_centroid_coords, centroid_coords.shape[0], n_random=N_RANDOM)
    p_emp = empirical_pvalue(Ms_obs, ms_rand)
    zscore = (Ms_obs - np.nanmean(ms_rand)) / (np.nanstd(ms_rand, ddof=1) if np.nanstd(ms_rand, ddof=1) > 0 else np.nan)

    try:
        mw_p = mannwhitneyu(np.array([Ms_obs]), ms_rand, alternative='greater').pvalue
    except Exception:
        mw_p = np.nan

    print(f"Ms_obs = {Ms_obs:.6e}, ms_rand_mean = {np.nanmean(ms_rand):.6e}, empirical p = {p_emp:.4f}, z = {zscore:.3f}, mannwhitney p = {mw_p}")

    # 6. Subclustering with dynamic threshold
    dynamic_threshold = calculate_dynamic_threshold(atom_coords_list)
    subcluster_labels = hierarchical_subclusters_centroid_based(centroid_coords, threshold=dynamic_threshold,
                                                                min_size=MIN_CLUSTER_SIZE)
    subset['subcluster'] = subcluster_labels

    from sklearn.metrics import silhouette_score
    valid = subset['subcluster'] > 0
    coords_valid = centroid_coords[valid.values]
    labels_valid = subset.loc[valid, 'subcluster'].values
    if len(np.unique(labels_valid)) > 1 and coords_valid.shape[0] > len(np.unique(labels_valid)):
        sil = silhouette_score(coords_valid, labels_valid, metric='euclidean')
        print(f"Silhouette score (non-noise clusters): {sil:.3f}")
    else:
        print("Too few clusters or samples to compute Silhouette score")

    from sklearn.metrics import davies_bouldin_score, calinski_harabasz_score
    if len(np.unique(labels_valid)) > 1:
        db = davies_bouldin_score(coords_valid, labels_valid)
        ch = calinski_harabasz_score(coords_valid, labels_valid)
        print(f"Davies-Bouldin index: {db:.3f} (lower is better)")
        print(f"Calinski-Harabasz index: {ch:.3f} (higher is better)")

    # 7. Dendrogram
    Z = plot_dendrogram_with_residue_labels(
        subset_df=subset,
        centroid_coords=centroid_coords,
        threshold=dynamic_threshold,
        out_png=f"{OUT_PREFIX}_dendrogram.tif"
    )

    labels = fcluster(Z, t=dynamic_threshold, criterion='distance')
    labels = pd.Series(labels)
    counts = labels.value_counts()
    small = counts[counts < MIN_CLUSTER_SIZE].index
    labels[labels.isin(small)] = 0
    n_clusters = labels[labels > 0].nunique()
    print(f"Number of subclusters at threshold {dynamic_threshold:.2f} Å: {n_clusters}")

    # 8. Save results
    out_csv = f"{OUT_PREFIX}.csv"
    subset.to_csv(out_csv, index=False)
    print(f"Saved results to {out_csv}")

    with open(f"{OUT_PREFIX}_ms_test.txt", "w") as fh:
        fh.write(f"Ms_obs = {Ms_obs}\n")
        fh.write(f"ms_rand_mean = {np.nanmean(ms_rand)}\n")
        fh.write(f"ms_rand_std = {np.nanstd(ms_rand, ddof=1)}\n")
        fh.write(f"empirical_p = {p_emp}\n")
        fh.write(f"zscore = {zscore}\n")
        fh.write(f"mannwhitney_p = {mw_p}\n")
        fh.write(f"N_random = {N_RANDOM}\n")
        fh.write(f"Dynamic_threshold_used = {dynamic_threshold}\n")
        fh.write(f"Distance_metric = Euclidean_distance_based_on_residue_centroids\n")
    print(f"Saved Ms random test summary to {OUT_PREFIX}_ms_test.txt")


if __name__ == "__main__":
    main()