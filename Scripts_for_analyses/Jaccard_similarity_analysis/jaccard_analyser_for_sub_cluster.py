#!/usr/bin/env python3
# encoding: utf-8
"""
Compute Jaccard similarities between subclusters with expandable neighborhoods.
Includes overlap (intersection) details for each comparison.
"""

import pandas as pd
import os

# ==================== 1. HARD-CODED SUBCLUSTER RESIDUES ====================
freq_subclusters = {
    1: ["H2A K74", "H2A G37", "H2B G53"],
    2: ["H2A R29", "H2A G37", "H2B E35", "H2B S36", "H2B D68", "H2B F70", "H2B E71"],
    3: ["H2A S19", "H2A E56", "H2B E113"],
    4: ["H3 R40", "H3 E50", "H3 R52"],
    5: ["H4 D68", "H2B E76", "H2B T88", "H2B E93", "H2B G104"],
    6: ["H3 R72", "H3 E73", "H3 D77", "H3 D81", "H3 R83"],
    7: ["H3 E97", "H3 E105", "H3 D106", "H3 A127", "H3 R131", "H3 E133", "H3 R134"]
}

ppi_subclusters = {
    1: ["H2A L63", "H2A I102", "H2B L45", "H2B I61", "H2B M62", "H2B F65"],
    2: ["H2A F25", "H2A L34", "H2A Y39", "H2A L51", "H2A L55", "H2B Y37", "H2B F70", "H2B I89"],
    3: ["H3 D106", "H3 L126"],
    4: ["H3 I51", "H3 Y54", "H3 F104", "H3 I119", "H4 I34", "H4 L37", "H4 R40"],
    5: ["H4 D68", "H4 Y72", "H4 D85", "H4 Y88", "H2B E76"],
    6: ["H3 F67", "H3 L92", "H3 Y99", "H3 L100", "H4 F61", "H4 L97", "H4 Y98"]
}

pdi_subclusters = {
    1: ["H3 P38", "H3 Y41", "H3 R49", "H3 R53"],
    2: ["H3 K115", "H3 R116", "H3 V117", "H3 T118", "H3 K122"],
    3: ["H3 P66", "H4 R36", "H4 Y51"],
    4: ["H3 R83", "H4 K77", "H4 R78", "H4 K79", "H4 T80", "H2B K85", "H2B R86", "H2B T88", "H2B T90", "H2B R92"],
    5: ["H2A K74", "H2A R77"],
    6: ["H2A R29", "H2B Y40", "H2B Y42", "H2B K43"]
}

# ==================== 2. HELPER FUNCTIONS ====================
def clean_residue_string(res_str):
    if pd.isna(res_str):
        return ""
    return str(res_str).strip()

def build_residue_neighbors_map_dynamic(file_path, sheet_name, cutoff_label):
    """
    Build a dictionary mapping center residue -> (spatial_list, seq_list).
    cutoff_label: e.g., '3A', '5A', '7A' to construct the spatial column name.
    """
    df = pd.read_excel(file_path, sheet_name=sheet_name)
    neighbors_map = {}
    spatial_col = f'Spatial_Neighbors_{cutoff_label}'
    seq_col = 'Sequence_Neighbors_±1'
    for _, row in df.iterrows():
        center = clean_residue_string(row['Center_Residue'])
        if not center:
            continue
        spatial_str = clean_residue_string(row.get(spatial_col, ''))
        spatial_list = [n.strip() for n in spatial_str.split(',')] if spatial_str else []
        seq_str = clean_residue_string(row.get(seq_col, ''))
        seq_list = [n.strip() for n in seq_str.split(',')] if seq_str else []
        neighbors_map[center] = (spatial_list, seq_list)
    return neighbors_map

def get_center_set(subcluster_dict, subcluster_id):
    return set(subcluster_dict[subcluster_id])

def get_expanded_set(center_set, neighbors_map, use_spatial=True):
    expanded = set(center_set)
    for residue in center_set:
        if residue in neighbors_map:
            spatial_list, seq_list = neighbors_map[residue]
            if use_spatial:
                expanded.update(spatial_list)
            else:
                expanded.update(seq_list)
    expanded.discard("")
    return expanded

def jaccard_similarity(set1, set2):
    if not set1 and not set2:
        return 0.0
    inter = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return inter / union if union != 0 else 0.0

def asymmetric_jaccard(expanded_set, ref_center_set):
    """Asymmetric Jaccard: intersection / (expanded_set ∪ ref_center_set)"""
    union = expanded_set.union(ref_center_set)
    if not union:
        return 0.0
    inter = len(expanded_set.intersection(ref_center_set))
    return inter / len(union)

def intersection_string(set1, set2):
    inter = set1.intersection(set2)
    if not inter:
        return ""
    return ", ".join(sorted(inter))

def compute_symmetric_matrix(exp_subclusters, ref_subclusters):
    exp_ids = sorted(exp_subclusters.keys())
    ref_ids = sorted(ref_subclusters.keys())
    exp_center_sets = {sid: get_center_set(exp_subclusters, sid) for sid in exp_ids}
    ref_center_sets = {sid: get_center_set(ref_subclusters, sid) for sid in ref_ids}
    matrix = {}
    inter_matrix = {}
    for eid in exp_ids:
        matrix[eid] = {}
        inter_matrix[eid] = {}
        for rid in ref_ids:
            set1 = exp_center_sets[eid]
            set2 = ref_center_sets[rid]
            matrix[eid][rid] = jaccard_similarity(set1, set2)
            inter_matrix[eid][rid] = intersection_string(set1, set2)
    return matrix, inter_matrix, exp_center_sets, ref_center_sets

def compute_sequence_asymmetric_matrix(exp_subclusters, ref_subclusters, seq_neighbors_map):
    exp_ids = sorted(exp_subclusters.keys())
    ref_ids = sorted(ref_subclusters.keys())
    exp_center_sets = {sid: get_center_set(exp_subclusters, sid) for sid in exp_ids}
    ref_center_sets = {sid: get_center_set(ref_subclusters, sid) for sid in ref_ids}
    exp_seq_sets = {sid: get_expanded_set(exp_center_sets[sid], seq_neighbors_map, use_spatial=False)
                    for sid in exp_ids}
    matrix = {}
    inter_matrix = {}
    for eid in exp_ids:
        matrix[eid] = {}
        inter_matrix[eid] = {}
        for rid in ref_ids:
            matrix[eid][rid] = asymmetric_jaccard(exp_seq_sets[eid], ref_center_sets[rid])
            inter_matrix[eid][rid] = intersection_string(exp_seq_sets[eid], ref_center_sets[rid])
    return matrix, inter_matrix

def compute_spatial_asymmetric_matrix(exp_subclusters, ref_subclusters, spatial_neighbors_map):
    exp_ids = sorted(exp_subclusters.keys())
    ref_ids = sorted(ref_subclusters.keys())
    exp_center_sets = {sid: get_center_set(exp_subclusters, sid) for sid in exp_ids}
    ref_center_sets = {sid: get_center_set(ref_subclusters, sid) for sid in ref_ids}
    exp_spatial_sets = {sid: get_expanded_set(exp_center_sets[sid], spatial_neighbors_map, use_spatial=True)
                        for sid in exp_ids}
    matrix = {}
    inter_matrix = {}
    for eid in exp_ids:
        matrix[eid] = {}
        inter_matrix[eid] = {}
        for rid in ref_ids:
            matrix[eid][rid] = asymmetric_jaccard(exp_spatial_sets[eid], ref_center_sets[rid])
            inter_matrix[eid][rid] = intersection_string(exp_spatial_sets[eid], ref_center_sets[rid])
    return matrix, inter_matrix

def matrix_to_dataframe(matrix, exp_name, ref_name):
    if not matrix:
        return pd.DataFrame()
    exp_ids = sorted(matrix.keys())
    ref_ids = sorted(matrix[exp_ids[0]].keys())
    df = pd.DataFrame(matrix).T
    df.index = [f"{exp_name}_{i}" for i in df.index]
    df.columns = [f"{ref_name}_{i}" for i in df.columns]
    return df

# ==================== 3. MAIN ====================
def main():
    # File list matching the output from the neighbor identification script
    files = [
        ("Nuc_NeighborSets_AtomBased_Jaccard_3A.xlsx", "3A"),
        ("Nuc_NeighborSets_AtomBased_Jaccard_5A.xlsx", "5A"),
        ("Nuc_NeighborSets_AtomBased_Jaccard_7A.xlsx", "7A")
    ]

    # Analysis configurations: (exp_dict, ref_dict, exp_name, ref_name, neighbor_sheet)
    # ref_name is changed to 'hotspot' to reflect frequency subclusters
    analyses = [
        (pdi_subclusters, freq_subclusters, "PDI", "hotspot", "PDI"),
        (ppi_subclusters, freq_subclusters, "PPI", "hotspot", "PPI")
    ]

    for exp_dict, ref_dict, exp_name, ref_name, neighbor_sheet in analyses:
        print(f"\nProcessing {exp_name} vs {ref_name}...")

        # 1. Symmetric matrix (center residues only)
        sym_matrix, sym_inter_matrix, exp_centers, ref_centers = compute_symmetric_matrix(exp_dict, ref_dict)
        sym_df = matrix_to_dataframe(sym_matrix, exp_name, ref_name)
        sym_inter_df = matrix_to_dataframe(sym_inter_matrix, exp_name, ref_name)

        # 2. Sequence asymmetric matrix
        seq_neighbors_map = None
        for file_path, label in files:
            if os.path.exists(file_path):
                try:
                    seq_neighbors_map = build_residue_neighbors_map_dynamic(file_path, neighbor_sheet, label)
                    print(f"  Sequence neighbor info loaded from {file_path}")
                    break
                except:
                    continue
        if seq_neighbors_map is None:
            print(f"  Failed to load sequence neighbor info, skipping {exp_name}")
            continue

        seq_matrix, seq_inter_matrix = compute_sequence_asymmetric_matrix(exp_dict, ref_dict, seq_neighbors_map)
        seq_df = matrix_to_dataframe(seq_matrix, exp_name, ref_name)
        seq_inter_df = matrix_to_dataframe(seq_inter_matrix, exp_name, ref_name)

        # 3. Statistics (cluster sizes)
        stats = []
        for sid in sorted(exp_centers.keys()):
            stats.append({'Type': exp_name, 'Subcluster': sid, 'Centers_Only': len(exp_centers[sid])})
        for sid in sorted(ref_centers.keys()):
            stats.append({'Type': ref_name, 'Subcluster': sid, 'Centers_Only': len(ref_centers[sid])})
        stats_df = pd.DataFrame(stats)

        # 4. Spatial asymmetric matrices for each cutoff
        spatial_dfs = {}
        spatial_inter_dfs = {}
        for file_path, label in files:
            if not os.path.exists(file_path):
                print(f"  File not found, skipping spatial matrix: {file_path}")
                continue
            try:
                spatial_map = build_residue_neighbors_map_dynamic(file_path, neighbor_sheet, label)
                spatial_matrix, spatial_inter_matrix = compute_spatial_asymmetric_matrix(exp_dict, ref_dict, spatial_map)
                spatial_df = matrix_to_dataframe(spatial_matrix, exp_name, ref_name)
                spatial_inter_df = matrix_to_dataframe(spatial_inter_matrix, exp_name, ref_name)
                spatial_dfs[label] = spatial_df
                spatial_inter_dfs[label] = spatial_inter_df
                print(f"  Computed spatial matrix for {label}")
            except Exception as e:
                print(f"  Failed to compute spatial matrix for {label}: {e}")

        if not spatial_dfs:
            print(f"  No spatial matrices computed, skipping {exp_name}")
            continue

        # Write to Excel with descriptive sheet names
        output_file = f"{exp_name}_vs_{ref_name}_jaccard_summary.xlsx"
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            sym_df.to_excel(writer, sheet_name="jaccard_vs_centers")
            sym_inter_df.to_excel(writer, sheet_name="overlap_vs_centers")
            for label, df in spatial_dfs.items():
                sheet_name = f"jaccard_vs_spatial_{label}"
                if len(sheet_name) > 31:
                    sheet_name = sheet_name[:31]
                df.to_excel(writer, sheet_name=sheet_name)
                inter_sheet = f"overlap_vs_spatial_{label}"
                if len(inter_sheet) > 31:
                    inter_sheet = inter_sheet[:31]
                spatial_inter_dfs[label].to_excel(writer, sheet_name=inter_sheet)
            seq_df.to_excel(writer, sheet_name="jaccard_vs_seq_pm1")
            seq_inter_df.to_excel(writer, sheet_name="overlap_vs_seq_pm1")

        print(f"  Results saved to {output_file}")

    print("\nAll summary calculations completed.")

if __name__ == "__main__":
    main()