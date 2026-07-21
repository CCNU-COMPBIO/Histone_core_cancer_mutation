"""
Analyze histone binding sites from an Excel file.
Compute Jaccard similarity between two site sets (frequency >= 8 vs. PPI criteria)
with optional positional expansion, and calculate interface ratios.
"""

import pandas as pd
import numpy as np


def calculate_jaccard_similarity(set1, set2, offset=0):
    """
    Standard Jaccard: J(A,B) = |A ∩ B_exp| / |A ∪ B_exp|
    where B_exp is set2 expanded by ±offset positions.
    """
    def parse_site(s):
        h, pos = s.split("_")
        return h, int(pos)

    # Expand set2 (B)
    expanded_set2 = set()
    expansion_map = {}  # expanded site -> list of original sites
    for s in set2:
        h, p = parse_site(s)
        for dp in range(-offset, offset + 1):
            new_p = p + dp
            if new_p >= 1:
                new_s = f"{h}_{new_p}"
                expanded_set2.add(new_s)
                expansion_map.setdefault(new_s, []).append(s)

    intersection_set = set1.intersection(expanded_set2)
    union_set = set1.union(expanded_set2)

    inter = len(intersection_set)
    union = len(union_set)
    sim = inter / union if union > 0 else 0.0

    details = [(site, expansion_map[site]) for site in intersection_set]

    return sim, inter, union, len(expanded_set2), details


def count_interface_sites(file_path, site_set):
    """Count how many sites in site_set are present in the 'percentage2' binding interface sheets."""
    histone_names = ['H2A', 'H2B', 'H3', 'H4']
    interface_count = 0
    for histone in histone_names:
        try:
            df = pd.read_excel(file_path, sheet_name=histone)
            df['percentage2'] = pd.to_numeric(df['percentage2'], errors='coerce')
            binding_sites = df[df['percentage2'] > 0]
            for _, row in binding_sites.iterrows():
                site_name = f"{histone}_{int(row['site'])}"
                if site_name in site_set:
                    interface_count += 1
        except Exception as e:
            print(f"Error processing {histone} interface data: {e}")
    return interface_count


def process_histone_data(file_path):
    """
    Extract site sets from each histone sheet:
      - Set A: sites with frequency >= 8
      - Set B: sites with |PPI| < 0.5
    """
    histone_names = ['H2A', 'H2B', 'H3', 'H4']
    ppi_set = set()
    freq_set = set()

    for histone in histone_names:
        try:
            df = pd.read_excel(file_path, sheet_name=histone)
            df['DDG_PPI_avg'] = pd.to_numeric(df['DDG_PPI_avg'], errors='coerce')
            df['frequency'] = pd.to_numeric(df['frequency'], errors='coerce')

            # PPI condition: |PPI| < 0.5
            ppi_high = df[df['DDG_PPI_avg'].abs() < 0.5]
            for _, row in ppi_high.iterrows():
                ppi_set.add(f"{histone}_{int(row['site'])}")

            # Frequency condition: >= 8
            freq_high = df[df['frequency'] >= 8]
            for _, row in freq_high.iterrows():
                freq_set.add(f"{histone}_{int(row['site'])}")

        except Exception as e:
            print(f"Error processing {histone} data: {e}")

    return ppi_set, freq_set


def main():
    file_path = "AlanineScan_DDG_Nucleosome_Averages.xlsx"

    ppi_set, freq_set = process_histone_data(file_path)
    set_A = freq_set   # Frequency >= 8
    set_B = ppi_set    # PPI criteria

    combined_set = set_A.union(set_B)

    print("=" * 50)
    print(f"Set A (Frequency ≥ 8): {len(set_A)} sites")
    print(f"Set B (PPI criteria):  {len(set_B)} sites")
    print(f"Union (original):      {len(combined_set)} sites")

    # Interface ratios
    inter_A = count_interface_sites(file_path, set_A)
    inter_B = count_interface_sites(file_path, set_B)
    inter_union = count_interface_sites(file_path, combined_set)
    print(f"\nInterface ratio (A): {inter_A/len(set_A):.2%}")
    print(f"Interface ratio (B): {inter_B/len(set_B):.2%}")
    print(f"Interface ratio (Union): {inter_union/len(combined_set):.2%}")

    print("\nJaccard similarity J(A, B_exp) with B expansion:")
    for offset, label in [(0, "exact"), (1, "±1"), (2, "±2")]:
        sim, inter, union, exp_size, details = calculate_jaccard_similarity(set_A, set_B, offset)
        print(f"  Offset {label}: J = {sim:.3f}  (intersection = {inter}, |B_exp| = {exp_size}, |A∪B_exp| = {union})")
        if offset == 2:
            print("  Detailed matches (A site -> original B sites):")
            for site, orig_list in details:
                print(f"    {site} -> {', '.join(orig_list)}")


if __name__ == "__main__":
    main()