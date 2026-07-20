import pandas as pd
import numpy as np


def calculate_jaccard_similarity(set1, set2, offset=0):
    """
    标准 Jaccard: J(A,B) = |A ∩ B_exp| / |A ∪ B_exp|
    其中 B_exp 为 set2 每个位点扩展 ±offset 后的集合
    """
    def parse_site(s):
        h, pos = s.split("_")
        return h, int(pos)

    # 扩展 set2 (B)
    expanded_set2 = set()
    expansion_map = {}  # 扩展位点 -> 原始位点列表
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

    # 详情：每个重合的 A 位点对应哪些原始 B 位点
    details = [(site, expansion_map[site]) for site in intersection_set]

    return sim, inter, union, len(expanded_set2), details


def count_interface_sites(file_path, site_set):
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
            print(f"处理 {histone} 结合界面统计出错: {e}")
    return interface_count


def process_histone_data(file_path):
    histone_names = ['H2A', 'H2B', 'H3', 'H4']
    ppi_high_set = set()
    frequency_top_set = set()

    for histone in histone_names:
        try:
            df = pd.read_excel(file_path, sheet_name=histone)
            df['PPI'] = pd.to_numeric(df['PPI'], errors='coerce')
            df['frequency'] = pd.to_numeric(df['frequency'], errors='coerce')

            # PPI 条件（根据您的规则）
            ppi_high = df[(df['PPI'].abs() < 0.5)]
            for _, row in ppi_high.iterrows():
                ppi_high_set.add(f"{histone}_{int(row['site'])}")

            # Frequency ≥ 8
            freq_high = df[(df['frequency'] >= 8)]
            for _, row in freq_high.iterrows():
                frequency_top_set.add(f"{histone}_{int(row['site'])}")

        except Exception as e:
            print(f"处理 {histone} 数据出错: {e}")

    return ppi_high_set, frequency_top_set


def main():
    file_path = "3afa.xlsx"

    ppi_set, freq_set = process_histone_data(file_path)   # 注意返回值顺序：ppi_set, freq_set
    # 但我们要把 freq_set 作为 Set A，ppi_set 作为 Set B
    set_A = freq_set   # Frequency ≥ 8
    set_B = ppi_set    # PPI 相关

    combined_set = set_A.union(set_B)

    print("=" * 50)
    print("Set A (Frequency ≥ 8):", len(set_A))
    print("Set B (PPI criteria):", len(set_B))
    print("Union (original):", len(combined_set))

    # 界面占比
    inter_A = count_interface_sites(file_path, set_A)
    inter_B = count_interface_sites(file_path, set_B)
    inter_union = count_interface_sites(file_path, combined_set)
    print(f"Set A interface ratio: {inter_A/len(set_A):.2%}")
    print(f"Set B interface ratio: {inter_B/len(set_B):.2%}")
    print(f"Union interface ratio: {inter_union/len(combined_set):.2%}")

    print("\nJaccard similarity J(A, B_exp) with expansion of B:")
    for offset, label in [(0, "exact"), (1, "±1"), (2, "±2")]:
        sim, inter, union, exp_size, details = calculate_jaccard_similarity(set_A, set_B, offset)
        print(f"  Offset {label}: J = {sim:.3f} (重合个数 = {inter}, 扩展后 |B_exp| = {exp_size}, |A∪B_exp| = {union})")
        if offset == 2:
            print("  匹配详情 (A位点 -> 对应原始B位点):")
            for site, orig_list in details:
                print(f"    {site} -> {', '.join(orig_list)}")

if __name__ == "__main__":
    main()