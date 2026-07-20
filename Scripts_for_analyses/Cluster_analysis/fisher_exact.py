"""
Fisher's exact test for association between mutation hotspots and energetically important sites.
Input: Excel file with four histone sheets (H2A, H2B, H3, H4).
"""

import pandas as pd
import scipy.stats as stats

EXCEL_FILE = "AlanineScan_DDG_Nucleosome_Averages.xlsx"
SHEET_NAMES = ["H2A", "H2B", "H3", "H4"]

total_a = 0
total_b = 0
total_c = 0
total_d = 0

print(f"{'Histone':<10} | {'a':<5} {'b':<5} {'c':<5} {'d':<5} | {'Odds Ratio':<12} | {'P-value'}")
print("-" * 75)

for sheet in SHEET_NAMES:
    try:
        df = pd.read_excel(EXCEL_FILE, sheet_name=sheet)
        df.fillna(0, inplace=True)

        # Hotspot: frequency >= 8 and percentage3 > 0
        is_hotspot = (df['frequency'] >= 8) & (df['percentage3'] > 0)
        # Non-hotspot: 0 < frequency < 8 and percentage3 > 0
        non_hotspot = (df['frequency'] > 0) & (df['frequency'] < 8) & (df['percentage3'] > 0)

        # Binding interface: DDG_MutPNI between 0.5 and 1.0, and percentage3 > 0
        is_binding = (df['DDG_MutPNI'] >= 0.5) & (df['percentage3'] > 0) & (df['DDG_MutPNI'] < 1.0)

        # Contingency table counts
        a = len(df[is_hotspot & is_binding])
        b = len(df[is_hotspot & ~is_binding])
        c = len(df[non_hotspot & is_binding])
        d = len(df[non_hotspot & ~is_binding])

        total_a += a
        total_b += b
        total_c += c
        total_d += d

        odds_ratio, p_value = stats.fisher_exact([[a, b], [c, d]])
        print(f"{sheet:<10} | {a:<5} {b:<5} {c:<5} {d:<5} | {odds_ratio:<12.4f} | {p_value:.4e}")

    except Exception as e:
        print(f"Error processing sheet {sheet}: {e}")

print("-" * 75)

combined_odds_ratio, combined_p_value = stats.fisher_exact([[total_a, total_b], [total_c, total_d]])

print(
    f"{'COMBINED':<10} | {total_a:<5} {total_b:<5} {total_c:<5} {total_d:<5} | {combined_odds_ratio:<12.4f} | {combined_p_value:.4e}")

print("\nResults:")
print(f"a (Hotspot & Binding):     {total_a}")
print(f"b (Hotspot & Non-binding): {total_b}")
print(f"c (Non-hotspot & Binding): {total_c}")
print(f"d (Non-hotspot & Non-binding): {total_d}")
print(f"Odds Ratio: {combined_odds_ratio:.4f}")
print(f"P-value: {combined_p_value:.4e}")

if combined_p_value < 0.05:
    print("Conclusion: P < 0.05, significant association between hotspots and binding interfaces.")
else:
    print("Conclusion: P >= 0.05, no significant association found.")