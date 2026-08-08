# Source Data

This directory contains raw and processed data files generated from nucleosome molecular dynamics simulations.

## Subdirectories

-   [Energetically_important_sites_data](./Energetically_important_sites_data/) – The raw and processed data for the alanine scanning to identify the energetically important sites
-   [Gapping_distance_analysis](./Gapping_distance_analysis/) – Nucleosomal DNA gapping distance measurement data
-   [Histone_mutation_data](./Histone_mutation_data/) – Curated histone mutation datasets for core regions
-   [Jaccard_similarity_analysis](./Jaccard_similarity_analysis/) – Jaccard similarity matrices and sub-cluster overlap data
-   [MMGBSA_calculation](./MMGBSA_calculation/) – Raw and post-processed MM/GBSA binding free energy calculation outputs
-   [RMSD_analysis](./RMSD_analysis/) – RMSD trajectory time-series data
-   [RMSF_analysis](./RMSF_analysis/) – Raw and computed RMSF-related data
-   [SASA_calculation](./SASA_calculation/) – Per-residue SASA (Solvent-accessible surface area) calculated values for each histone
-   [hbond_analysis](./hbond_analysis/) – Hydrogen bond calculations and statistical data

## Notes

-   All data files in this directory are intended to be consumed by the corresponding analysis scripts in [`../Scripts_for_Analyses/`](../Scripts_for_Analyses/).
