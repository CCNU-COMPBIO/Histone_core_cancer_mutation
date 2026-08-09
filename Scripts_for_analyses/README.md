# Scripts for Analyses

This directory contains modular analysis scripts organized by computational task. Each subdirectory holds scripts and auxiliary files for a specific type of analysis.

## Subdirectories

-   [Cluster_analysis](./Cluster_analysis/) – Clustering analysis targeting energetically important residues and mutation hotspots
-   [Contact_number_analysis](./Contact_number_analysis/) – Scripts for analyzing heavy atom contact numbers between two chains
-   [Energetically_important_sites_analysis](./Energetically_important_sites_analysis/) – Identification and distribution analysis of energetically important residues
-   [Histone_mutation_analysis](./Histone_mutation_analysis/) – Analysis of histone mutations in core regions, including mutation frequency and charge changes
-   [Jaccard_similarity_analysis](./Jaccard_similarity_analysis/) – Jaccard similarity analysis for residue sets and sub-clusters
-   [MMGBSA_calculation](./MMGBSA_calculation/) – Scripts for running and post-processing MM/GBSA binding free energy calculations
-   [RMSD_analysis](./RMSD_analysis/) – Root-mean-square deviation (RMSD) trajectory analysis and plotting
-   [RMSF_analysis](./RMSF_analysis/) – Root-mean-square fluctuation (RMSF) per-residue flexibility analysis
-   [SASA_calculation](./SASA_calculation/) – Solvent-accessible surface area (SASA) calculation
-   [hbond_analysis](./hbond_analysis/) – Hydrogen bond analysis

## Notes

-   Input data are expected to reside in [`../Source_data/`](../Source_data/).
