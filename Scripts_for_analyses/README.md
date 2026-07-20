# Scripts for Analyses

This directory contains modular analysis scripts organized by computational task. Each subdirectory holds scripts and auxiliary files for a specific type of analysis.

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `Cluster_analysis` | Clustering of MD conformations (e.g., by RMSD or dihedral space) |
| `Energetically_important_sites_analysis` | Identification of key residues contributing to binding energy (e.g., via per-residue decomposition) |
| `Histone_mutation_analysis` | Analysis of histone variant/mutation effects on dynamics or interactions |
| `Jaccard_similarity_analysis` | Quantification of residue contact similarity between systems |
| `MMGBSA_calculation` | Scripts for running and post-processing MM/GBSA binding free energy calculations |
| `RMSD_analysis` | Root-mean-square deviation (RMSD) trajectory analysis and plotting |
| `RMSF_analysis` | Root-mean-square fluctuation (RMSF) per-residue flexibility analysis |
| `SASA_calculation` | Solvent-accessible surface area (SASA) calculation and interpretation |
| `hbond_analysis` | Hydrogen bond detection, lifetime, and occupancy analysis |

## Notes

- Input data (trajectories, topologies, etc.) are expected to reside in [`../Source_data/`](../Source_data/).
