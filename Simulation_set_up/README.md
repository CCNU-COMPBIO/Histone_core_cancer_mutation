# Simulation Set Up

This directory contains initial structures and simulation setup templates for nucleosome molecular dynamics simulations.

## Contents

| File / Folder             | Description                                                                 |
| :------------------------ | :-------------------------------------------------------------------------- |
| [NUC_models](./NUC_models/) | Directory storing initial structure files for all nucleosome (NUC) systems  |
| [WT_simulation_templates](./WT_simulation_templates/)| Wild-type (WT) reference templates for simulation setup, job submission, and analysis preprocessing. Copy and modify these files when setting up mutant or variant systems. |

> **Note**: All non-structure setup files (AMBER input scripts, tleap/VMD scripts) are organized inside `WT_simulation_templates/`. The `NUC_models/` directory exclusively stores initial PDB/coordinate files for all simulated systems.
