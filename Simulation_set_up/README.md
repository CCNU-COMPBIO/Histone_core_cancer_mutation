# Simulation Set Up

This directory contains input files, scripts, and initial structures for nucleosome molecular dynamics simulations.

## Contents

| File / Folder             | Description                                                                 |
| :------------------------ | :-------------------------------------------------------------------------- |
| `NUC_models/`             | Directory storing initial structure files for all nucleosome (NUC) systems  |
| `Equil_pt.in`             | Input file for NPT equilibration simulation                                 |
| `Equil_v.in`              | Input file for NVT equilibration simulation                                 |
| `Min.in`                  | Input file for system energy minimization                                   |
| `Prod.in`                 | Input file for production MD simulation                                     |
| `sub_job_GPU.sh`          | Slurm job submission script for GPU-accelerated simulations                 |
| `wt_gen_nucl.leap`        | tleap script for WT nucleosome modeling, topology & coordinate generation   |
| `wt_mmpbsa_topology.leap` | tleap script to generate topology specifically for MMPBSA free energy calculation |
| `wt_sep_bp.tcl`           | VMD TCL script for base pair splitting and nucleosome structure preprocessing |

> **Note**: Files prefixed with `wt_` are wild-type specific setup scripts. The `NUC_models/` directory contains initial structures for all simulated systems.
