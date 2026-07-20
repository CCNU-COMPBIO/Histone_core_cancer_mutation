# Source Data

This directory contains molecular dynamics trajectories, per-frame MM/GBSA calculation results, and summary statistics used in this study.

## File List

| Filename | Type | Description |
| :--- | :--- | :--- |
| `md_0_1.xtc` | Trajectory | MD simulation trajectory file (GROMACS format) |
| `mmpbsa_perframe.csv` | Analysis | Per-frame MM/GBSA binding free energy calculations |
| `summary_mmpbsa.xlsx` | Summary | Mean and standard deviation of MM/GBSA results across systems |

## Usage Notes

- **Trajectory files**: Must be loaded with the corresponding topology file (`.tpr` / `.gro`). Recommended viewers: VMD or PyMOL.
- **CSV data**: Can be read directly with pandas / R. Column definitions are documented in the analysis scripts under `../Scripts_for_analyses/`.
- **Excel summary**: Contains multiple sheets, each corresponding to a different protein–ligand system.

## Related Resources

- Analysis scripts: [Scripts_for_analyses](../Scripts_for_analyses/)
- Simulation parameters: [Simulations_set_up](../Simulations_set_up/)
