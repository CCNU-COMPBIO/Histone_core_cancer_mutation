#!/bin/bash
# run_rmsf_cumulative.sh
# Purpose: time-convergence analysis. For each run, compute cumulative RMSF over
#          0–100, 0–200, 0–300, 0–400, 0–500 ns.
# Differences from run_rmsf_block.sh:
#   1) No segments are excluded; the full trajectory is used (contiguous time
#      windows, easier to compare across runs);
#   2) Output suffixes are _cum100ns / _cum200ns / _cum300ns / _cum400ns / _cum500ns;
#   3) cpptraj is called only once per (system, run, cutoff), all eight chains
#      at once, to reduce I/O.
# Usage: run from the BASE_DIR (the directory containing WT_rw/, H2A_R29Q_rw/, ...):
#        bash run_rmsf_cumulative.sh

set -euo pipefail

FRAMES_PER_NS=50
START_NS=0                 # start from the first frame of the production trajectory
CUTOFFS=(100 200 300 400 500)

SYSTEMS=(
  "WT:WT_rw"
  "R29Q:H2A_R29Q_rw"
  "E56K:H2A_E56K_rw"
  "K74N:H2A_K74N_rw"
  "H2BE105K:H2B_E105K_rw"
  "F70L:H2B_F70L_rw"
  "E50K:H3_E50K_rw"
  "E73K:H3_E73K_rw"
  "E105K:H3_E105K_rw"
)

declare -A TOPO_MAP=(
  [WT_rw]="WT_rw.prmtop"
  [H2A_R29Q_rw]="R29Q_rw.prmtop"
  [H2A_E56K_rw]="E56K_rw.prmtop"
  [H2A_K74N_rw]="K74N_rw.prmtop"
  [H2B_E105K_rw]="E105K_rw.prmtop"
  [H2B_F70L_rw]="F70L_rw.prmtop"
  [H3_E50K_rw]="E50K_rw.prmtop"
  [H3_E73K_rw]="E73K_rw.prmtop"
  [H3_E105K_rw]="E105K_rw.prmtop"
)

TRAJ_TEMPLATE="{run}/{traj_label}_rw_{run}_500ns.nc"

# Chain definitions and residue ranges, matching the main analysis
FRAGMENTS=(
  "1-97@CA:H3_chain_A"
  "98-175@CA:H4_chain_B"
  "176-280@CA:H2A_chain_C"
  "281-376@CA:H2B_chain_D"
  "377-475@CA:H3_chain_E"
  "476-557@CA:H4_chain_F"
  "558-661@CA:H2A_chain_G"
  "662-753@CA:H2B_chain_H"
)

RUNS=(run1 run2 run3)

first_frame=$(( START_NS * FRAMES_PER_NS + 1 ))

for entry in "${SYSTEMS[@]}"; do
  label="${entry%%:*}"
  sys="${entry##*:}"
  topo="${TOPO_MAP[$sys]}"

  if [[ ! -f "$sys/$topo" ]]; then
    echo "[WARN] missing topology, skip system: $sys/$topo"
    continue
  fi

  for run in "${RUNS[@]}"; do
    traj_rel="${TRAJ_TEMPLATE//\{run\}/$run}"
    traj_rel="${traj_rel//\{traj_label\}/$label}"
    traj_file="$(basename "$traj_rel")"

    if [[ ! -f "$sys/$traj_rel" ]]; then
      echo "[WARN] missing trajectory, skip: $sys/$traj_rel"
      continue
    fi

    for T in "${CUTOFFS[@]}"; do
      last_frame=$(( T * FRAMES_PER_NS + 1 ))
      inp="$sys/$run/.rmsf_cum${T}_tmp.in"

      {
        echo "parm ../$topo"
        echo "trajin $traj_file $first_frame $last_frame $FRAMES_PER_NS"
        echo "autoimage"
        for frag in "${FRAGMENTS[@]}"; do
          range_part="${frag%%:*}"      # e.g. 176-280@CA
          filename="${frag##*:}"        # e.g. H2A_chain_C
          res_range="${range_part%@*}"  # e.g. 176-280
          echo "rms first :${res_range}@C,CA,N,O"
          echo "atomicfluct rmsf_${filename}_cum${T} :${range_part} byres out rmsf-${filename}_cum${T}ns.dat"
        done
        echo "run"
        echo "quit"
      } > "$inp"

      ( cd "$sys/$run" && cpptraj -i ".rmsf_cum${T}_tmp.in" >/dev/null )
      rm -f "$sys/$run/.rmsf_cum${T}_tmp.in"

      echo "[done] $sys/$run  cum${T}ns"
    done
  done
done

echo ""
echo "All done. Output: <sys>/run*/rmsf-<histone>_chain_<chain>_cum{100,200,300,400,500}ns.dat"