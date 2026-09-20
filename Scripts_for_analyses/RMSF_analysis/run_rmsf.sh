#!/bin/bash
set -euo pipefail

# ========== Configuration ==========
# Run this in the directory containing run1/ run2/ run3/ and the prmtop file
TOPO="WT_rw.prmtop"

# Trajectory frame rate: 500 ns total, 25000 frames -> 50 frames/ns
FRAMES_PER_NS=50

OFFSET=50

declare -a runs=(run1 run2 run3)

# Trajectories live inside their respective run subdirectories
TRAJ_TEMPLATE="{run}/WT_rw_{run}_500ns.nc"

# Fragment definitions: residue range@atom:output prefix
declare -a fragments=(
    "1-97@CA:H3_chain_A"
    "98-175@CA:H4_chain_B"
    "176-280@CA:H2A_chain_C"
    "281-376@CA:H2B_chain_D"
    "377-475@CA:H3_chain_E"
    "476-557@CA:H4_chain_F"
    "558-661@CA:H2A_chain_G"
    "662-753@CA:H2B_chain_H"
)

# ========== Main loop ==========
for run in "${runs[@]}"; do
    TRAJ="${TRAJ_TEMPLATE//\{run\}/$run}"

    if [[ ! -f "$TRAJ" ]]; then
        echo "[WARN] Missing trajectory, skip: $TRAJ"
        continue
    fi

    echo "##############################################################"
    echo "#  Run: $run   Traj: $TRAJ   OFFSET=${OFFSET}"
    echo "##############################################################"

    for frag in "${fragments[@]}"; do
        range_part="${frag%%:*}"          # e.g. 1-97@CA
        filename="${frag##*:}"            # e.g. H3_chain_A
        res_range="${range_part%@*}"      # e.g. 1-97

        # ---------- (A) Two halves: 100-300 ns and 300-500 ns ----------
        declare -a half_windows=(
            "100 300"
            "300 500"
        )
        for win in "${half_windows[@]}"; do
            read -r start_ns end_ns <<< "$win"
            start_frame=$(( start_ns * FRAMES_PER_NS + 1 ))
            end_frame=$(( end_ns   * FRAMES_PER_NS     ))
            n_frames=$(( (end_frame - start_frame) / OFFSET + 1 ))

            out="${run}/rmsf-${filename}_${start_ns}to${end_ns}ns.dat"

            echo "  [half] ${filename}  ${start_ns}-${end_ns} ns  (frames ${start_frame}-${end_frame}, step ${OFFSET}, ~${n_frames} frames)  -> ${out}"
            cpptraj <<EOF
parm ${TOPO}
trajin ${TRAJ} ${start_frame} ${end_frame} ${OFFSET}
autoimage
rms first :${res_range}@C,CA,N,O
atomicfluct RMSF :${range_part} byres out ${out}
run
EOF
        done

        # ---------- (B) Full window: 100-500 ns, filename without time suffix ----------
        start_frame=$(( 100 * FRAMES_PER_NS + 1 ))   # 5001
        end_frame=$(( 500   * FRAMES_PER_NS     ))   # 25000
        n_frames=$(( (end_frame - start_frame) / OFFSET + 1 ))

        out="${run}/rmsf-${filename}.dat"

        echo "  [full] ${filename}  100-500 ns  (frames ${start_frame}-${end_frame}, step ${OFFSET}, ~${n_frames} frames)  -> ${out}"
        cpptraj <<EOF
parm ${TOPO}
trajin ${TRAJ} ${start_frame} ${end_frame} ${OFFSET}
autoimage
rms first :${res_range}@C,CA,N,O
atomicfluct RMSF :${range_part} byres out ${out}
run
EOF
    done
done

echo ""
echo "All done.  OFFSET=${OFFSET}"
echo "  Full-window files : run*/rmsf-<histone>_chain_<chain>.dat"
echo "  Half-window files : run*/rmsf-<histone>_chain_<chain>_100to300ns.dat"
echo "                      run*/rmsf-<histone>_chain_<chain>_300to500ns.dat"
