#!/bin/bash
set -euo pipefail

TOPO="WT_rw.prmtop"
FRAMES_PER_NS=50
OFFSET=50
START_FRAME=$((100 * FRAMES_PER_NS + 1))   # 5001
END_FRAME=$((500 * FRAMES_PER_NS))         # 25000

OUT_ROOT="unwrapping_out"
mkdir -p "${OUT_ROOT}"

declare -A TRAJ_MAP
TRAJ_MAP[run1]="run1/WT_rw_run1_500ns.nc"
TRAJ_MAP[run2]="run2/WT_rw_run2_500ns.nc"
TRAJ_MAP[run3]="run3/WT_rw_run3_500ns.nc"

# DNA end definitions
END1_I="754-763"
END1_J="1036-1045"
END2_I="890-899"
END2_J="900-909"
DNA_ATOMS="@P,O5',C5',C4',C3'"

for run in run1 run2 run3; do
    TRAJ="${TRAJ_MAP[$run]}"
    [[ -f "$TRAJ" ]] || { echo "[WARN] skip $TRAJ"; continue; }

    OUT_DIR="${OUT_ROOT}/${run}"
    mkdir -p "${OUT_DIR}"

    echo "### ${run}: ${TRAJ}"

    cpptraj <<EOF
parm ${TOPO}
trajin ${TRAJ} ${START_FRAME} ${END_FRAME} ${OFFSET}
autoimage
rms first :1-753@C,CA,N,O
rms end1_rmsd first :${END1_I},${END1_J}${DNA_ATOMS} nofit out ${OUT_DIR}/end1_rmsd.dat
rms end2_rmsd first :${END2_I},${END2_J}${DNA_ATOMS} nofit out ${OUT_DIR}/end2_rmsd.dat
run
EOF
done

echo "Done. Output: ${OUT_ROOT}/runX/end{1,2}_rmsd.dat"