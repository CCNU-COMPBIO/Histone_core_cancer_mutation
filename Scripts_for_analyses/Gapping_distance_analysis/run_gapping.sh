#!/bin/bash
# =============================================================
# DNA Gapping Distance Analysis using cpptraj
# System: no water/ions (already stripped from topology/trajectory)
# 500 ns / 25000 frames / 20 ps per frame
# Analysis range: 100 ns – 500 ns, sampled every 1 ns (50 frames)
# SHL sites based on corrected shl_ranges (only adjust 3rd/4th columns ±1, keep pairing)
# =============================================================

mkdir -p distance_results

RUNS=("run1" "run2" "run3")

FIRST=5001   # start at 100 ns (frame 5001)
LAST=25000   # end at 500 ns (frame 25000)
STEP=50      # every 1 ns

for RUN in "${RUNS[@]}"; do
    echo "=========================================="
    echo "Processing: ${RUN}"
    echo "=========================================="

    cat > cpptraj_${RUN}.in << EOF
parm WT_rw.prmtop
trajin ${RUN}/WT_rw_${RUN}_500ns.nc ${FIRST} ${LAST} ${STEP}

# Periodic boundary handling (comment out if trajectory already imaged)
autoimage

# RMSD alignment: reference to histone core CA (residues 1-753)
rms first :1-753@CA

# ----- 7 SHL gapping distances (final paired version) -----

# SHL -1 / +7
distance gap1 (:817-819|:980-982) (:894-896|:903-905) type center out distance_results/WT_${RUN}_SHL_1_7_gap.dat

# SHL -2 / +6
distance gap2 (:807-809|:990-992) (:884-886|:913-915) type center out distance_results/WT_${RUN}_SHL_2_6_gap.dat

# SHL -3 / +5
distance gap3 (:797-799|:1000-1002) (:874-876|:923-925) type center out distance_results/WT_${RUN}_SHL_3_5_gap.dat

# SHL -4 / +4
distance gap4 (:787-789|:1010-1012) (:864-866|:933-935) type center out distance_results/WT_${RUN}_SHL_4_4_gap.dat

# SHL -5 / +3
distance gap5 (:777-779|:1020-1022) (:854-856|:943-945) type center out distance_results/WT_${RUN}_SHL_5_3_gap.dat

# SHL -6 / +2
distance gap6 (:767-769|:1030-1032) (:844-846|:953-955) type center out distance_results/WT_${RUN}_SHL_6_2_gap.dat

# SHL -7 / +1
distance gap7 (:757-759|:1040-1042) (:834-836|:961-963) type center out distance_results/WT_${RUN}_SHL_7_1_gap.dat

run
quit
EOF

    cpptraj -i cpptraj_${RUN}.in | tee cpptraj_${RUN}.log
    echo "${RUN} done."
done

echo "All runs completed."