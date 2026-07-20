#!/bin/bash

# Define fragment ranges and output filenames
# Format: residue_range@atom_selection:output_filename
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

for frag in "${fragments[@]}"; do
    # Split: extract the range part (e.g., "1-97@CA") and the output filename
    range_part=$(echo $frag | cut -d':' -f1)
    filename=$(echo $frag | cut -d':' -f2)
    
    # Extract the pure residue range (remove the "@CA" suffix) for alignment mask
    res_range=$(echo $range_part | sed 's/@.*//')   # yields "1-97", etc.
    
    cpptraj << EOF
parm ../WT_rw.prmtop
trajin WT_rw_run1_500ns.nc
autoimage
# Align to the first frame using backbone atoms (C,CA,N,O) of the current chain
rms first :${res_range}@C,CA,N,O
# Compute RMSF for the specified residue range and atom selection (e.g., CA)
atomicfluct RMSF :${range_part} byres out rmsf-${filename}.dat
run
EOF
done