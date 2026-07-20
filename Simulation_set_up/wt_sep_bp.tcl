mol load pdb wt.pdb

set bp [atomselect top "protein and (not chain A B C D E F G H) and noh"]
set nuc [atomselect top "(nucleic or protein) and (chain A B C D E F G H I J) and noh"]
set com [atomselect top "(protein or nucleic) and noh"]

$bp writepdb wt_bp.pdb
$nuc writepdb wt_nuc.pdb
$com writepdb wt_com.pdb
exit
