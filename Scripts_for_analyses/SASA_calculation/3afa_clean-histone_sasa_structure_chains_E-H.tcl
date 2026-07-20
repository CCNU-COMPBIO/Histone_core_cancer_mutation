mol load pdb ./3afa_clean.pdb

set outfile1 [open ./3afa_clean-H3_sasa_nuc_structure.dat w]
set outfile2 [open ./3afa_clean-H4_sasa_nuc_structure.dat w]
set outfile3 [open ./3afa_clean-H2A_sasa_nuc_structure.dat w]
set outfile4 [open ./3afa_clean-H2B_sasa_nuc_structure.dat w]

#sasa for H3
for { set r1 1 } { $r1<=135 } { incr r1 } {

set sel1 [atomselect top "(chain E and resid $r1) and noh and protein"]
set sel2 [atomselect top "(chain E) and noh and protein"]
set sel3 [atomselect top "(chain A B C D E F G H) and noh and protein"]
set sel4 [atomselect top "(chain A B C D E F G H I J) and noh and (protein or nucleic)"]

set sasa1 [measure sasa 1.4  $sel1]
set sasa2 [measure sasa 1.4  $sel2 -restrict $sel1]
set sasa3 [measure sasa 1.4  $sel3 -restrict $sel1]
set sasa4 [measure sasa 1.4  $sel4 -restrict $sel1]

puts -nonewline $outfile1 "$r1\t$sasa1\t$sasa2\t$sasa3\t$sasa4\n"

unset sel1 
unset sel2 
unset sel3 
unset sel4
unset sasa1 
unset sasa2 
unset sasa3
unset sasa4
#unset ratio 

}

#sasa for H4
for { set r2 1 } { $r2<=102 } { incr r2 } {

set sel1 [atomselect top "(chain F and resid $r2) and noh and protein"]
set sel2 [atomselect top "(chain F) and noh and protein"]
set sel3 [atomselect top "(chain A B C D E F G H) and noh and protein"]
set sel4 [atomselect top "(chain A B C D E F G H I J) and noh and (protein or nucleic)"]

set sasa1 [measure sasa 1.4  $sel1]
set sasa2 [measure sasa 1.4  $sel2 -restrict $sel1]
set sasa3 [measure sasa 1.4  $sel3 -restrict $sel1]
set sasa4 [measure sasa 1.4  $sel4 -restrict $sel1]

puts -nonewline $outfile2 "$r2\t$sasa1\t$sasa2\t$sasa3\t$sasa4\n"

unset sel1 
unset sel2 
unset sel3 
unset sel4
unset sasa1 
unset sasa2 
unset sasa3
unset sasa4
#unset ratio 
}

#sasa for H2A
for { set r3 1 } { $r3<=130 } { incr r3 } {

set sel1 [atomselect top "(chain G and resid $r3) and noh and protein"]
set sel2 [atomselect top "(chain G) and noh and protein"]
set sel3 [atomselect top "(chain A B C D E F G H) and noh and protein"]
set sel4 [atomselect top "(chain A B C D E F G H I J) and noh and (protein or nucleic)"]

set sasa1 [measure sasa 1.4  $sel1]
set sasa2 [measure sasa 1.4  $sel2 -restrict $sel1]
set sasa3 [measure sasa 1.4  $sel3 -restrict $sel1]
set sasa4 [measure sasa 1.4  $sel4 -restrict $sel1]

puts -nonewline $outfile3 "$r3\t$sasa1\t$sasa2\t$sasa3\t$sasa4\n"

unset sel1 
unset sel2 
unset sel3 
unset sel4
unset sasa1 
unset sasa2 
unset sasa3
unset sasa4
#unset ratio 
}

#sasa for H2B
for { set r4 1 } { $r4<=125 } { incr r4 } {

set sel1 [atomselect top "(chain H and resid $r4) and noh and protein"]
set sel2 [atomselect top "(chain H) and noh and protein"]
set sel3 [atomselect top "(chain A B C D E F G H) and noh and protein"]
set sel4 [atomselect top "(chain A B C D E F G H I J) and noh and (protein or nucleic)"]

set sasa1 [measure sasa 1.4  $sel1]
set sasa2 [measure sasa 1.4  $sel2 -restrict $sel1]
set sasa3 [measure sasa 1.4  $sel3 -restrict $sel1]
set sasa4 [measure sasa 1.4  $sel4 -restrict $sel1]

puts -nonewline $outfile4 "$r4\t$sasa1\t$sasa2\t$sasa3\t$sasa4\n"

unset sel1 
unset sel2 
unset sel3 
unset sel4
unset sasa1 
unset sasa2 
unset sasa3
unset sasa4
#unset ratio 
}

close $outfile1 
close $outfile2
close $outfile3
close $outfile4

exit


