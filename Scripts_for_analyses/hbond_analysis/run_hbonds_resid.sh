#!/bin/bash
# ============================================================
# Script: compute_residue_hbonds.sh
# Purpose: Compute hydrogen bonds between a specified residue
#          (e.g. residue 70 of chain D) and the rest of the system
#          (other protein residues + DNA).
# Directory layout:
#   Place this script in the same directory as WT_rw.prmtop.
#   Trajectories in ./run1/WT_rw_run1_500ns.nc, ./run2/..., ./run3/...
#   Results are written to the corresponding run subdirectories.
# Dependencies: VMD + hbonds.tcl
# ============================================================
# Run from the script directory:
# ./compute_residue_hbonds.sh -s D:70 -r 1
# ./compute_residue_hbonds.sh -s D:70 -r all

# ========== Default parameters ==========
PRMTOP="WT_rw.prmtop"                                    # Topology file (same level as script)
TRAJ_TEMPLATE="./run{RUN}/WT_rw_run{RUN}_500ns.nc"      # Trajectory template; {RUN} is replaced
FIRST=5000                                               # First frame
LAST=25000                                               # Last frame
STEP=50                                                  # Frame stride
DIST_CUTOFF=3.5                                          # H-bond distance cutoff (Å)
ANGLE_CUTOFF=30                                          # H-bond angle cutoff (degrees)
VMD_PLUGIN="/home/lenovo/vmd-1.9.4a57/plugins/noarch/tcl/hbonds1.3/hbonds.tcl"
OUTDIR_BASE="./"                                         # Base output directory

# ========== Global residue ranges for chains (VMD resid numbering) ==========
declare -A CHAIN_RANGE=(
  [A]="1 to 97"
  [B]="98 to 175"
  [C]="176 to 280"
  [D]="281 to 376"
  [E]="377 to 475"
  [F]="476 to 557"
  [G]="558 to 661"
  [H]="662 to 753"
  [I]="754 to 899"   # DNA chain I
  [J]="900 to 1045"  # DNA chain J
)

# ========== Local starting residue number for each chain ==========
declare -A CHAIN_LOCAL_START=(
  [A]=38
  [B]=25
  [C]=14
  [D]=30
  [E]=37
  [F]=21
  [G]=15
  [H]=33
)

# ========== Help ==========
usage() {
  cat <<EOF
Usage: $0 -s CHAIN:LOCAL_RESID -r RUN [options]

Required arguments:
  -s CHAIN:LOCAL_RESID   Chain and local residue number, e.g. -s D:70
  -r RUN                 Run number (1,2,3 or all)

Optional arguments:
  -i PRMTOP              Topology file path (default: $PRMTOP in current directory)
  -t TEMPLATE            Trajectory template; {RUN} is replaced (default: $TRAJ_TEMPLATE)
  --first N              First frame (default: $FIRST)
  --last N               Last frame (default: $LAST)
  --step N               Frame stride (default: $STEP)
  --dist D               H-bond distance cutoff (default: $DIST_CUTOFF)
  --ang A                H-bond angle cutoff (default: $ANGLE_CUTOFF)
  --vmd-plugin P         Path to hbonds.tcl (default: $VMD_PLUGIN)
  --outdir DIR           Base output directory; output goes to DIR/run{RUN}/ (default: current directory)
  -h                     Show this help

Examples:
  $0 -s D:70 -r 1
  $0 -s D:70 -r all
  $0 -s D:70 -r 1 -i ../WT_rw.prmtop -t ../run{RUN}/WT_rw_run{RUN}_500ns.nc
EOF
  exit 0
}

# ========== Parse command-line arguments ==========
single_spec=""
runs=()
all_runs=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -s) single_spec="$2"; shift 2 ;;
    -r)
      if [ "$2" = "all" ]; then
        all_runs=true
      else
        runs+=("$2")
      fi
      shift 2 ;;
    -i) PRMTOP="$2"; shift 2 ;;
    -t) TRAJ_TEMPLATE="$2"; shift 2 ;;
    --first) FIRST="$2"; shift 2 ;;
    --last) LAST="$2"; shift 2 ;;
    --step) STEP="$2"; shift 2 ;;
    --dist) DIST_CUTOFF="$2"; shift 2 ;;
    --ang) ANGLE_CUTOFF="$2"; shift 2 ;;
    --vmd-plugin) VMD_PLUGIN="$2"; shift 2 ;;
    --outdir) OUTDIR_BASE="$2"; shift 2 ;;
    -h) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# Check required arguments
if [ -z "$single_spec" ]; then
  echo "Error: -s must be specified, e.g. -s D:70"
  usage
fi

# Determine run list
if $all_runs; then
  runs=(1 2 3)
fi
if [ ${#runs[@]} -eq 0 ]; then
  echo "Error: -r must specify one or more run numbers, or use -r all"
  exit 1
fi

# Check VMD and plugin
if ! command -v vmd &> /dev/null; then
  echo "Error: vmd command not found in PATH"
  exit 1
fi
if [ ! -f "$VMD_PLUGIN" ]; then
  echo "Error: hbonds plugin not found: $VMD_PLUGIN"
  exit 1
fi

# Check topology file
if [ ! -f "$PRMTOP" ]; then
  echo "Error: topology file '$PRMTOP' does not exist"
  exit 1
fi

TOPOBASE=$(basename "$PRMTOP" .prmtop)

# ========== Helper: local residue number -> global residue number ==========
local_to_global() {
  local chain="$1"
  local local_res="$2"
  local range="${CHAIN_RANGE[$chain]}"
  local local_start="${CHAIN_LOCAL_START[$chain]}"
  if [ -z "$range" ]; then
    echo "Error: unknown chain identifier '$chain'" >&2
    exit 1
  fi
  if [ -z "$local_start" ]; then
    echo "Error: chain '$chain' has no local start number defined; please update CHAIN_LOCAL_START" >&2
    exit 1
  fi
  # Extract global start number (first number in the range string)
  local global_start=$(echo "$range" | awk '{print $1}')
  echo $(( global_start + local_res - local_start ))
}

# ========== Main: single-residue H-bond analysis ==========
overall_success=true

# Parse -s argument, e.g. "D:70"
if [[ "$single_spec" != *:* ]]; then
  echo "Error: -s format must be CHAIN:LOCAL_RESID, e.g. D:70"
  exit 1
fi
chain="${single_spec%%:*}"
local_resid="${single_spec##*:}"
global_resid=$(local_to_global "$chain" "$local_resid")
echo ">>> Single-residue mode: chain $chain, local residue $local_resid -> global residue $global_resid"

# DNA selection (chains I and J from CHAIN_RANGE)
dna_sel="residue ${CHAIN_RANGE[I]} or residue ${CHAIN_RANGE[J]}"
# sel2: all protein residues except the target residue, plus DNA
sel2_expr="(not (same residue as within 0 of resid $global_resid)) or ($dna_sel)"

for run in "${runs[@]}"; do
  # Build trajectory file path
  traj_file="${TRAJ_TEMPLATE//\{RUN\}/$run}"
  if [ ! -f "$traj_file" ]; then
    echo "Warning: trajectory file '$traj_file' does not exist, skipping run $run"
    continue
  fi

  # Output directory: $OUTDIR_BASE/run{RUN}/
  outdir_run="${OUTDIR_BASE}/run${run}"
  mkdir -p "$outdir_run"

  # Output file names use topology name and residue label
  outbase="${TOPOBASE}_run${run}_resid_${chain}${local_resid}"
  detail="${outdir_run}/${outbase}-detail.dat"
  log="${outdir_run}/${outbase}.log"
  outfile="${outdir_run}/${outbase}-hbonds.dat"

  echo ">>> Processing run $run, residue $single_spec -> output: $outfile"

  vmd_script=$(mktemp /tmp/vmd_hbonds_XXXXXX.tcl)
  cat > "$vmd_script" <<EOF
mol delete all
mol load parm7 $PRMTOP
mol addfile $traj_file first $FIRST last $LAST step $STEP waitfor all

set sel1 [atomselect top "resid $global_resid"]
set sel2 [atomselect top "$sel2_expr"]

source $VMD_PLUGIN
hbonds -sel1 \$sel1 -sel2 \$sel2 -dist $DIST_CUTOFF -ang $ANGLE_CUTOFF -plot no -type pair -writefile yes -detailout $detail -log $log -outfile $outfile

exit
EOF

  vmd -dispdev text -e "$vmd_script"
  if [ $? -ne 0 ]; then
    echo "Error: VMD run failed, run $run, residue $single_spec"
    overall_success=false
  fi
  rm -f "$vmd_script"
done

if $overall_success; then
  echo "All single-residue H-bond analyses completed."
else
  echo "Some analyses failed. Please check the error messages above."
  exit 1
fi