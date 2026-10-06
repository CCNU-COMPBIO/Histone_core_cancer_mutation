#!/bin/bash
# run_rmsf_block.sh
# Usage: run this script from the level containing WT_rw/, R29Q_rw/, ... (i.e. BASE_DIR):
#     bash run_rmsf_block.sh
#
# Directory layout:
#   BASE_DIR/
#     WT_rw/
#       WT_rw.prmtop              <- topology in the system directory
#       run1/
#         WT_rw_run1_500ns.nc     <- trajectory in the run subdirectory
#         rmsf-*.dat              <- outputs are written here
#       run2/ run3/ ...
#
# Outputs (two sets are generated for each run):
#   Regular (unwrapping segments excluded; identical to _noexcl for runs with no excluded segments):
#     rmsf-<histone>_chain_<chain>.dat
#   Extra (generated for all runs, with no segments excluded):
#     rmsf-<histone>_chain_<chain>_noexcl.dat

set -euo pipefail

# ==================== Configuration ====================
FRAMES_PER_NS=50                 # 500 ns = 25000 frames; sample 1 frame per ns
ANALYSIS_START_NS=100
ANALYSIS_END_NS=500

SYSTEMS=(
  "WT:WT_rw"
  "R29Q:R29Q_rw"
  "E56K:E56K_rw"
  "K74N:K74N_rw"
  "E105K:H2BE105K_rw"
  "F70L:F70L_rw"
  "E50K:E50K_rw"
  "E73K:E73K_rw"
)

declare -A TOPO_MAP=(
  [WT_rw]="WT_rw.prmtop"
  [R29Q_rw]="R29Q_rw.prmtop"
  [E56K_rw]="E56K_rw.prmtop"
  [K74N_rw]="K74N_rw.prmtop"
  [H2BE105K_rw]="E105K_rw.prmtop"
  [F70L_rw]="F70L_rw.prmtop"
  [E50K_rw]="E50K_rw.prmtop"
  [E73K_rw]="E73K_rw.prmtop"
)

# Trajectory filename template: relative to BASE_DIR
TRAJ_TEMPLATE="{run}/{sys}_{run}_500ns.nc"

declare -A EXCLUDE
EXCLUDE[K74N_rw/run1]="253-298,427-469"
EXCLUDE[K74N_rw/run3]="322-379,383-416,422-499"
EXCLUDE[F70L_rw/run3]="262-276"
EXCLUDE[WT_rw/run1]=""
EXCLUDE[WT_rw/run2]=""
EXCLUDE[WT_rw/run3]=""
EXCLUDE[R29Q_rw/run1]=""
EXCLUDE[R29Q_rw/run2]=""
EXCLUDE[R29Q_rw/run3]=""
EXCLUDE[E56K_rw/run1]=""
EXCLUDE[E56K_rw/run2]=""
EXCLUDE[E56K_rw/run3]=""
EXCLUDE[K74N_rw/run2]=""
EXCLUDE[F70L_rw/run1]=""
EXCLUDE[F70L_rw/run2]=""
EXCLUDE[H2BE105K_rw/run1]="455-474,476-499"
EXCLUDE[H2BE105K_rw/run2]=""
EXCLUDE[H2BE105K_rw/run3]=""
EXCLUDE[E50K_rw/run1]=""
EXCLUDE[E50K_rw/run2]=""
EXCLUDE[E50K_rw/run3]=""
EXCLUDE[E73K_rw/run1]=""
EXCLUDE[E73K_rw/run2]=""
EXCLUDE[E73K_rw/run3]=""

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

# ==================== Utility functions ====================
range_to_frames () {
  local s=$1 e=$2
  local out=()
  for (( ns=s; ns<=e; ns++ )); do
    out+=( $(( ns * FRAMES_PER_NS + 1 )) )
  done
  printf '%s\n' "${out[@]}"
}

# Convert a frame-number array into cpptraj trajin lines (grouped by contiguous segments).
# cpptraj trajin accepts only start/stop/offset, not an arbitrary frame list.
emit_cpptraj_frames () {
  local -n _fr=$1
  local fname=$2
  local n=${#_fr[@]}
  if (( n == 0 )); then return; fi
  local seg_start=${_fr[0]}
  local prev=${_fr[0]}
  local cur i
  for (( i=1; i<n; i++ )); do
    cur=${_fr[$i]}
    if (( cur != prev + FRAMES_PER_NS )); then
      if (( seg_start == prev )); then
        echo "trajin ${fname} ${seg_start} ${seg_start} 1"
      else
        echo "trajin ${fname} ${seg_start} ${prev} ${FRAMES_PER_NS}"
      fi
      seg_start=$cur
    fi
    prev=$cur
  done
  if (( seg_start == prev )); then
    echo "trajin ${fname} ${seg_start} ${seg_start} 1"
  else
    echo "trajin ${fname} ${seg_start} ${prev} ${FRAMES_PER_NS}"
  fi
}

# Generate full-segment RMSF files for the given frame array.
# Arguments:
#   $1=sys        system directory
#   $2=run        run directory
#   $3=topo       topology filename
#   $4=traj_file  trajectory filename (relative to the run directory)
#   $5=suffix     output filename suffix ("" or "_noexcl")
#   $6=frames     frame-number array variable name (nameref)
generate_rmsf_set () {
  local sys=$1 run=$2 topo=$3 traj_file=$4 suffix=$5
  local -n _frames=$6
  local n=${#_frames[@]}
  if (( n < 2 )); then
    echo "[WARN] too few frames for ${sys}/${run} (suffix='${suffix}'), skip"
    return
  fi

  local frag range_part filename res_range
  for frag in "${FRAGMENTS[@]}"; do
    range_part="${frag%%:*}"
    filename="${frag##*:}"
    res_range="${range_part%@*}"

    # ---------- Full segment ----------
    {
      echo "parm ../$topo"
      emit_cpptraj_frames _frames "$traj_file"
      echo "autoimage"
      echo "rms first :${res_range}@C,CA,N,O"
      echo "atomicfluct RMSF :${range_part} byres out rmsf-${filename}${suffix}.dat"
      echo "run"
    } > "$sys/$run/.rmsf_tmp_full.in"
    ( cd "$sys/$run" && cpptraj -i .rmsf_tmp_full.in >/dev/null )
    rm -f "$sys/$run/.rmsf_tmp_full.in"

    echo "  [done] ${filename}${suffix}"
  done
}

# ==================== Main loop ====================
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
    traj_rel="${traj_rel//\{sys\}/$sys}"
    traj_file="$(basename "$traj_rel")"
    key="$sys/$run"

    if [[ ! -f "$sys/$traj_rel" ]]; then
      echo "[WARN] missing trajectory, skip: $sys/$traj_rel"
      continue
    fi

    excl_list="${EXCLUDE[$key]:-}"

    # ---------- Build frame set ----------
    mapfile -t all_frames < <(range_to_frames "$ANALYSIS_START_NS" "$ANALYSIS_END_NS")

    declare -A drop=()
    if [[ -n "$excl_list" ]]; then
      IFS=',' read -ra segs <<< "$excl_list"
      for seg in "${segs[@]}"; do
        seg="${seg// /}"
        [[ -z "$seg" ]] && continue
        es="${seg%-*}"; ee="${seg#*-}"
        while IFS= read -r f; do drop[$f]=1; done < <(range_to_frames "$es" "$ee")
      done
    fi

    keep=()
    for f in "${all_frames[@]}"; do
      [[ -n "${drop[$f]:-}" ]] || keep+=( "$f" )
    done

    n_keep=${#keep[@]}
    if (( n_keep < 2 )); then
      echo "[WARN] too few retained frames for $key, skip"
      continue
    fi

    echo "##############################################################"
    echo "#  System: $sys   Run: $run   Topo: $sys/$topo"
    echo "#  Traj  : $sys/$traj_rel"
    echo "#  Excluded segments (ns): ${excl_list:-none}"
    echo "#  All frames     : ${#all_frames[@]}"
    echo "#  Retained frames: ${n_keep}  (dropped $(( ${#all_frames[@]} - n_keep )))"
    echo "##############################################################"

    # ---------- Regular: with exclusions ----------
    generate_rmsf_set "$sys" "$run" "$topo" "$traj_file" "" keep

    # ---------- Extra: no-exclusion version for all runs (no special cases downstream) ----------
    echo "  -- also computing NO-EXCLUSION version (suffix=_noexcl) --"
    generate_rmsf_set "$sys" "$run" "$topo" "$traj_file" "_noexcl" all_frames
  done
done

echo ""
echo "All done."
echo "  Regular (exclusions applied):"
echo "    Full : <sys>/run*/rmsf-<histone>_chain_<chain>.dat"
echo "  Extra (no exclusions; available for all runs):"
echo "    Full : <sys>/run*/rmsf-<histone>_chain_<chain>_noexcl.dat"