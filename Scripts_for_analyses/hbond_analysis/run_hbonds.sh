#!/bin/bash
# 功能: 批量运行 VMD 氢键分析，支持多 run、任意链组间分析（蛋白-蛋白、蛋白-DNA等）
# 目录结构要求:
#   脚本与拓扑文件（默认 WT_rw.prmtop）放在同一目录下
#   轨迹文件位于 ./run1/WT_rw_run1_500ns.nc, ./run2/..., ./run3/...
#   结果文件也输出到对应的 run 子目录中
# 用法示例:
#   ./run_hbonds.sh -p CD -r 1 -r 2            # 链C与链D
#   ./run_hbonds.sh -p A:IJ -r all              # 链A与DNA双链IJ
#   ./run_hbonds.sh -p AB:IJ -r 1 -r 2          # 链A+B与DNA双链IJ
#   ./run_hbonds.sh -p A:B -p C:IJ -r all       # 多组分析
#   ./run_hbonds.sh -i ../F70L_rw.prmtop -t ../run{RUN}/F70L_rw_run{RUN}_500ns.nc -p A:IJ -r 1

# ========== 默认参数 ==========
PRMTOP="WT_rw.prmtop"                                    # 拓扑文件（与脚本同级）
TRAJ_TEMPLATE="./run{RUN}/WT_rw_run{RUN}_500ns.nc"      # 轨迹模板，{RUN} 会被替换
FIRST=5000
LAST=25000
STEP=50
DIST_CUTOFF=3.5
ANGLE_CUTOFF=30
VMD_PLUGIN="/home/lenovo/vmd-1.9.4a57/plugins/noarch/tcl/hbonds1.3/hbonds.tcl"
OUTDIR_BASE="./"                                         # 输出基本目录（将在每个 run 子目录中生成文件）

# ========== 链残基范围 (VMD atomselect 语法) ==========
declare -A CHAIN_RANGE=(
  [A]="1 to 97"
  [B]="98 to 175"
  [C]="176 to 280"
  [D]="281 to 376"
  [E]="377 to 475"
  [F]="476 to 557"
  [G]="558 to 661"
  [H]="662 to 753"
  [I]="754 to 899"
  [J]="900 to 1045"
)

usage() {
  cat <<EOF
Usage: $0 [options]
Options:
  -p SPEC       Chain selection pair. Two formats:
                  - "AB" (shorthand for A:B)  —  two single chains
                  - "A:IJ"                    —  group A vs group I+J
                Repeatable. Group letters allowed: A-J.
  -r RUN        Run number (1, 2, 3). Repeatable, or '-r all' for all 3 runs.
  -i PRMTOP     Topology file (default: $PRMTOP)
  -t TEMPLATE   Trajectory template with {RUN} placeholder (default: $TRAJ_TEMPLATE)
  --first N     First frame (default: $FIRST)
  --last N      Last frame (default: $LAST)
  --step N      Frame step (default: $STEP)
  --dist D      H-bond distance cutoff (default: $DIST_CUTOFF)
  --ang A       H-bond angle cutoff (default: $ANGLE_CUTOFF)
  --vmd-plugin P Path to hbonds.tcl (default: $VMD_PLUGIN)
  --outdir DIR  Output base directory (default: $OUTDIR_BASE)
  -h            Show this help.
Examples:
  $0 -p CD -r 1 -r 2               # Chain C vs D
  $0 -p A:IJ -r all                # Chain A vs DNA (I+J)
  $0 -p AB:IJ -r 1                 # Chains A+B vs DNA
  $0 -p A:B -p C:IJ -r all         # Multiple analyses in one run
EOF
  exit 0
}

# ========== 解析命令行 ==========
specs=()
runs=()
all_runs=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -p) specs+=("$2"); shift 2 ;;
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

# 确定 run 列表
if $all_runs; then
  runs=(1 2 3)
fi
if [ ${#runs[@]} -eq 0 ]; then
  echo "Error: No run specified. Use -r <num> or -r all."
  exit 1
fi

# 确定分析需求
if [ ${#specs[@]} -eq 0 ]; then
  echo "Error: No pair specification. Use -p."
  exit 1
fi

# 检查 VMD 及插件
if ! command -v vmd &> /dev/null; then
  echo "Error: 'vmd' not found in PATH."
  exit 1
fi
if [ ! -f "$VMD_PLUGIN" ]; then
  echo "Error: hbonds plugin not found at $VMD_PLUGIN"
  exit 1
fi

# 检查拓扑文件是否存在
if [ ! -f "$PRMTOP" ]; then
  echo "Error: Topology file '$PRMTOP' not found."
  exit 1
fi

TOPOBASE=$(basename "$PRMTOP" .prmtop)

# 解析一个规范字符串（如 "A:IJ" 或 "AB"）并输出 group1_letters group2_letters
parse_spec() {
  local spec="$1"
  local g1 g2
  if [[ "$spec" == *:* ]]; then
    g1="${spec%%:*}"
    g2="${spec##*:}"
  else
    # 简写格式，如 "AB" → A:B
    if [ ${#spec} -ne 2 ]; then
      echo "Error: Shorthand pair must be exactly 2 letters, got '$spec'. Use 'A:B' for multi-letter groups."
      exit 1
    fi
    g1="${spec:0:1}"
    g2="${spec:1:1}"
  fi
  echo "$g1" "$g2"
}

# 将链字母串转换为 VMD 残基选择表达式，如 "A" → "residue 1 to 97"
# 多链用 or 连接，如 "IJ" → "residue 754 to 899 or residue 900 to 1045"
chainletters_to_vmdsel() {
  local letters="$1"
  local sel=""
  local first=true
  for (( i=0; i<${#letters}; i++ )); do
    chain="${letters:$i:1}"
    range="${CHAIN_RANGE[$chain]}"
    if [ -z "$range" ]; then
      echo "Error: Invalid chain label '$chain'. Allowed: A-J."
      exit 1
    fi
    if $first; then
      sel="residue $range"
      first=false
    else
      sel="$sel or residue $range"
    fi
  done
  echo "($sel)"
}

overall_success=true

for run in "${runs[@]}"; do
  # 轨迹文件路径
  traj_file="${TRAJ_TEMPLATE//\{RUN\}/$run}"
  if [ ! -f "$traj_file" ]; then
    echo "Warning: Trajectory '$traj_file' not found, skipping run $run."
    continue
  fi

  # 输出目录为 $OUTDIR_BASE/run{RUN}/
  outdir_run="${OUTDIR_BASE}/run${run}"
  mkdir -p "$outdir_run"

  for spec in "${specs[@]}"; do
    # 解析出两个组的字母
    read g1_letters g2_letters < <(parse_spec "$spec")

    # 生成 VMD 选择字符串
    sel1=$(chainletters_to_vmdsel "$g1_letters")
    sel2=$(chainletters_to_vmdsel "$g2_letters")

    # 生成文件标签，如 chainA-chainIJ
    label1="chain"
    for (( i=0; i<${#g1_letters}; i++ )); do
      label1+="${g1_letters:$i:1}"
    done
    label2="chain"
    for (( i=0; i<${#g2_letters}; i++ )); do
      label2+="${g2_letters:$i:1}"
    done
    tag="${label1}-${label2}"

    outbase="${TOPOBASE}_run${run}_${tag}"
    detail="${outdir_run}/${outbase}-detail.dat"
    log="${outdir_run}/${outbase}.log"
    outfile="${outdir_run}/${outbase}-hbonds.dat"

    echo ">>> Run $run, $spec -> $outfile"

    vmd_script=$(mktemp /tmp/vmd_hbonds_XXXXXX.tcl)
    cat > "$vmd_script" <<EOF
mol delete all
mol load parm7 $PRMTOP
mol addfile $traj_file first $FIRST last $LAST step $STEP waitfor all

set sel1 [atomselect top "$sel1"]
set sel2 [atomselect top "$sel2"]

source $VMD_PLUGIN
hbonds -sel1 \$sel1 -sel2 \$sel2 -dist $DIST_CUTOFF -ang $ANGLE_CUTOFF -plot no -type pair -writefile yes -detailout $detail -log $log -outfile $outfile

exit
EOF

    vmd -dispdev text -e "$vmd_script"
    if [ $? -ne 0 ]; then
      echo "Error: VMD failed for run $run, spec $spec."
      overall_success=false
    fi

    rm -f "$vmd_script"
  done
done

if $overall_success; then
  echo "All hydrogen bond analyses completed successfully."
else
  echo "Some analyses failed. Check the messages above."
  exit 1
fi
