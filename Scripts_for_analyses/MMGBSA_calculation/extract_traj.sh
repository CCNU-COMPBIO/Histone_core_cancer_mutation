#!/bin/bash
# 使用说明
# 提取特定链对应的轨迹：-p AB
# 提取run1对应的轨迹：-r 1
# ./extract_traj.sh -p AB -p EF -r 1 -r 2
# 提取 DNA + 组蛋白（所有 run）
# ./extract_traj.sh -p AIJ -r all

# ===== 默认配置 =====
PRMTOP="WT_rw.prmtop"
TRAJ_TEMPLATE="./run{RUN}/WT_rw_run{RUN}_500ns.nc"
STRIDE=50                                 # 提取步长，默认每 50 帧取 1 帧（即每 1 ns 一帧）
FIRST_FRAME=5000                          # 起始帧（默认 100 ns，对应 5000 帧）
LAST_FRAME=25000                          # 终止帧（默认 500 ns，对应 25000 帧）
OUT_TEMPLATE="{PREFIX}_run{RUN}_500ns_1ns_{PAIR}.nc"

PREFIX="${PRMTOP%.prmtop}"

declare -A CHAIN_RANGE=(
  [A]="1-97"
  [B]="98-175"
  [C]="176-280"
  [D]="281-376"
  [E]="377-475"
  [F]="476-557"
  [G]="558-661"
  [H]="662-753"
  [I]="754-899"
  [J]="900-1045"
)

usage() {
  cat <<EOF
Usage: $0 [options]

Options:
  -i <prmtop>        Topology file (default: $PRMTOP)
  -t <template>      Trajectory path template, use {RUN} as placeholder
                     (default: $TRAJ_TEMPLATE)
  -s <stride>        Frame stride for trajin (default: $STRIDE)
  --first <frame>    First frame to read (default: $FIRST_FRAME)
  --last <frame>     Last frame to read (default: $LAST_FRAME)
  -o <template>      Output name template, use {PREFIX}, {RUN}, {PAIR}
                     (default: $OUT_TEMPLATE)
  -p <group>         Chain combination to extract (e.g., AB, IJ, AIJ).
                     Repeatable. Allowed chain letters: A-J.
  -r <run>           Run number (1, 2, 3). Repeatable, or '-r all' for all.
  -a                 Extract all predefined pairs (AB, CD, EF, GH, IJ).
  -h                 Show this help.

Examples:
  $0 -p AB -p EF -r 1 -r 2        # Protein pairs AB & EF
  $0 -p IJ -r all                  # DNA double strand
  $0 -p AIJ -r 1 -r 2              # Protein chain A + DNA
  $0 -a -r all                     # All pairs (including IJ) from all runs
EOF
}

# ===== 解析命令行 =====
groups=()
runs=()
all_groups=false
all_runs=false

while getopts "i:t:s:o:p:r:ah-:" opt; do
  case $opt in
    -)
      case "${OPTARG}" in
        first) FIRST_FRAME="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 )) ;;
        last) LAST_FRAME="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 )) ;;
        first=*) FIRST_FRAME="${OPTARG#*=}" ;;
        last=*) LAST_FRAME="${OPTARG#*=}" ;;
        *) usage; exit 1 ;;
      esac
      ;;
    i) PRMTOP="$OPTARG"
       PREFIX="${PRMTOP%.prmtop}"
       ;;
    t) TRAJ_TEMPLATE="$OPTARG" ;;
    s) STRIDE="$OPTARG" ;;
    o) OUT_TEMPLATE="$OPTARG" ;;
    p) groups+=("$OPTARG") ;;
    r) if [ "$OPTARG" = "all" ]; then
         all_runs=true
       else
         runs+=("$OPTARG")
       fi ;;
    a) all_groups=true ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# 如果 -a，则添加所有预设对（包括新增的 IJ）
if $all_groups; then
  groups=("AB" "CD" "EF" "GH" "IJ")
fi
if [ ${#groups[@]} -eq 0 ]; then
  echo "Error: No chain group specified. Use -p or -a."
  usage
  exit 1
fi

# 确定 run 列表
if $all_runs; then
  runs=(1 2 3)
elif [ ${#runs[@]} -eq 0 ]; then
  echo "Error: No run specified. Use -r or -r all."
  usage
  exit 1
fi

# 检查输入拓扑
if [ ! -f "$PRMTOP" ]; then
  echo "Error: Topology file '$PRMTOP' not found."
  exit 1
fi

# 检查 cpptraj 是否可用
if ! command -v cpptraj &> /dev/null; then
  echo "Error: 'cpptraj' command not found in PATH."
  exit 1
fi

# ===== 执行任务 =====
overall_success=true

for run in "${runs[@]}"; do
  traj_file="${TRAJ_TEMPLATE//\{RUN\}/$run}"
  if [ ! -f "$traj_file" ]; then
    echo "Warning: Trajectory file '$traj_file' not found, skipping run $run."
    continue
  fi

  for group in "${groups[@]}"; do
    # 构建残基范围字符串，例如 "1-97,98-175"
    range_parts=()
    # 构建文件名中的链标签部分，例如 "chainA_chainB"
    name_parts=()
    
    # 逐字符检查并收集
    for (( i=0; i<${#group}; i++ )); do
      chain="${group:$i:1}"
      if [ -z "${CHAIN_RANGE[$chain]}" ]; then
        echo "Error: Invalid chain label '$chain' in group '$group'. Allowed: A-J."
        exit 1
      fi
      range_parts+=("${CHAIN_RANGE[$chain]}")
      name_parts+=("chain$chain")
    done

    # 拼接范围（用逗号分隔）
    range_str=$(IFS=,; echo "${range_parts[*]}")
    # 拼接文件名部分（用下划线分隔）
    name_str=$(IFS=_; echo "${name_parts[*]}")

    # 生成输出文件名，替换 {PAIR} 为 name_str
    outfile="$OUT_TEMPLATE"
    outfile="${outfile//\{PREFIX\}/$PREFIX}"
    outfile="${outfile//\{RUN\}/$run}"
    outfile="${outfile//\{PAIR\}/$name_str}"

    echo ">>> Run $run, group $group -> $outfile (frames $FIRST_FRAME to $LAST_FRAME, stride $STRIDE)"

    cpptraj <<EOF
parm $PRMTOP
trajin $traj_file $FIRST_FRAME $LAST_FRAME $STRIDE
strip !(:$range_str)
trajout $outfile
EOF

    if [ $? -ne 0 ]; then
      echo "Error: Failed to process run $run, group $group."
      overall_success=false
    fi
  done
done

if $overall_success; then
  echo "All trajectory extractions completed successfully."
else
  echo "Some extractions failed. Check the messages above."
  exit 1
fi
