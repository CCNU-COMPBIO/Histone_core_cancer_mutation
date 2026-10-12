#!/bin/bash
# 功能: 从原始拓扑中提取单链或任意链组合的拓扑文件
# 用法示例:
#   ./extract_chains.sh -c A -c I          # 单链A和I
#   ./extract_chains.sh -p AB -p EF        # 蛋白链对
#   ./extract_chains.sh -p IJ              # DNA双链
#   ./extract_chains.sh -p AIJ             # 蛋白A + DNA双链
#   ./extract_chains.sh -a                 # 所有常用链组 (AB, CD, EF, GH, IJ)

INPUT="WT_rw.prmtop"
PREFIX="${INPUT%.prmtop}"

# 链残基范围 (新增 I, J)
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
  -c CHAIN      Extract single chain (A-J). Repeatable.
  -p GROUP      Extract chain combination (e.g., AB, IJ, AIJ). Repeatable.
  -a            Extract all predefined groups (AB, CD, EF, GH, IJ).
  -h            Show this help.

Examples:
  $0 -c A -c E                     # Single chains A and E
  $0 -p AB -p EF                   # Protein pairs AB and EF
  $0 -p IJ -p AIJ                  # DNA and protein+DNA
  $0 -a                            # All common groups
  $0 -a -c A                       # All groups plus single chain A
EOF
  exit 0
}

# 解析参数
chains=()
groups=()
all_groups=false

while getopts "c:p:ah" opt; do
  case $opt in
    c) chains+=("$OPTARG") ;;
    p) groups+=("$OPTARG") ;;
    a) all_groups=true ;;
    h) usage ;;
    *) usage ;;
  esac
done

# 如果 -a，预设所有常用链组
if $all_groups; then
  groups=("AB" "CD" "EF" "GH" "IJ")
fi

# 检查是否至少要求了一种提取
if [ ${#chains[@]} -eq 0 ] && [ ${#groups[@]} -eq 0 ]; then
  echo "Error: No chain or group specified. Use -c, -p, or -a."
  usage
fi

# 检查输入文件
if [ ! -f "$INPUT" ]; then
  echo "Error: Topology file '$INPUT' not found."
  exit 1
fi

# 检查 cpptraj
if ! command -v cpptraj &> /dev/null; then
  echo "Error: 'cpptraj' command not found. Please source AmberTools."
  exit 1
fi

overall_success=true

# ---------- 处理所有单链 ----------
for chain in "${chains[@]}"; do
  range="${CHAIN_RANGE[$chain]}"
  if [ -z "$range" ]; then
    echo "Error: Invalid chain label '$chain'. Allowed: A-J."
    exit 1
  fi
  out="${PREFIX}_chain${chain}.prmtop"
  echo ">>> Chain $chain (residues $range) -> $out"

  cpptraj <<EOF
parm $INPUT
parmstrip !(:$range)
parmwrite out $out
EOF

  if [ $? -ne 0 ]; then
    echo "Error generating $out"
    overall_success=false
  fi
done

# ---------- 处理所有链组 ----------
for group in "${groups[@]}"; do
  # 构建残基范围字符串与文件名部分
  range_parts=()
  name_parts=()

  for (( i=0; i<${#group}; i++ )); do
    chain="${group:$i:1}"
    if [ -z "${CHAIN_RANGE[$chain]}" ]; then
      echo "Error: Invalid chain label '$chain' in group '$group'. Allowed: A-J."
      exit 1
    fi
    range_parts+=("${CHAIN_RANGE[$chain]}")
    name_parts+=("chain${chain}")
  done

  range_str=$(IFS=,; echo "${range_parts[*]}")
  name_str=$(IFS=_; echo "${name_parts[*]}")
  out="${PREFIX}_${name_str}.prmtop"

  echo ">>> Group $group (residues $range_str) -> $out"

  cpptraj <<EOF
parm $INPUT
parmstrip !(:$range_str)
parmwrite out $out
EOF

  if [ $? -ne 0 ]; then
    echo "Error generating $out"
    overall_success=false
  fi
done

if $overall_success; then
  echo "All requested topology files have been generated successfully."
else
  echo "Some files failed. Check the messages above."
  exit 1
fi
