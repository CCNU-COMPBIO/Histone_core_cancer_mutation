#!/bin/bash
# 功能：批量运行 MMPBSA.py 计算，针对不同 run 依次执行
# 用法：./mmgbsa_batch.sh -p WT_rw -r A -l IJ [-n "1 2 3"] [-i mmgbsa.in]

set -e  # 遇到错误立即退出（可选）

# ========== 默认配置 ==========
PREFIX="WT_rw"                         # 系统前缀
RECEPTOR=""                            # 受体链（如 A）
LIGAND=""                              # 配体链组合（如 IJ）
RUNS=(1 2 3)                           # 默认三个 run
MMGBSA_IN="mmgbsa.in"                  # 输入参数文件
MMPBSA_CMD="MMPBSA.py"                 # 可执行文件
TOP_DIR=".."                           # 拓扑文件所在目录（相对路径）
TRAJ_DIR=".."                          # 轨迹文件所在目录

usage() {
  cat <<EOF
Usage: $0 -p <prefix> -r <receptor_chain> -l <ligand_chains> [options]

Required:
  -p PREFIX        System prefix (e.g., WT_rw)
  -r CHAIN         Receptor chain letter (e.g., A)
  -l CHAINS        Ligand chain letters (e.g., IJ)

Optional:
  -n "RUNS"        Run numbers, quoted space-separated (default: "1 2 3")
  -i INPUT         MMPBSA input file (default: mmgbsa.in)
  --topdir DIR     Directory containing topology files (default: ..)
  --trajdir DIR    Directory containing trajectory files (default: ..)
  --mpirun         Use mpirun (e.g., 'mpirun -np 4')
  --dry-run        Print commands without executing
  -h               Show this help

Examples:
  $0 -p WT_rw -r A -l IJ
  $0 -p E73K_rw -r E -l F -n "1 2"
  $0 -p WT_rw -r A -l IJ --mpirun "mpirun -np 8"
EOF
  exit 0
}

# ========== 解析命令行 ==========
MPIRUN=""
DRYRUN=false

while [[ $# -gt 0 ]]; do
  case $1 in
    -p) PREFIX="$2"; shift 2 ;;
    -r) RECEPTOR="$2"; shift 2 ;;
    -l) LIGAND="$2"; shift 2 ;;
    -n) IFS=' ' read -r -a RUNS <<< "$2"; shift 2 ;;
    -i) MMGBSA_IN="$2"; shift 2 ;;
    --topdir) TOP_DIR="$2"; shift 2 ;;
    --trajdir) TRAJ_DIR="$2"; shift 2 ;;
    --mpirun) MPIRUN="$2"; shift 2 ;;
    --dry-run) DRYRUN=true; shift ;;
    -h) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# 检查必要参数
if [ -z "$RECEPTOR" ] || [ -z "$LIGAND" ]; then
  echo "Error: Both receptor (-r) and ligand (-l) must be specified."
  usage
fi

# 构建链标识字符串 (用于文件名和输出目录)
COMPLEX_NAME="chain${RECEPTOR}_chain${LIGAND}"
RECEPTOR_NAME="chain${RECEPTOR}"
LIGAND_NAME="chain${LIGAND}"

# 输出主目录
MAIN_OUTDIR="mmgbsa_${COMPLEX_NAME}"

# 构建文件名模板
COMPLEX_PRM="${PREFIX}_${COMPLEX_NAME}.prmtop"
RECEPTOR_PRM="${PREFIX}_${RECEPTOR_NAME}.prmtop"
LIGAND_PRM="${PREFIX}_${LIGAND_NAME}.prmtop"
TRAJ_TEMPLATE="${PREFIX}_run{RUN}_500ns_1ns_${COMPLEX_NAME}.nc"

# 检查输入文件
if [ ! -f "${TOP_DIR}/${COMPLEX_PRM}" ]; then
  echo "Error: Complex topology not found: ${TOP_DIR}/${COMPLEX_PRM}"
  exit 1
fi
if [ ! -f "${TOP_DIR}/${RECEPTOR_PRM}" ]; then
  echo "Error: Receptor topology not found: ${TOP_DIR}/${RECEPTOR_PRM}"
  exit 1
fi
if [ ! -f "${TOP_DIR}/${LIGAND_PRM}" ]; then
  echo "Error: Ligand topology not found: ${TOP_DIR}/${LIGAND_PRM}"
  exit 1
fi
if [ ! -f "$MMGBSA_IN" ]; then
  echo "Error: MMPBSA input file not found: $MMGBSA_IN"
  exit 1
fi

# 创建主输出目录
mkdir -p "$MAIN_OUTDIR"

echo "=============================================="
echo "  MMPBSA.py Batch Run"
echo "  System:   $PREFIX"
echo "  Receptor: $RECEPTOR  Ligand: $LIGAND"
echo "  Runs:     ${RUNS[*]}"
echo "  Output:   $MAIN_OUTDIR/"
echo "=============================================="

# ========== 执行每个 run ==========
overall_success=true

for run in "${RUNS[@]}"; do
  # 轨迹文件
  TRAJ_FILE="${TRAJ_DIR}/${TRAJ_TEMPLATE//\{RUN\}/$run}"
  if [ ! -f "$TRAJ_FILE" ]; then
    echo "Warning: Trajectory file '$TRAJ_FILE' not found, skipping run $run."
    continue
  fi

  # 为当前 run 创建子目录
  RUN_DIR="${MAIN_OUTDIR}/run${run}"
  mkdir -p "$RUN_DIR"

  echo ">>> Running MMPBSA.py for run $run ..."

  # 构建 MMPBSA.py 命令
  CMD="$MMPBSA_CMD -O -i ../${MMGBSA_IN} \
    -cp ../${TOP_DIR}/${COMPLEX_PRM} \
    -rp ../${TOP_DIR}/${RECEPTOR_PRM} \
    -lp ../${TOP_DIR}/${LIGAND_PRM} \
    -y ../${TRAJ_FILE}"

  # 如果使用 mpirun，则加上
  if [ -n "$MPIRUN" ]; then
    CMD="$MPIRUN $CMD"
  fi

  echo "  Directory: $RUN_DIR"
  echo "  Command:   $CMD"

  if $DRYRUN; then
    echo "  (Dry run, skipping execution)"
  else
    # 进入子目录执行，使输出文件留在 run 子目录内
    ( cd "$RUN_DIR" && eval "$CMD" )
    if [ $? -ne 0 ]; then
      echo "Error: MMPBSA.py failed for run $run."
      overall_success=false
    else
      echo "  Run $run completed successfully."
    fi
  fi
done

if $overall_success; then
  echo "All MMPBSA calculations finished."
else
  echo "Some calculations failed. Check the output above."
  exit 1
fi
