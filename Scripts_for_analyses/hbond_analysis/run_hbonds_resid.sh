#!/bin/bash
# ============================================================
# 脚本名: run_hbonds_resid.sh
# 功能:   计算指定残基（例如 D 链的 70 号残基）与整个系统其余部分
#         （蛋白其他残基 + DNA）之间的氢键。
# 目录结构要求:
#   脚本放在与 WT_rw.prmtop 相同的目录下
#   轨迹文件位于 ./run1/WT_rw_run1_500ns.nc, ./run2/..., ./run3/...
#   结果文件也输出到对应的 run 子目录中
# 依赖:   VMD + hbonds.tcl
# ============================================================
# 在脚本所在目录执行
# ./run_hbonds_resid.sh -s D:70 -r 1
# ./run_hbonds_resid.sh -s D:70 -r all

# ========== 默认参数 ==========
PRMTOP="WT_rw.prmtop"                                    # 拓扑文件（与脚本同级）
TRAJ_TEMPLATE="./run{RUN}/WT_rw_run{RUN}_500ns.nc"      # 轨迹模板，{RUN} 会被替换
FIRST=5000                                               # 起始帧
LAST=25000                                               # 结束帧
STEP=50                                                 # 步长
DIST_CUTOFF=3.5                                          # 氢键距离 cutoff (Å)
ANGLE_CUTOFF=30                                          # 氢键角度 cutoff (度)
VMD_PLUGIN="/home/lenovo/vmd-1.9.4a57/plugins/noarch/tcl/hbonds1.3/hbonds.tcl"
OUTDIR_BASE="./"                                         # 输出基本目录（将在每个 run 子目录中生成文件）

# ========== 链残基全局范围（VMD resid 编号）==========
declare -A CHAIN_RANGE=(
  [A]="1 to 97"
  [B]="98 to 175"
  [C]="176 to 280"
  [D]="281 to 376"
  [E]="377 to 475"
  [F]="476 to 557"
  [G]="558 to 661"
  [H]="662 to 753"
  [I]="754 to 899"   # DNA 链 I
  [J]="900 to 1045"  # DNA 链 J
)

# ========== 每条链的局部起始编号（用于局部→全局转换）==========
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

# ========== 帮助信息 ==========
usage() {
  cat <<EOF
Usage: $0 -s CHAIN:LOCAL_RESID -r RUN [options]

必需参数:
  -s CHAIN:LOCAL_RESID   指定链和局部残基号，例如 -s D:70
  -r RUN                 模拟序号 (1,2,3 或 all)

可选参数:
  -i PRMTOP              拓扑文件路径 (默认: 当前目录下的 $PRMTOP)
  -t TEMPLATE            轨迹模板，其中 {RUN} 会被替换 (默认: $TRAJ_TEMPLATE)
  --first N              起始帧 (默认: $FIRST)
  --last N               结束帧 (默认: $LAST)
  --step N               帧步长 (默认: $STEP)
  --dist D               氢键距离 cutoff (默认: $DIST_CUTOFF)
  --ang A                氢键角度 cutoff (默认: $ANGLE_CUTOFF)
  --vmd-plugin P         hbonds.tcl 路径 (默认: $VMD_PLUGIN)
  --outdir DIR           输出基本目录，实际输出到 DIR/run{RUN}/ (默认: 当前目录)
  -h                     显示本帮助

示例:
  $0 -s D:70 -r 1
  $0 -s D:70 -r all
  $0 -s D:70 -r 1 -i ../WT_rw.prmtop -t ../run{RUN}/WT_rw_run{RUN}_500ns.nc
EOF
  exit 0
}

# ========== 解析命令行参数 ==========
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
    *) echo "未知选项: $1"; usage ;;
  esac
done

# 检查必需参数
if [ -z "$single_spec" ]; then
  echo "错误: 必须使用 -s 指定残基，例如 -s D:70"
  usage
fi

# 确定 run 列表
if $all_runs; then
  runs=(1 2 3)
fi
if [ ${#runs[@]} -eq 0 ]; then
  echo "错误: 必须使用 -r 指定一个或多个 run 序号，或使用 -r all"
  exit 1
fi

# 检查 VMD 及插件
if ! command -v vmd &> /dev/null; then
  echo "错误: 未在 PATH 中找到 vmd 命令"
  exit 1
fi
if [ ! -f "$VMD_PLUGIN" ]; then
  echo "错误: 找不到 hbonds 插件: $VMD_PLUGIN"
  exit 1
fi

# 检查拓扑文件是否存在
if [ ! -f "$PRMTOP" ]; then
  echo "错误: 拓扑文件 '$PRMTOP' 不存在"
  exit 1
fi

TOPOBASE=$(basename "$PRMTOP" .prmtop)

# ========== 辅助函数：局部残基编号 → 全局残基编号 ==========
local_to_global() {
  local chain="$1"
  local local_res="$2"
  local range="${CHAIN_RANGE[$chain]}"
  local local_start="${CHAIN_LOCAL_START[$chain]}"
  if [ -z "$range" ]; then
    echo "错误: 未知链标识 '$chain'" >&2
    exit 1
  fi
  if [ -z "$local_start" ]; then
    echo "错误: 链 '$chain' 未定义局部起始编号，请在脚本中补充 CHAIN_LOCAL_START 数组" >&2
    exit 1
  fi
  # 提取全局起始号 (范围字符串的第一个数字)
  local global_start=$(echo "$range" | awk '{print $1}')
  echo $(( global_start + local_res - local_start ))
}

# ========== 主程序：单残基氢键分析 ==========
overall_success=true

# 解析 -s 参数，例如 "D:70"
if [[ "$single_spec" != *:* ]]; then
  echo "错误: -s 格式必须为 CHAIN:LOCAL_RESID，例如 D:70"
  exit 1
fi
chain="${single_spec%%:*}"
local_resid="${single_spec##*:}"
global_resid=$(local_to_global "$chain" "$local_resid")
echo ">>> 单残基模式: 链 $chain, 局部残基 $local_resid -> 全局残基编号 $global_resid"

# DNA 选择（根据 CHAIN_RANGE 中的 I 和 J 链）
dna_sel="residue ${CHAIN_RANGE[I]} or residue ${CHAIN_RANGE[J]}"
# sel2: 除自身以外的所有蛋白残基 + DNA
sel2_expr="(not (same residue as within 0 of resid $global_resid)) or ($dna_sel)"

for run in "${runs[@]}"; do
  # 生成轨迹文件路径
  traj_file="${TRAJ_TEMPLATE//\{RUN\}/$run}"
  if [ ! -f "$traj_file" ]; then
    echo "警告: 轨迹文件 '$traj_file' 不存在, 跳过 run $run"
    continue
  fi

  # 输出目录为 $OUTDIR_BASE/run{RUN}/
  outdir_run="${OUTDIR_BASE}/run${run}"
  mkdir -p "$outdir_run"

  # 输出文件名使用拓扑名和残基标识
  outbase="${TOPOBASE}_run${run}_resid_${chain}${local_resid}"
  detail="${outdir_run}/${outbase}-detail.dat"
  log="${outdir_run}/${outbase}.log"
  outfile="${outdir_run}/${outbase}-hbonds.dat"

  echo ">>> 正在处理 Run $run, 残基 $single_spec -> 输出: $outfile"

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
    echo "错误: VMD 运行失败, run $run, 残基 $single_spec"
    overall_success=false
  fi
  rm -f "$vmd_script"
done

if $overall_success; then
  echo "所有单残基氢键分析已完成。"
else
  echo "部分分析失败，请检查上述错误信息。"
  exit 1
fi
