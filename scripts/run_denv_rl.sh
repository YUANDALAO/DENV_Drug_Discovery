#!/bin/bash
# DENV抑制剂AI生成 - 强化学习启动脚本
# 作者: Your Name
# 日期: 2025-10-02

set -e  # 遇到错误立即退出

echo "=========================================="
echo "DENV抑制剂AI生成 - 强化学习任务"
echo "=========================================="

# 定义项目路径
PROJECT_ROOT="/mnt/c/Users/ucsaheu/python_projects/DENV_Drug_Discovery"
REINVENT_DIR="${PROJECT_ROOT}/05_Generative_AI_REINVENT/REINVENT4-main"
CONFIG_FILE="${REINVENT_DIR}/configs/final_rl_task.toml"

echo "项目根目录: ${PROJECT_ROOT}"
echo "REINVENT目录: ${REINVENT_DIR}"
echo "配置文件: ${CONFIG_FILE}"
echo ""

# 检查conda是否已安装
if ! command -v conda &> /dev/null; then
    echo "错误: conda未安装或未在PATH中"
    exit 1
fi

# 检查环境是否存在
if ! conda env list | grep -q "reinvent4-gpu"; then
    echo "错误: reinvent4-gpu环境不存在"
    echo "请先运行环境安装命令"
    exit 1
fi

# 激活conda环境
echo "激活conda环境: reinvent4-gpu"
eval "$(conda shell.bash hook)"
conda activate reinvent4-gpu

# 检查必要文件是否存在
echo ""
echo "检查必要文件..."

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "错误: 配置文件不存在: ${CONFIG_FILE}"
    exit 1
fi
echo "✓ 配置文件存在"

if [ ! -f "${REINVENT_DIR}/priors/denv_ultra_clean_model.model" ]; then
    echo "错误: 专家模型不存在"
    exit 1
fi
echo "✓ 专家模型存在"

if [ ! -f "${PROJECT_ROOT}/02_ML_Modeling/models/random_forest_champion.joblib" ]; then
    echo "错误: QSAR冠军模型不存在"
    exit 1
fi
echo "✓ QSAR模型存在"

# 创建必要的目录
echo ""
echo "创建必要目录..."
mkdir -p "${REINVENT_DIR}/scaffolds"
mkdir -p "${REINVENT_DIR}/results"
mkdir -p "${REINVENT_DIR}/tb_denv_rl"

# 创建骨架文件（如果不存在）
SCAFFOLD_FILE="${REINVENT_DIR}/scaffolds/denv_scaffolds.smi"
if [ ! -f "${SCAFFOLD_FILE}" ]; then
    echo "创建骨架文件: ${SCAFFOLD_FILE}"
    cat > "${SCAFFOLD_FILE}" << 'EOF'
C1CCNC1* pyrrolidine
C1CCC1* cyclobutane
EOF
    echo "✓ 骨架文件已创建"
else
    echo "✓ 骨架文件已存在"
fi

# 显示GPU信息
echo ""
echo "GPU信息:"
nvidia-smi --query-gpu=index,name,memory.total,memory.free --format=csv,noheader

# 进入工作目录
cd "${REINVENT_DIR}"
echo ""
echo "当前工作目录: $(pwd)"

# 显示配置文件内容摘要
echo ""
echo "配置文件摘要:"
echo "----------------------------------------"
grep -E "^(run_type|device|prior_file|batch_size|max_steps)" "${CONFIG_FILE}" || true
echo "----------------------------------------"

# 询问是否继续
echo ""
read -p "是否开始运行强化学习任务? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "任务已取消"
    exit 0
fi

# 运行强化学习
echo ""
echo "=========================================="
echo "开始强化学习训练..."
echo "=========================================="
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# 记录开始时间
START_TIME=$(date +%s)

# 运行REINVENT4
python reinvent/main.py "${CONFIG_FILE}" 2>&1 | tee "logs/rl_training_$(date +%Y%m%d_%H%M%S).log"

# 计算运行时间
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(((DURATION % 3600) / 60))
SECONDS=$((DURATION % 60))

echo ""
echo "=========================================="
echo "训练完成!"
echo "=========================================="
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "总耗时: ${HOURS}小时 ${MINUTES}分钟 ${SECONDS}秒"
echo ""
echo "结果文件位置:"
echo "  - CSV结果: ${REINVENT_DIR}/denv_rl_results_*.csv"
echo "  - 检查点: ${REINVENT_DIR}/denv_rl_final.chkpt"
echo "  - TensorBoard日志: ${REINVENT_DIR}/tb_denv_rl/"
echo ""
echo "查看TensorBoard:"
echo "  tensorboard --logdir=${REINVENT_DIR}/tb_denv_rl --port=6006"
echo ""