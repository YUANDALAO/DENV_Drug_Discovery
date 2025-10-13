#!/bin/bash
#SBATCH --job-name=run3_scaffold
#SBATCH --output=logs/run3_%j.out
#SBATCH --error=logs/run3_%j.err
#SBATCH --time=06:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

# ========================================
# Run3: 骨架约束强化验证
# 目标：生成100个符合严格骨架定义的分子
# ========================================

echo "=========================================="
echo "Run3 启动时间: $(date)"
echo "=========================================="

# 激活环境
source ~/.bashrc
conda activate reinvent4_env

# 工作目录
cd /mnt/c/Users/ucsaheu/python_projects/DENV_Drug_Discovery/05_Generative_AI_REINVENT/REINVENT4-main/experiments/runs/run3

# 运行REINVENT4
echo "开始骨架约束生成..."
reinvent \
    configs/config.toml \
    -l logs/run3_reinvent.log

# 检查运行状态
if [ $? -eq 0 ]; then
    echo "✅ Run3 成功完成"
    echo "结束时间: $(date)"
else
    echo "❌ Run3 运行失败"
    exit 1
fi

# 输出文件统计
echo ""
echo "=========================================="
echo "生成文件统计："
echo "=========================================="
ls -lh results/csv/run3_scaffold_strict*.csv 2>/dev/null || echo "未找到CSV文件"
ls -lh checkpoints/*.chkpt 2>/dev/null || echo "未找到检查点文件"

echo ""
echo "=========================================="
echo "Run3 任务完成"
echo "=========================================="
