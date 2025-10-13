#!/bin/bash
# Run3 训练启动脚本

cd /mnt/c/Users/ucsaheu/python_projects/DENV_Drug_Discovery/05_Generative_AI_REINVENT/REINVENT4-main

# 激活环境
source activate reinvent4

# 启动训练
nohup python input.py experiments/runs/run3/config.toml \
    > experiments/runs/run3/logs/training_$(date +%Y%m%d_%H%M%S).log 2>&1 &

TRAIN_PID=$!
echo "✅ 训练已启动 (PID: $TRAIN_PID)"
echo "📋 日志: experiments/runs/run3/logs/training_*.log"
echo ""
echo "监控命令:"
echo "  tail -f experiments/runs/run3/logs/training_*.log"
echo ""
echo "停止训练:"
echo "  kill $TRAIN_PID"
