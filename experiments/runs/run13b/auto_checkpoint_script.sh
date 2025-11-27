#!/bin/bash
# 自动checkpoint监控脚本
# 功能：每2000步自动Ctrl+C保存，然后自动继续训练

WORK_DIR="experiments/runs/run13b"
CONFIG_FILE="$WORK_DIR/config.toml"
LOG_FILE="$WORK_DIR/training.log"
PID_FILE="$WORK_DIR/training.pid"
CHECKPOINT_INTERVAL=2000  # 每2000步保存一次

echo "=== 自动Checkpoint监控启动 ==="
echo "工作目录: $WORK_DIR"
echo "Checkpoint间隔: ${CHECKPOINT_INTERVAL}步"
echo "时间: $(date)"
echo ""

# 阶段1: 0-2000步
echo "[Stage 1/3] 启动训练: 0-2000步"
cd "$WORK_DIR"
nohup reinvent config.toml > training.log 2>&1 &
TRAIN_PID=$!
echo $TRAIN_PID > training.pid
echo "训练PID: $TRAIN_PID"

# 监控到2000步
while true; do
    sleep 30
    
    # 检查进程是否还在
    if ! ps -p $TRAIN_PID > /dev/null 2>&1; then
        echo "训练进程意外停止！"
        exit 1
    fi
    
    # 获取当前step
    CURRENT_STEP=$(tail -50 training.log | grep -oP 'Step: \K\d+' | tail -1)
    
    if [ -n "$CURRENT_STEP" ] && [ "$CURRENT_STEP" -ge 2000 ]; then
        echo ""
        echo "[Stage 1完成] 达到2000步，保存checkpoint..."
        
        # 发送SIGINT保存checkpoint
        kill -INT $TRAIN_PID
        sleep 10
        
        # 重命名checkpoint
        if [ -f "$WORK_DIR/checkpoint_6k.chkpt" ]; then
            cp "$WORK_DIR/checkpoint_6k.chkpt" "$WORK_DIR/checkpoint_2k_backup.chkpt"
            echo "✓ Checkpoint已保存: checkpoint_2k_backup.chkpt"
        fi
        break
    fi
    
    echo -ne "\r当前Step: $CURRENT_STEP / 2000"
done

echo ""
echo "=== 阶段1完成，准备继续训练 ==="
sleep 5

# 阶段2: 2000-4000步
echo ""
echo "[Stage 2/3] 继续训练: 2000-4000步"

# 修改配置使用checkpoint
sed -i.bak 's/agent_file = "priors\/denv_libinvent_model_v2.model"/agent_file = "experiments\/runs\/run13b\/checkpoint_6k.chkpt"/' config.toml

nohup reinvent config.toml >> training.log 2>&1 &
TRAIN_PID=$!
echo $TRAIN_PID > training.pid
echo "训练PID: $TRAIN_PID"

# 监控到4000步
TOTAL_STEPS=0
while true; do
    sleep 30
    
    if ! ps -p $TRAIN_PID > /dev/null 2>&1; then
        echo "训练进程意外停止！"
        exit 1
    fi
    
    CURRENT_STEP=$(tail -50 training.log | grep -oP 'Step: \K\d+' | tail -1)
    TOTAL_STEPS=$((2000 + CURRENT_STEP))
    
    if [ -n "$CURRENT_STEP" ] && [ "$CURRENT_STEP" -ge 2000 ]; then
        echo ""
        echo "[Stage 2完成] 达到4000步（总计），保存checkpoint..."
        
        kill -INT $TRAIN_PID
        sleep 10
        
        if [ -f "$WORK_DIR/checkpoint_6k.chkpt" ]; then
            cp "$WORK_DIR/checkpoint_6k.chkpt" "$WORK_DIR/checkpoint_4k_backup.chkpt"
            echo "✓ Checkpoint已保存: checkpoint_4k_backup.chkpt"
        fi
        break
    fi
    
    echo -ne "\r当前Step: $TOTAL_STEPS / 4000"
done

echo ""
echo "=== 阶段2完成，最后冲刺 ==="
sleep 5

# 阶段3: 4000-6000步
echo ""
echo "[Stage 3/3] 最终训练: 4000-6000步"

nohup reinvent config.toml >> training.log 2>&1 &
TRAIN_PID=$!
echo $TRAIN_PID > training.pid
echo "训练PID: $TRAIN_PID"

# 监控到6000步（或自然结束）
while true; do
    sleep 30
    
    if ! ps -p $TRAIN_PID > /dev/null 2>&1; then
        echo ""
        echo "✓ 训练自然完成！"
        break
    fi
    
    CURRENT_STEP=$(tail -50 training.log | grep -oP 'Step: \K\d+' | tail -1)
    TOTAL_STEPS=$((4000 + CURRENT_STEP))
    
    echo -ne "\r当前Step: $TOTAL_STEPS / 6000"
    
    if [ -n "$CURRENT_STEP" ] && [ "$CURRENT_STEP" -ge 2000 ]; then
        echo ""
        echo "✓ 达到6000步目标！"
        break
    fi
done

echo ""
echo "==="
echo "🎉 Run13b 完整训练完成！"
echo "==="
echo "结果文件:"
echo "  - results_1.csv"
echo "  - checkpoint_6k.chkpt (最终)"
echo "  - checkpoint_2k_backup.chkpt"
echo "  - checkpoint_4k_backup.chkpt"
echo ""
echo "时间: $(date)"
