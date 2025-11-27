#!/bin/bash
echo "启动TensorBoard..."
pkill tensorboard 2>/dev/null || true
sleep 2

nohup tensorboard --logdir experiments/runs/run14a/tensorboard_0 \
    --port 6008 \
    --host 0.0.0.0 \
    --load_fast=false > tensorboard.log 2>&1 &

sleep 3

IP=$(hostname -I | awk '{print $1}')
echo ""
echo "✓ TensorBoard已启动"
echo "访问地址: http://$IP:6008"
echo ""
