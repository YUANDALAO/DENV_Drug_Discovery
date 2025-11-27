#!/bin/bash

echo "=== TensorBoard 局域网分享设置 ==="
echo ""

# 停止现有TensorBoard
echo "1. 停止现有TensorBoard..."
pkill tensorboard
sleep 2

# 启动TensorBoard
echo "2. 启动TensorBoard（允许外部访问）..."
tensorboard --logdir experiments/runs/run14a/tensorboard_0 --port 6008 --host 0.0.0.0 > /dev/null 2>&1 &
sleep 3

# 获取IP地址
echo "3. 获取IP地址..."
IP=$(hostname -I | awk '{print $1}')

# 检查端口
echo "4. 检查端口状态..."
if netstat -tlnp 2>/dev/null | grep -q 6008; then
    echo "✓ 端口6008正在监听"
else
    echo "✗ 端口6008未监听，可能启动失败"
    exit 1
fi

# 显示分享信息
echo ""
echo "================================"
echo "✓ TensorBoard已启动！"
echo "================================"
echo ""
echo "本地访问："
echo "  http://localhost:6008"
echo ""
echo "局域网分享给同事："
echo "  http://${IP}:6008"
echo ""
echo "按 Ctrl+C 停止分享"
echo "================================"

# 保持运行
tail -f /dev/null
