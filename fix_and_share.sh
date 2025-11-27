#!/bin/bash
set -e

echo "=== 修复并启动TensorBoard分享 ==="
echo ""

# 1. 停止现有TensorBoard
echo "1. 停止现有TensorBoard..."
pkill tensorboard 2>/dev/null || true
sleep 2

# 2. 启动TensorBoard（正确配置）
echo "2. 启动TensorBoard..."
nohup tensorboard --logdir experiments/runs/run14a/tensorboard_0 \
    --port 6008 \
    --host 0.0.0.0 \
    --load_fast=false > tensorboard.log 2>&1 &

sleep 5

# 3. 检查
echo "3. 检查端口..."
if netstat -tlnp 2>/dev/null | grep -q "0.0.0.0:6008" || ss -tlnp 2>/dev/null | grep -q "0.0.0.0:6008"; then
    echo "✓ 端口绑定成功"
else
    echo "✗ 端口绑定失败！"
    echo "查看日志："
    tail -20 tensorboard.log
    exit 1
fi

# 4. 获取IP
IP=$(hostname -I | awk '{print $1}')

# 5. 测试连接
echo "4. 测试连接..."
if curl -s -o /dev/null -w "%{http_code}" http://$IP:6008 | grep -q 200; then
    echo "✓ 连接测试成功"
else
    echo "⚠️  本地可能需要配置防火墙"
fi

echo ""
echo "================================"
echo "✓✓✓ TensorBoard已启动 ✓✓✓"
echo "================================"
echo ""
echo "您访问（本机）："
echo "  http://localhost:6008"
echo ""
echo "同事访问（同一网络）："
echo "  http://$IP:6008"
echo ""
echo "如果同事无法访问，请运行："
echo "  cat tensorboard.log"
echo "  查看错误日志"
echo "================================"
