#!/bin/bash
echo "=== 完整网络诊断 ==="
echo ""

IP=$(hostname -I | awk '{print $1}')

echo "【WSL网络信息】"
echo "IP地址: $IP"
ip route | grep default
echo ""

echo "【TensorBoard状态】"
ps aux | grep tensorboard | grep -v grep
netstat -tlnp 2>/dev/null | grep 6008
echo ""

echo "【测试本地访问】"
curl -I http://localhost:6008 2>&1 | head -5
echo ""

echo "【测试IP访问】"
curl -I http://$IP:6008 2>&1 | head -5
echo ""

echo "【WSL防火墙】"
sudo iptables -L -n | grep 6008 || echo "无相关规则"
echo ""

echo "================================"
echo "Windows端检查（在PowerShell中运行）："
echo "================================"
echo "1. 检查防火墙规则："
echo "   Get-NetFirewallRule -DisplayName 'TensorBoard-6008'"
echo ""
echo "2. 检查端口转发："
echo "   netsh interface portproxy show all"
echo ""
echo "3. 检查网络配置文件："
echo "   Get-NetConnectionProfile"
echo ""
echo "4. 测试端口："
echo "   Test-NetConnection -ComputerName $IP -Port 6008"
