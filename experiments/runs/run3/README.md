# REINVENT4 Run3 - DENV抑制剂生成

## 基本信息
- **日期**: 2025-10-05
- **目标**: 生成含吡咯烷型或环丁烷型骨架的DENV抑制剂
- **评分函数**: arithmetic_mean (OR逻辑)

## 骨架定义

### 吡咯烷型
```
SMILES: CC1CC(C(N1)c1ccccc1)C(N)=O
SMARTS: [#6]~[#6]1~[#6]~[#6](-[#6](=[#8])[#7])~[#7]~[#6]~1[c,n,o,s]
```

### 环丁烷型
```
SMILES: CC1CC(C1c1ccccc1)C(N)=O
SMARTS: [#6]~[#6]1~[#6]~[#6](-[#6](=[#8])[#7])~[#6]~1[c,n,o,s]
```

## 预期结果
- **理论生成**: 128,000个分子 (64 × 2000)
- **预计输出**: 60,000-90,000个唯一分子
- **骨架匹配率目标**: >80%
- **运行时间**: 2-3小时 (GPU)

## 监控命令
```bash
# 查看日志
tail -f experiments/runs/run3/logs/training_*.log

# TensorBoard
tensorboard --logdir=experiments/runs/run3/tensorboard
```
