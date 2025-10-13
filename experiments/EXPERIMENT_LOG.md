# REINVENT4 实验记录

## 实验列表

### scaffold_v1_original (2025-10-04)
- **状态**: 已完成
- **配置**: 原始骨架约束（可能有SMARTS问题）
- **结果**: Score 0.847, pIC50 avg 7.65
- **问题**: 生成桥环而非简单环
- **路径**: `runs/scaffold_v1_original/`

### scaffold_v2_corrected_smarts (计划中)
- **状态**: 待运行
- **改进**: 修正SMARTS，arithmetic_mean
- **路径**: `runs/scaffold_v2_corrected_smarts/`

## 配置版本历史
- v1: scaffold_FIXED.toml - 原始版本，SMARTS不含羰基
- v2: scaffold_PI_CONFIRMED.toml - PI确认版本，添加羰基约束
