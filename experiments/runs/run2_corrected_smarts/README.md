# Run 2: Corrected SMARTS - Scaffold-Constrained Generation

**日期**: 2025-10-04  
**状态**: 准备运行

## 改进点（相比Run1）

### 1. SMARTS修正
- **Run1**: `c1ccccc1C1CNCC1` (缺少酰胺基)
- **Run2**: `C1NC(c2ccccc2)C(C(=O)N)C1` (包含酰胺基)

### 2. 评分逻辑
- **Run1**: geometric_mean (AND逻辑 - 需同时满足)
- **Run2**: arithmetic_mean (OR逻辑 - 满足任一即可)

### 3. 权重调整
- 活性权重: 0.4 → 1.5
- 骨架权重: 保持1.0（但在arithmetic下影响不同）

## 预期结果
- 生成简单的5元吡咯烷或4元环丁烷结构
- 两个骨架各占约50%
- R基团多样化
- 避免Run1中的桥环问题

## 文件说明
- `config.toml` - 训练配置
- `results.csv` - 训练结果（待生成）
- `logs/` - 训练日志
- `checkpoints/` - 训练检查点
