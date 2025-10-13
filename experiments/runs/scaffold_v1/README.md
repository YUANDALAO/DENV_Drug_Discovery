# Scaffold-Constrained Generation v1

## 实验信息
- **日期**: 2025-10-04
- **目标**: 生成含吡咯烷或环丁烷骨架的DENV抑制剂
- **配置文件**: config.toml

## 关键参数
- Batch size: 64
- Max steps: 1000
- 评分函数: geometric_mean
- 骨架权重: Pyrrolidine(0.6) + Cyclobutane(0.6)
- 活性权重: 0.4

## 结果摘要
- 总生成分子: 64,000
- 最高Score: 0.847
- 平均pIC50: 7.65
- 骨架匹配问题: 35000+分子被识别为同时含两个骨架（可能是SMARTS误判）

## 主要发现
1. 训练稳定收敛
2. 活性集中在7.5-8.0，缺少高活性分子
3. 生成的主要是桥环系统，不是简单的5元/4元环
4. SMARTS定义需要修正

## 下一步
- [ ] 修正SMARTS表达式（添加羰基约束）
- [ ] 改用arithmetic_mean实现OR逻辑
- [ ] 提高活性权重到1.5
