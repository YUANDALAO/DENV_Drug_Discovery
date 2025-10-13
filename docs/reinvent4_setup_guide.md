# REINVENT4 DENV项目设置指南

## 成功的配置
- 数据：100个简单分子作为测试集
- 模型：priors/reinvent.prior
- 参数：batch_size=32, sample_batch_size=100, num_epochs=2

## 重要发现
- 需要严格的数据清理，过滤包含B、复杂离子的分子
- sample_batch_size 必须 ≥ 100
- 使用简单分子结构可以避免token兼容性问题

## 下一步
- 扩展到更大的数据集
- 进行完整的100 epoch训练
- 测试生成新分子
