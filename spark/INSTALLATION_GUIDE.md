# REINVENT4 on NVIDIA DGX Spark - 完整部署指南

## 硬件配置
- **GPU**: GB10 (Grace Blackwell架构)
- **内存**: 128GB unified memory
- **CPU**: 20 ARM cores (Grace CPU)
- **存储**: 4TB
- **系统**: Linux ARM架构 + CUDA 13.0

## 部署步骤

### 第一步：QSAR插件安装

```bash
# 1. 进入REINVENT4目录
cd ~/aidd/REINVENT4

# 2. 运行安装脚本（从Claude提供的文件）
bash qsar_plugin_setup.sh

# 3. 验证插件安装
ls -lh reinvent/reinvent_plugins/components/comp_qsar_scorer.py
```

### 第二步：复制必需文件

```bash
# 创建必要的目录
mkdir -p ~/aidd/REINVENT4/models
mkdir -p ~/aidd/REINVENT4/data
mkdir -p ~/aidd/REINVENT4/priors
mkdir -p ~/aidd/REINVENT4/experiments/runs/spark_run1

# 从旧电脑复制以下文件到新机器（使用scp或其他方式）:
# 1. QSAR模型:
#    random_forest_champion.joblib → ~/aidd/REINVENT4/models/

# 2. Prior和Agent模型:
#    libinvent.prior → ~/aidd/REINVENT4/priors/
#    denv_libinvent_model.model → ~/aidd/REINVENT4/priors/

# 3. Scaffold文件:
#    pyrrolidine_dual_aryl.smi → ~/aidd/REINVENT4/data/
```

### 第三步：验证环境

```bash
# 1. 激活conda环境
conda activate reinvent4

# 2. 检查PyTorch和CUDA
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}')"

# 3. 检查REINVENT4
python -c "import reinvent; print(f'REINVENT4: {reinvent.__version__}')"

# 4. 检查RDKit
python -c "from rdkit import Chem; print('RDKit OK')"
```

### 第四步：运行方式选择

#### 选项A：使用Jupyter Notebook（推荐）

```bash
# 1. 复制notebook到REINVENT目录
cp REINVENT4_DGX_Spark.ipynb ~/aidd/REINVENT4/

# 2. 启动Jupyter
cd ~/aidd/REINVENT4
jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser

# 3. 在浏览器中打开显示的URL
# 4. 打开 REINVENT4_DGX_Spark.ipynb
# 5. 点击 "Run All" 即可开始训练！
```

#### 选项B：使用命令行

```bash
# 1. 复制配置文件
cp spark_optimized_config.toml ~/aidd/REINVENT4/experiments/runs/spark_run1/config.toml

# 2. 进入REINVENT目录
cd ~/aidd/REINVENT4

# 3. 运行REINVENT（后台运行）
nohup reinvent experiments/runs/spark_run1/config.toml > experiments/runs/spark_run1/training.log 2>&1 &

# 4. 查看运行状态
tail -f experiments/runs/spark_run1/training.log

# 5. 查看进程
ps aux | grep reinvent
```

### 第五步：实时监控（可选）

```bash
# 启动TensorBoard
cd ~/aidd/REINVENT4
tensorboard --logdir experiments/runs/spark_run1/tensorboard --bind_all --port 6006

# 在浏览器打开: http://<your-server-ip>:6006
```

## 配置参数说明

### 关键优化参数（vs T1200）

| 参数 | T1200旧配置 | DGX Spark新配置 | 提升倍数 |
|------|-------------|-----------------|----------|
| batch_size | 32 | 512 | 16x |
| max_steps | 1000-2000 | 5000 | 2.5-5x |
| learning_rate | 0.0001 | 0.00015 | 1.5x |
| inception memory | 30 | 100 | 3.3x |
| inception sample | 3 | 20 | 6.7x |
| diversity bucket | 25 | 15 | - |

### 预期性能

```
生成速度: ~30-50 molecules/second (vs T1200的 ~12-15 molecules/second)
总分子数: 2,560,000 (5000 steps × 512 batch)
运行时间: 8-12 hours (vs T1200的 24+ hours)
内存使用: ~80-100GB (峰值)
GPU利用率: 85-95%
```

## 输出文件

训练完成后，在 `experiments/runs/spark_run1/` 目录下会生成：

```
├── config.toml                      # 配置文件
├── _config.json                     # JSON格式配置
├── training.log                     # 训练日志
├── results_1.csv                    # 所有生成分子的详细数据
├── gold_standard_candidates.csv     # 金标准候选分子
├── analysis_plots.png               # 分析图表
├── REPORT.txt                       # 最终报告
├── checkpoint_stage1.chkpt          # 训练检查点
└── tensorboard/                     # TensorBoard数据
    └── events.out.tfevents.*
```

## 故障排查

### 问题1: CUDA内存不足
```bash
# 解决方案: 降低batch_size
# 在config.toml中修改:
batch_size = 256  # 从512降低到256
```

### 问题2: QSAR插件加载失败
```bash
# 检查插件文件
ls -lh reinvent/reinvent_plugins/components/comp_qsar_scorer.py

# 检查模型文件
ls -lh models/random_forest_champion.joblib

# 重新安装插件
bash qsar_plugin_setup.sh
```

### 问题3: 训练速度慢
```bash
# 检查GPU利用率
nvidia-smi -l 1

# 如果GPU利用率低，可能是数据加载瓶颈
# 尝试增加数据预加载
```

### 问题4: 进程意外终止
```bash
# 查看日志最后的错误信息
tail -100 experiments/runs/spark_run1/training.log

# 如果是OOM错误，降低batch_size
# 如果是其他错误，检查配置文件语法
```

## 下一步工作

1. **结果分析**: 使用Jupyter notebook的分析cells
2. **候选筛选**: 从gold_standard_candidates.csv中选择候选
3. **虚拟筛选**: 使用LigUnity进行分子对接
4. **结构优化**: 对top候选进行进一步优化
5. **实验验证**: 选择10-50个候选进行合成和测试

## 性能基准测试结果（预期）

基于硬件规格的理论估算：

```
每步时间: ~6-10秒 (vs T1200的 ~40-60秒)
总训练时间: 8-12小时 (vs T1200的 24-36小时)
速度提升: 3-4倍
生成分子数: 2,560,000 (vs T1200的 ~1,000,000)
金标准候选: 3,000-6,000 (vs T1200的 100-300)
化学空间覆盖: ~90-95% (vs T1200的 ~30-50%)
```

## 联系信息

如有问题，请查看：
- REINVENT4文档: https://github.com/MolecularAI/REINVENT4
- 项目README: ~/aidd/REINVENT4/README.md
- 训练日志: experiments/runs/spark_run1/training.log

## 版本信息

- REINVENT4: 4.x
- PyTorch: 2.x (ARM + CUDA 13.0)
- RDKit: 2023.x
- Python: 3.10+
- CUDA: 13.0

---

**最后更新**: 2024
**作者**: Jason's AI Drug Discovery Pipeline
**硬件**: NVIDIA DGX Spark
