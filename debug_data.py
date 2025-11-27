import os
import pandas as pd

# 检查一个run的所有文件
run_dir = "experiments/runs/run9_t1200"  # 改成你的run名称

print(f"检查目录: {run_dir}\n")

print("=== 所有CSV文件 ===")
for f in os.listdir(run_dir):
    if f.endswith('.csv'):
        path = os.path.join(run_dir, f)
        df = pd.read_csv(path)
        print(f"\n{f}:")
        print(f"  行数: {len(df)}")
        print(f"  列: {df.columns.tolist()}")
        print(f"  前3行:")
        print(df.head(3))

print("\n=== config.toml内容 ===")
config_path = os.path.join(run_dir, "config.toml")
if os.path.exists(config_path):
    with open(config_path, 'r') as f:
        print(f.read()[:500])
