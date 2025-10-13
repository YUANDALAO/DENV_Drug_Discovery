import os
import glob
from datetime import datetime

# 查找所有可能的模型文件
print("=== 查找模型文件 ===")
for pattern in ["*.chkpt", "*.model", "*.prior", "*.ckpt"]:
    files = glob.glob(pattern)
    for f in files:
        mtime = os.path.getmtime(f)
        size = os.path.getsize(f)
        print(f"{f}: {datetime.fromtimestamp(mtime)}, {size/1024/1024:.2f} MB")

# 检查CSV输出
print("\n=== 检查输出文件 ===")
csv_files = glob.glob("*.csv")
for csv in csv_files:
    lines = sum(1 for _ in open(csv))
    print(f"{csv}: {lines} 行")