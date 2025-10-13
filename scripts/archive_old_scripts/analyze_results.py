#!/usr/bin/env python3
import pandas as pd
import glob
import sys
import os

# 获取运行目录
if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"

# 读取结果文件
files = glob.glob(f"{run_folder}/results_*.csv")
if not files:
    print(f"未找到结果文件")
    sys.exit(1)

latest_file = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest_file, low_memory=False)

print("=" * 80)
print(f"{os.path.basename(run_folder)} 结果统计")
print("=" * 80)

print(f"\n总分子数: {len(df):,}")
print(f"唯一分子数: {df['SMILES'].nunique():,}")

if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    pyrr = (df['Pyrrolidine_Scaffold (raw)'] > 0).sum()
    cyclo = (df['Cyclobutane_Scaffold (raw)'] > 0).sum()
    any_match = ((df['Pyrrolidine_Scaffold (raw)'] > 0) | (df['Cyclobutane_Scaffold (raw)'] > 0)).sum()
    
    print(f"\n骨架匹配:")
    print(f"  吡咯烷: {pyrr:,} ({pyrr/len(df)*100:.1f}%)")
    print(f"  环丁烷: {cyclo:,} ({cyclo/len(df)*100:.1f}%)")
    print(f"  任一: {any_match:,} ({any_match/len(df)*100:.1f}%)")

print(f"\nScore: 平均={df['Score'].mean():.4f}, 最大={df['Score'].max():.4f}")

if 'DENV_Activity_pIC50 (raw)' in df.columns:
    print(f"pIC50: 平均={df['DENV_Activity_pIC50 (raw)'].mean():.2f}, 最大={df['DENV_Activity_pIC50 (raw)'].max():.2f}")
    print(f"  >8.0: {(df['DENV_Activity_pIC50 (raw)'] > 8.0).sum():,}")

print("\nTop 10:")
top = df.nlargest(10, 'Score')[['SMILES', 'Score', 'DENV_Activity_pIC50 (raw)']].drop_duplicates('SMILES')
print(top.to_string(index=False))

# 保存
matched = df[((df['Pyrrolidine_Scaffold (raw)'] > 0) | (df['Cyclobutane_Scaffold (raw)'] > 0))]
matched.to_csv(f"{run_folder}/matched.csv", index=False)
top.to_csv(f"{run_folder}/top10.csv", index=False)
print(f"\n保存: {run_folder}/matched.csv, {run_folder}/top10.csv")
