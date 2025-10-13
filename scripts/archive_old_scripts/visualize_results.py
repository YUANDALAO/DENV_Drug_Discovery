#!/usr/bin/env python3
import pandas as pd
import glob
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.sans-serif'] = ['Arial']
sns.set_style("whitegrid")

if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"

files = glob.glob(f"{run_folder}/results_*.csv")
latest_file = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest_file, low_memory=False)

print(f"生成可视化: {run_folder}")

plot_dir = f"{run_folder}/plots"
os.makedirs(plot_dir, exist_ok=True)

# 图1: Score分布
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['Score'], bins=50, edgecolor='black', alpha=0.7)
ax.axvline(df['Score'].mean(), color='red', linestyle='--', label=f'Mean: {df["Score"].mean():.3f}')
ax.set_xlabel('Score')
ax.set_ylabel('Frequency')
ax.set_title('Score Distribution')
ax.legend()
plt.tight_layout()
plt.savefig(f'{plot_dir}/score_distribution.png', dpi=300)
print(f"  {plot_dir}/score_distribution.png")
plt.close()

# 图2: pIC50分布
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(df['DENV_Activity_pIC50 (raw)'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(df['DENV_Activity_pIC50 (raw)'].mean(), color='red', linestyle='--')
    ax.axvline(8.0, color='green', linestyle='--', label='Target: 8.0')
    ax.set_xlabel('pIC50')
    ax.set_ylabel('Frequency')
    ax.set_title('pIC50 Distribution')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/pic50_distribution.png', dpi=300)
    print(f"  {plot_dir}/pic50_distribution.png")
    plt.close()

# 图3: Score vs pIC50
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(df['DENV_Activity_pIC50 (raw)'], df['Score'], alpha=0.3, s=10)
    ax.set_xlabel('pIC50')
    ax.set_ylabel('Score')
    ax.set_title('Score vs pIC50')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/score_vs_pic50.png', dpi=300)
    print(f"  {plot_dir}/score_vs_pic50.png")
    plt.close()

# 图4: 骨架匹配饼图
if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    pyrr_only = ((df['Pyrrolidine_Scaffold (raw)'] > 0) & (df['Cyclobutane_Scaffold (raw)'] == 0)).sum()
    cyclo_only = ((df['Cyclobutane_Scaffold (raw)'] > 0) & (df['Pyrrolidine_Scaffold (raw)'] == 0)).sum()
    both = ((df['Pyrrolidine_Scaffold (raw)'] > 0) & (df['Cyclobutane_Scaffold (raw)'] > 0)).sum()
    neither = ((df['Pyrrolidine_Scaffold (raw)'] == 0) & (df['Cyclobutane_Scaffold (raw)'] == 0)).sum()
    
    fig, ax = plt.subplots(figsize=(8, 8))
    sizes = [pyrr_only, cyclo_only, both, neither]
    labels = [f'Pyrr ({pyrr_only:,})', f'Cyclo ({cyclo_only:,})', f'Both ({both:,})', f'None ({neither:,})']
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax.set_title('Scaffold Matching')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/scaffold_pie.png', dpi=300)
    print(f"  {plot_dir}/scaffold_pie.png")
    plt.close()

print(f"\n所有图表: {plot_dir}/")
