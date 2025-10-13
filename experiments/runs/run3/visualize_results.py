cat > visualize_results.py << 'EOF'
#!/usr/bin/env python3
"""
REINVENT结果可视化脚本
用法: python visualize_results.py <run_folder>
"""

import pandas as pd
import glob
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# 设置中文字体（可选）
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("whitegrid")

# 获取运行目录
if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"

# 读取数据
files = glob.glob(f"{run_folder}/results_*.csv")
latest_file = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest_file, low_memory=False)

print(f"📊 生成可视化图表: {run_folder}")

# 创建图表目录
plot_dir = f"{run_folder}/plots"
os.makedirs(plot_dir, exist_ok=True)

# ============================================
# 图1: Score分布
# ============================================
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['Score'], bins=50, edgecolor='black', alpha=0.7)
ax.axvline(df['Score'].mean(), color='red', linestyle='--', 
           label=f'Mean: {df["Score"].mean():.3f}')
ax.set_xlabel('Score', fontsize=12)
ax.set_ylabel('Frequency', fontsize=12)
ax.set_title('Score Distribution', fontsize=14, fontweight='bold')
ax.legend()
plt.tight_layout()
plt.savefig(f'{plot_dir}/score_distribution.png', dpi=300)
print(f"  ✅ {plot_dir}/score_distribution.png")
plt.close()

# ============================================
# 图2: pIC50分布
# ============================================
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(df['DENV_Activity_pIC50 (raw)'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(df['DENV_Activity_pIC50 (raw)'].mean(), color='red', 
               linestyle='--', label=f'Mean: {df["DENV_Activity_pIC50 (raw)"].mean():.2f}')
    ax.axvline(8.0, color='green', linestyle='--', label='Target: 8.0')
    ax.set_xlabel('pIC50', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title('pIC50 Distribution', fontsize=14, fontweight='bold')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/pic50_distribution.png', dpi=300)
    print(f"  ✅ {plot_dir}/pic50_distribution.png")
    plt.close()

# ============================================
# 图3: Score vs pIC50 散点图
# ============================================
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(df['DENV_Activity_pIC50 (raw)'], df['Score'], 
                        alpha=0.3, s=10)
    ax.set_xlabel('pIC50', fontsize=12)
    ax.set_ylabel('Score', fontsize=12)
    ax.set_title('Score vs pIC50', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/score_vs_pic50.png', dpi=300)
    print(f"  ✅ {plot_dir}/score_vs_pic50.png")
    plt.close()

# ============================================
# 图4: 骨架匹配饼图
# ============================================
if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    pyrr_only = ((df['Pyrrolidine_Scaffold (raw)'] > 0) & 
                 (df['Cyclobutane_Scaffold (raw)'] == 0)).sum()
    cyclo_only = ((df['Cyclobutane_Scaffold (raw)'] > 0) & 
                  (df['Pyrrolidine_Scaffold (raw)'] == 0)).sum()
    both = ((df['Pyrrolidine_Scaffold (raw)'] > 0) & 
            (df['Cyclobutane_Scaffold (raw)'] > 0)).sum()
    neither = ((df['Pyrrolidine_Scaffold (raw)'] == 0) & 
               (df['Cyclobutane_Scaffold (raw)'] == 0)).sum()
    
    fig, ax = plt.subplots(figsize=(8, 8))
    sizes = [pyrr_only, cyclo_only, both, neither]
    labels = [f'Pyrrolidine only\n({pyrr_only:,})', 
              f'Cyclobutane only\n({cyclo_only:,})',
              f'Both\n({both:,})', 
              f'Neither\n({neither:,})']
    colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99']
    
    ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
    ax.set_title('Scaffold Matching Distribution', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/scaffold_pie.png', dpi=300)
    print(f"  ✅ {plot_dir}/scaffold_pie.png")
    plt.close()

# ============================================
# 图5: 分子量vs LogP
# ============================================
if 'Molecular_Weight (raw)' in df.columns and 'LogP (raw)' in df.columns:
    fig, ax = plt.subplots(figsize=(10, 8))
    scatter = ax.scatter(df['Molecular_Weight (raw)'], df['LogP (raw)'], 
                        c=df['Score'], cmap='viridis', alpha=0.5, s=10)
    ax.set_xlabel('Molecular Weight', fontsize=12)
    ax.set_ylabel('LogP', fontsize=12)
    ax.set_title('Molecular Weight vs LogP (colored by Score)', 
                 fontsize=14, fontweight='bold')
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Score', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/mw_vs_logp.png', dpi=300)
    print(f"  ✅ {plot_dir}/mw_vs_logp.png")
    plt.close()

# ============================================
# 图6: 训练进度（如果有Step列）
# ============================================
if 'Step' in df.columns:
    step_stats = df.groupby('Step').agg({
        'Score': 'mean',
        'DENV_Activity_pIC50 (raw)': 'mean'
    }).reset_index()
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    ax1.plot(step_stats['Step'], step_stats['Score'], linewidth=2)
    ax1.set_xlabel('Step', fontsize=12)
    ax1.set_ylabel('Average Score', fontsize=12)
    ax1.set_title('Training Progress: Score', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(step_stats['Step'], step_stats['DENV_Activity_pIC50 (raw)'], 
             linewidth=2, color='orange')
    ax2.set_xlabel('Step', fontsize=12)
    ax2.set_ylabel('Average pIC50', fontsize=12)
    ax2.set_title('Training Progress: pIC50', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/training_progress.png', dpi=300)
    print(f"  ✅ {plot_dir}/training_progress.png")
    plt.close()

print(f"\n✅ 所有图表已保存到: {plot_dir}/")
EOF

chmod +x visualize_results.py