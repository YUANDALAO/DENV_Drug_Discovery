#!/usr/bin/env python3
"""
REINVENT完整分析可视化 - 9子图版本
用法: python visualize_full_analysis.py experiments/runs/run3
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import sys
import os

if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"

files = glob.glob(f"{run_folder}/results_*.csv")
latest_file = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest_file, low_memory=False)

print(f"生成9子图分析: {run_folder}")

# 转换数据类型并计算IC50
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    df['DENV_Activity_pIC50 (raw)'] = pd.to_numeric(df['DENV_Activity_pIC50 (raw)'], errors='coerce')
    df['IC50_nM'] = 10**(9 - df['DENV_Activity_pIC50 (raw)'])

if 'Drug_Likeness (raw)' in df.columns:
    df['Drug_Likeness (raw)'] = pd.to_numeric(df['Drug_Likeness (raw)'], errors='coerce')

if 'Synthetic_Accessibility (raw)' in df.columns:
    df['Synthetic_Accessibility (raw)'] = pd.to_numeric(df['Synthetic_Accessibility (raw)'], errors='coerce')

if 'Molecular_Weight (raw)' in df.columns:
    df['Molecular_Weight (raw)'] = pd.to_numeric(df['Molecular_Weight (raw)'], errors='coerce')

if 'LogP (raw)' in df.columns:
    df['LogP (raw)'] = pd.to_numeric(df['LogP (raw)'], errors='coerce')

# 创建图表
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.suptitle('DENV Inhibitor Generation Analysis', fontsize=16, fontweight='bold')

# 1. Learning Curve
if 'step' in df.columns:
    step_mean = df.groupby('step')['Score'].mean()
    axes[0,0].plot(step_mean.index, step_mean.values, linewidth=2)
    axes[0,0].fill_between(step_mean.index, step_mean.values, alpha=0.3)
axes[0,0].set_xlabel('step')
axes[0,0].set_ylabel('Mean Score')
axes[0,0].set_title('Learning Curve')
axes[0,0].grid(True, alpha=0.3)

# 2. Activity Distribution (pIC50)
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    axes[0,1].hist(df['DENV_Activity_pIC50 (raw)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[0,1].axvline(8.0, color='red', linestyle='--', linewidth=2, label='Target (IC50=32nM)')
    axes[0,1].legend()
axes[0,1].set_xlabel('pIC50')
axes[0,1].set_ylabel('Frequency')
axes[0,1].set_title('Activity Distribution')

# 3. IC50 Distribution
if 'IC50_nM' in df.columns:
    axes[0,2].hist(np.log10(df['IC50_nM'].dropna()), bins=50, edgecolor='black', alpha=0.7)
axes[0,2].set_xlabel('log10(IC50 in nM)')
axes[0,2].set_ylabel('Frequency')
axes[0,2].set_title('IC50 Distribution')

# 4. Drug-likeness (QED)
if 'Drug_Likeness (raw)' in df.columns:
    axes[1,0].hist(df['Drug_Likeness (raw)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[1,0].axvline(0.7, color='red', linestyle='--', linewidth=2, label='Threshold')
    axes[1,0].legend()
axes[1,0].set_xlabel('QED')
axes[1,0].set_ylabel('Frequency')
axes[1,0].set_title('Drug-likeness')

# 5. Synthetic Accessibility
if 'Synthetic_Accessibility (raw)' in df.columns:
    axes[1,1].hist(df['Synthetic_Accessibility (raw)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
axes[1,1].set_xlabel('SA Score (1=easy, 10=hard)')
axes[1,1].set_ylabel('Frequency')
axes[1,1].set_title('Synthetic Accessibility')

# 6. Molecular Weight Distribution
if 'Molecular_Weight (raw)' in df.columns:
    axes[1,2].hist(df['Molecular_Weight (raw)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[1,2].axvspan(250, 500, alpha=0.2, color='green', label='Ideal range')
    axes[1,2].legend()
axes[1,2].set_xlabel('Molecular Weight (Da)')
axes[1,2].set_ylabel('Frequency')
axes[1,2].set_title('Molecular Weight Distribution')

# 7. Lipophilicity (LogP)
if 'LogP (raw)' in df.columns:
    axes[2,0].hist(df['LogP (raw)'].dropna(), bins=50, edgecolor='black', alpha=0.7)
    axes[2,0].axvspan(1, 5, alpha=0.2, color='green', label='Ideal range')
    axes[2,0].legend()
axes[2,0].set_xlabel('LogP')
axes[2,0].set_ylabel('Frequency')
axes[2,0].set_title('Lipophilicity')

# 8. Activity vs Drug-likeness
if 'DENV_Activity_pIC50 (raw)' in df.columns and 'Drug_Likeness (raw)' in df.columns:
    valid_data = df[['DENV_Activity_pIC50 (raw)', 'Drug_Likeness (raw)', 'Score']].dropna()
    scatter = axes[2,1].scatter(valid_data['DENV_Activity_pIC50 (raw)'], 
                               valid_data['Drug_Likeness (raw)'],
                               c=valid_data['Score'], cmap='viridis', alpha=0.5, s=10)
    plt.colorbar(scatter, ax=axes[2,1], label='Total Score')
axes[2,1].set_xlabel('pIC50 (Activity)')
axes[2,1].set_ylabel('QED (Drug-likeness)')
axes[2,1].set_title('Activity vs Drug-likeness')

# 9. Activity vs Synthesizability
if 'DENV_Activity_pIC50 (raw)' in df.columns and 'Synthetic_Accessibility (raw)' in df.columns:
    valid_data = df[['DENV_Activity_pIC50 (raw)', 'Synthetic_Accessibility (raw)', 'Score']].dropna()
    scatter = axes[2,2].scatter(valid_data['DENV_Activity_pIC50 (raw)'], 
                               valid_data['Synthetic_Accessibility (raw)'],
                               c=valid_data['Score'], cmap='RdYlGn_r', alpha=0.5, s=10)
    plt.colorbar(scatter, ax=axes[2,2], label='Total Score')
axes[2,2].set_xlabel('pIC50 (higher=better)')
axes[2,2].set_ylabel('SA Score (lower=easier)')
axes[2,2].set_title('Activity vs Synthesizability')

plt.tight_layout()
output_file = f'{run_folder}/generation_analysis.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"9子图分析: {output_file}")
plt.close()

# 输出统计信息
print("\n=== 基本统计 ===")
print(f"总分子数: {len(df):,}")
print(f"唯一分子: {df['SMILES'].nunique():,}")
if 'IC50_nM' in df.columns:
    print(f"\n平均IC50: {df['IC50_nM'].mean():.2f} nM (被极端值拉高)")
    print(f"中位IC50: {df['IC50_nM'].median():.2f} nM (更准确)")
    print(f"最佳IC50: {df['IC50_nM'].min():.2f} nM")

# 筛选优质候选物 - 添加类型转换
try:
    good = df[
        (df['DENV_Activity_pIC50 (raw)'] > 8.0) &
        (df['Drug_Likeness (raw)'] > 0.7) &
        (df['Synthetic_Accessibility (raw)'] < 4.0) &
        (df['Molecular_Weight (raw)'] >= 250) &
        (df['Molecular_Weight (raw)'] <= 500) &
        (df['LogP (raw)'] >= 1) &
        (df['LogP (raw)'] <= 5)
    ]
    
    print(f"\n满足所有条件的分子: {len(good):,} ({len(good)/len(df)*100:.1f}%)")
    print("条件: pIC50>8.0, QED>0.7, SA<4.0, MW 250-500, LogP 1-5")
    
    if len(good) > 0:
        good_file = f'{run_folder}/promising_candidates.csv'
        good.to_csv(good_file, index=False)
        print(f"已保存: {good_file}")
except Exception as e:
    print(f"\n筛选优质候选物时出错: {e}")
