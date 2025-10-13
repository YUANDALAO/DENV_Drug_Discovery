#!/usr/bin/env python3
"""
REINVENT完整分析可视化 - 修复版
用法: python visualize_full_analysis_fixed.py experiments/runs/run7_optimized
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
    run_folder = "experiments/runs/run7_optimized"

files = glob.glob(f"{run_folder}/results_*.csv")
if not files:
    print(f"错误: 在 {run_folder} 中找不到 results_*.csv 文件")
    sys.exit(1)

latest_file = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest_file, low_memory=False)

print(f"数据文件: {latest_file}")
print(f"数据形状: {df.shape}")
print(f"\n可用列名:")
for i, col in enumerate(df.columns, 1):
    print(f"  {i}. {col}")

# 动态检测列名（适配不同的命名方式）
def find_column(patterns):
    """根据模式列表查找列名"""
    for pattern in patterns:
        matches = [col for col in df.columns if pattern.lower() in col.lower()]
        if matches:
            return matches[0]
    return None

# 查找关键列 - 优先使用 (raw) 列获取真实数值
activity_col = find_column(['DENV_Activity (raw)', 'Activity (raw)', 'pIC50 (raw)', 'DENV_Activity'])
qed_col = find_column(['QED (raw)', 'Drug_Likeness (raw)', 'QED'])
sa_col = find_column(['SA (raw)', 'Synthetic_Accessibility (raw)', 'SA'])
mw_col = find_column(['MW (raw)', 'Molecular_Weight (raw)', 'MolecularWeight (raw)', 'MW'])
logp_col = find_column(['LogP (raw)', 'SlogP (raw)', 'LogP'])
hba_col = find_column(['HBA (raw)', 'HBondAcceptors (raw)', 'HBA'])
hbd_col = find_column(['HBD (raw)', 'HBondDonors (raw)', 'HBD'])
tpsa_col = find_column(['TPSA (raw)', 'TPSA'])

print(f"\n检测到的关键列:")
print(f"  Activity: {activity_col}")
print(f"  QED: {qed_col}")
print(f"  SA: {sa_col}")
print(f"  MW: {mw_col}")
print(f"  LogP: {logp_col}")

# 转换数据类型
if activity_col:
    df[activity_col] = pd.to_numeric(df[activity_col], errors='coerce')
    df['IC50_nM'] = 10**(9 - df[activity_col])

if qed_col:
    df[qed_col] = pd.to_numeric(df[qed_col], errors='coerce')

if sa_col:
    df[sa_col] = pd.to_numeric(df[sa_col], errors='coerce')

if mw_col:
    df[mw_col] = pd.to_numeric(df[mw_col], errors='coerce')

if logp_col:
    df[logp_col] = pd.to_numeric(df[logp_col], errors='coerce')

# 创建图表
fig, axes = plt.subplots(3, 3, figsize=(18, 15))
fig.suptitle(f'DENV Inhibitor Generation Analysis - {os.path.basename(run_folder)}', fontsize=16, fontweight='bold')

# 1. Learning Curve
if 'step' in df.columns:
    step_mean = df.groupby('step')['Score'].mean()
    axes[0,0].plot(step_mean.index, step_mean.values, linewidth=2, color='#2E86AB')
    axes[0,0].fill_between(step_mean.index, step_mean.values, alpha=0.3, color='#2E86AB')
axes[0,0].set_xlabel('Step', fontsize=10)
axes[0,0].set_ylabel('Mean Score', fontsize=10)
axes[0,0].set_title('Learning Curve', fontsize=12, fontweight='bold')
axes[0,0].grid(True, alpha=0.3)

# 2. Activity Distribution (pIC50)
if activity_col:
    valid_activity = df[activity_col].dropna()
    axes[0,1].hist(valid_activity, bins=50, edgecolor='black', alpha=0.7, color='#A23B72')
    axes[0,1].axvline(8.0, color='red', linestyle='--', linewidth=2, label='Target (IC50=10nM)')
    axes[0,1].legend()
    axes[0,1].set_xlabel('pIC50', fontsize=10)
    axes[0,1].set_ylabel('Frequency', fontsize=10)
    axes[0,1].set_title('Activity Distribution', fontsize=12, fontweight='bold')
else:
    axes[0,1].text(0.5, 0.5, 'Activity data not found', ha='center', va='center')
    axes[0,1].set_title('Activity Distribution', fontsize=12)

# 3. IC50 Distribution
if 'IC50_nM' in df.columns:
    valid_ic50 = df['IC50_nM'].dropna()
    # 过滤极端值
    valid_ic50_filtered = valid_ic50[valid_ic50 < 1e6]
    axes[0,2].hist(np.log10(valid_ic50_filtered), bins=50, edgecolor='black', alpha=0.7, color='#F18F01')
    axes[0,2].axvline(np.log10(10), color='red', linestyle='--', linewidth=2, label='IC50=10nM')
    axes[0,2].legend()
    axes[0,2].set_xlabel('log10(IC50 in nM)', fontsize=10)
    axes[0,2].set_ylabel('Frequency', fontsize=10)
    axes[0,2].set_title('IC50 Distribution', fontsize=12, fontweight='bold')
else:
    axes[0,2].text(0.5, 0.5, 'IC50 data not found', ha='center', va='center')

# 4. Drug-likeness (QED)
if qed_col:
    valid_qed = df[qed_col].dropna()
    axes[1,0].hist(valid_qed, bins=50, edgecolor='black', alpha=0.7, color='#6A994E')
    axes[1,0].axvline(0.7, color='red', linestyle='--', linewidth=2, label='Threshold')
    axes[1,0].legend()
    axes[1,0].set_xlabel('QED', fontsize=10)
    axes[1,0].set_ylabel('Frequency', fontsize=10)
    axes[1,0].set_title('Drug-likeness', fontsize=12, fontweight='bold')
else:
    axes[1,0].text(0.5, 0.5, 'QED data not found', ha='center', va='center')

# 5. Synthetic Accessibility
if sa_col:
    valid_sa = df[sa_col].dropna()
    axes[1,1].hist(valid_sa, bins=50, edgecolor='black', alpha=0.7, color='#BC4749')
    axes[1,1].set_xlabel('SA Score (1=easy, 10=hard)', fontsize=10)
    axes[1,1].set_ylabel('Frequency', fontsize=10)
    axes[1,1].set_title('Synthetic Accessibility', fontsize=12, fontweight='bold')
else:
    axes[1,1].text(0.5, 0.5, 'SA data not found', ha='center', va='center')

# 6. Molecular Weight Distribution
if mw_col:
    valid_mw = df[mw_col].dropna()
    axes[1,2].hist(valid_mw, bins=50, edgecolor='black', alpha=0.7, color='#386641')
    axes[1,2].axvspan(300, 500, alpha=0.2, color='green', label='Ideal range')
    axes[1,2].legend()
    axes[1,2].set_xlabel('Molecular Weight (Da)', fontsize=10)
    axes[1,2].set_ylabel('Frequency', fontsize=10)
    axes[1,2].set_title('Molecular Weight Distribution', fontsize=12, fontweight='bold')
else:
    axes[1,2].text(0.5, 0.5, 'MW data not found', ha='center', va='center')

# 7. Lipophilicity (LogP)
if logp_col:
    valid_logp = df[logp_col].dropna()
    axes[2,0].hist(valid_logp, bins=50, edgecolor='black', alpha=0.7, color='#9D4EDD')
    axes[2,0].axvspan(1, 4, alpha=0.2, color='green', label='Ideal range')
    axes[2,0].legend()
    axes[2,0].set_xlabel('LogP', fontsize=10)
    axes[2,0].set_ylabel('Frequency', fontsize=10)
    axes[2,0].set_title('Lipophilicity', fontsize=12, fontweight='bold')
else:
    axes[2,0].text(0.5, 0.5, 'LogP data not found', ha='center', va='center')

# 8. Activity vs Drug-likeness
if activity_col and qed_col:
    valid_data = df[[activity_col, qed_col, 'Score']].dropna()
    if len(valid_data) > 0:
        scatter = axes[2,1].scatter(valid_data[activity_col], 
                                   valid_data[qed_col],
                                   c=valid_data['Score'], cmap='viridis', alpha=0.5, s=10)
        plt.colorbar(scatter, ax=axes[2,1], label='Total Score')
        axes[2,1].set_xlabel('pIC50 (Activity)', fontsize=10)
        axes[2,1].set_ylabel('QED (Drug-likeness)', fontsize=10)
        axes[2,1].set_title('Activity vs Drug-likeness', fontsize=12, fontweight='bold')
else:
    axes[2,1].text(0.5, 0.5, 'Insufficient data', ha='center', va='center')

# 9. Activity vs Synthesizability
if activity_col and sa_col:
    valid_data = df[[activity_col, sa_col, 'Score']].dropna()
    if len(valid_data) > 0:
        scatter = axes[2,2].scatter(valid_data[activity_col], 
                                   valid_data[sa_col],
                                   c=valid_data['Score'], cmap='RdYlGn_r', alpha=0.5, s=10)
        plt.colorbar(scatter, ax=axes[2,2], label='Total Score')
        axes[2,2].set_xlabel('pIC50 (higher=better)', fontsize=10)
        axes[2,2].set_ylabel('SA Score (lower=easier)', fontsize=10)
        axes[2,2].set_title('Activity vs Synthesizability', fontsize=12, fontweight='bold')
else:
    axes[2,2].text(0.5, 0.5, 'Insufficient data', ha='center', va='center')

plt.tight_layout()
output_file = f'{run_folder}/generation_analysis.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ 已保存9子图分析: {output_file}")
plt.close()

# 输出统计信息
print("\n" + "="*60)
print("基本统计".center(60))
print("="*60)
print(f"总分子数: {len(df):,}")
print(f"唯一分子: {df['SMILES'].nunique():,}")

if 'IC50_nM' in df.columns:
    valid_ic50 = df['IC50_nM'].dropna()
    print(f"\nIC50 统计:")
    print(f"  平均IC50: {valid_ic50.mean():.2f} nM")
    print(f"  中位IC50: {valid_ic50.median():.2f} nM")
    print(f"  最佳IC50: {valid_ic50.min():.2f} nM")
    print(f"  最差IC50: {valid_ic50.max():.2f} nM")

# 筛选优质候选物
try:
    conditions = []
    filter_desc = []
    
    if activity_col:
        conditions.append(df[activity_col] > 8.0)
        filter_desc.append("pIC50>8.0")
    
    if qed_col:
        conditions.append(df[qed_col] > 0.7)
        filter_desc.append("QED>0.7")
    
    if sa_col:
        conditions.append(df[sa_col] < 4.0)
        filter_desc.append("SA<4.0")
    
    if mw_col:
        conditions.append((df[mw_col] >= 300) & (df[mw_col] <= 500))
        filter_desc.append("MW 300-500")
    
    if logp_col:
        conditions.append((df[logp_col] >= 1) & (df[logp_col] <= 4))
        filter_desc.append("LogP 1-4")
    
    if conditions:
        # 组合所有条件
        combined_condition = conditions[0]
        for cond in conditions[1:]:
            combined_condition = combined_condition & cond
        
        good = df[combined_condition].copy()
        
        print("\n" + "="*60)
        print("优质候选物筛选".center(60))
        print("="*60)
        print(f"筛选条件: {', '.join(filter_desc)}")
        print(f"满足所有条件: {len(good):,} ({len(good)/len(df)*100:.2f}%)")
        
        if len(good) > 0:
            # 按Score排序
            good = good.sort_values('Score', ascending=False)
            
            # 保存结果
            good_file = f'{run_folder}/promising_candidates.csv'
            good.to_csv(good_file, index=False)
            print(f"✓ 已保存优质候选物: {good_file}")
            
            # 显示Top 5
            print("\nTop 5 候选物:")
            display_cols = ['SMILES', 'Score']
            if activity_col:
                display_cols.append(activity_col)
            if qed_col:
                display_cols.append(qed_col)
            if sa_col:
                display_cols.append(sa_col)
            
            available_cols = [c for c in display_cols if c in good.columns]
            print(good[available_cols].head().to_string(index=False))
        else:
            print("⚠ 没有分子满足所有条件，尝试放宽筛选标准")
    else:
        print("\n⚠ 缺少足够的列进行筛选")
        
except Exception as e:
    print(f"\n❌ 筛选优质候选物时出错: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*60)