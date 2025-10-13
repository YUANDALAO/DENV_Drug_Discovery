#!/usr/bin/env python3
"""
可视化Top结构 - 修复版
用法: python visualize_top_structures_fixed.py experiments/runs/run7_optimized
"""
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import sys
import os
import glob

if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run7_optimized"

# 1. 首先尝试读取 promising_candidates.csv（如果存在）
promising_file = f"{run_folder}/promising_candidates.csv"
results_files = glob.glob(f"{run_folder}/results_*.csv")

if os.path.exists(promising_file):
    print(f"✓ 找到优质候选物文件: {promising_file}")
    df = pd.read_csv(promising_file)
    source = "promising_candidates"
elif results_files:
    latest_file = max(results_files, key=lambda x: x.split('_')[-1])
    print(f"使用完整结果文件: {latest_file}")
    df = pd.read_csv(latest_file, low_memory=False)
    source = "results"
else:
    print(f"❌ 错误: 在 {run_folder} 中找不到数据文件")
    sys.exit(1)

print(f"数据形状: {df.shape}")
print(f"可用列: {df.columns.tolist()[:10]}...")  # 只显示前10列

# 动态查找活性列
def find_column(patterns):
    for pattern in patterns:
        matches = [col for col in df.columns if pattern.lower() in col.lower()]
        if matches:
            return matches[0]
    return None

activity_col = find_column(['DENV_Activity (raw)', 'Activity (raw)', 'pIC50 (raw)', 'DENV_Activity'])
qed_col = find_column(['QED (raw)', 'Drug_Likeness (raw)', 'QED'])
sa_col = find_column(['SA (raw)', 'Synthetic_Accessibility (raw)', 'SA'])

print(f"\n检测到的列:")
print(f"  Activity: {activity_col}")
print(f"  QED: {qed_col}")
print(f"  SA: {sa_col}")

# 去重并排序
df_unique = df.drop_duplicates(subset=['SMILES'])
print(f"\n去重后: {len(df_unique):,} 个唯一分子")

# 选择Top分子数量
n_top = min(20, len(df_unique))
top_molecules = df_unique.nlargest(n_top, 'Score')

print(f"选择Top {n_top}分子进行可视化")

# 准备分子和图例
mols = []
legends = []

for idx, row in top_molecules.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        
        # 构建图例
        legend_parts = [f"Score: {row['Score']:.4f}"]
        
        if activity_col and activity_col in row:
            activity_val = row[activity_col]
            if pd.notna(activity_val):
                legend_parts.append(f"pIC50: {activity_val:.2f}")
                # 计算IC50
                ic50 = 10**(9 - activity_val)
                if ic50 < 1000:
                    legend_parts.append(f"IC50: {ic50:.1f} nM")
                else:
                    legend_parts.append(f"IC50: {ic50/1000:.2f} μM")
        
        if qed_col and qed_col in row:
            qed_val = row[qed_col]
            if pd.notna(qed_val):
                legend_parts.append(f"QED: {qed_val:.3f}")
        
        if sa_col and sa_col in row:
            sa_val = row[sa_col]
            if pd.notna(sa_val):
                legend_parts.append(f"SA: {sa_val:.2f}")
        
        legend = "\n".join(legend_parts)
        legends.append(legend)

if not mols:
    print("❌ 没有有效的分子可以可视化")
    sys.exit(1)

print(f"\n✓ 成功解析 {len(mols)} 个分子")

# 生成网格图像
output_file = f"{run_folder}/top{n_top}_structures.png"

try:
    img = Draw.MolsToGridImage(
        mols, 
        molsPerRow=4, 
        subImgSize=(400, 400), 
        legends=legends,
        returnPNG=False
    )
    img.save(output_file)
    print(f"✓ 已保存结构图: {output_file}")
    
    # 同时保存分子详细信息
    info_file = f"{run_folder}/top{n_top}_info.csv"
    
    # 选择要保存的列
    cols_to_save = ['SMILES', 'Score']
    if activity_col:
        cols_to_save.append(activity_col)
    if qed_col:
        cols_to_save.append(qed_col)
    if sa_col:
        cols_to_save.append(sa_col)
    
    # 添加其他可能有用的列
    other_cols = ['MW', 'LogP', 'HBA', 'HBD', 'TPSA', 'NumRotBond']
    for col in other_cols:
        matched = find_column([col])
        if matched and matched not in cols_to_save:
            cols_to_save.append(matched)
    
    available_cols = [c for c in cols_to_save if c in top_molecules.columns]
    top_molecules[available_cols].to_csv(info_file, index=False)
    print(f"✓ 已保存分子信息: {info_file}")
    
    # 显示统计摘要
    print("\n" + "="*60)
    print(f"Top {n_top} 分子统计摘要".center(60))
    print("="*60)
    
    if activity_col and activity_col in top_molecules.columns:
        valid_activity = top_molecules[activity_col].dropna()
        if len(valid_activity) > 0:
            print(f"\npIC50:")
            print(f"  范围: {valid_activity.min():.2f} - {valid_activity.max():.2f}")
            print(f"  平均: {valid_activity.mean():.2f}")
            
            # IC50范围
            ic50s = 10**(9 - valid_activity)
            print(f"\nIC50:")
            print(f"  最佳: {ic50s.min():.2f} nM")
            print(f"  最差: {ic50s.max():.2f} nM")
    
    if qed_col and qed_col in top_molecules.columns:
        valid_qed = top_molecules[qed_col].dropna()
        if len(valid_qed) > 0:
            print(f"\nQED: {valid_qed.min():.3f} - {valid_qed.max():.3f} (平均: {valid_qed.mean():.3f})")
    
    if sa_col and sa_col in top_molecules.columns:
        valid_sa = top_molecules[sa_col].dropna()
        if len(valid_sa) > 0:
            print(f"SA:  {valid_sa.min():.2f} - {valid_sa.max():.2f} (平均: {valid_sa.mean():.2f})")
    
    print("\n" + "="*60)
    
except Exception as e:
    print(f"❌ 生成图像时出错: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)