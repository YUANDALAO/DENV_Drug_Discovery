# 创建通用分析脚本
cat > analyze_results.py << 'EOF'
#!/usr/bin/env python3
"""
REINVENT结果通用分析脚本
用法: python analyze_results.py <run_folder>
示例: python analyze_results.py experiments/runs/run3
"""

import pandas as pd
import glob
import sys
import os
import matplotlib.pyplot as plt
import seaborn as sns

# 获取运行目录
if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"  # 默认

# 读取结果文件
files = glob.glob(f"{run_folder}/results_*.csv")
if not files:
    print(f"❌ 未找到结果文件: {run_folder}/results_*.csv")
    sys.exit(1)

latest_file = max(files, key=lambda x: x.split('_')[-1])
print(f"📂 读取文件: {latest_file}")

df = pd.read_csv(latest_file, low_memory=False)

print("=" * 80)
print(f"REINVENT {os.path.basename(run_folder)} 结果统计")
print("=" * 80)

# ============================================
# 1. 基本统计
# ============================================
print(f"\n总分子数: {len(df):,}")
print(f"唯一分子数: {df['SMILES'].nunique():,}")
print(f"去重率: {(1 - df['SMILES'].nunique()/len(df))*100:.1f}%")

# ============================================
# 2. 骨架匹配率
# ============================================
if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    pyrr_match = (df['Pyrrolidine_Scaffold (raw)'] > 0).sum()
    cyclo_match = (df['Cyclobutane_Scaffold (raw)'] > 0).sum()
    any_match = ((df['Pyrrolidine_Scaffold (raw)'] > 0) | 
                 (df['Cyclobutane_Scaffold (raw)'] > 0)).sum()
    
    print(f"\n骨架匹配统计:")
    print(f"  吡咯烷型: {pyrr_match:,} ({pyrr_match/len(df)*100:.1f}%)")
    print(f"  环丁烷型: {cyclo_match:,} ({cyclo_match/len(df)*100:.1f}%)")
    print(f"  任一骨架: {any_match:,} ({any_match/len(df)*100:.1f}%)")

# ============================================
# 3. Score分布
# ============================================
print(f"\nScore统计:")
print(f"  平均: {df['Score'].mean():.4f}")
print(f"  中位数: {df['Score'].median():.4f}")
print(f"  最大: {df['Score'].max():.4f}")
print(f"  >0.3: {(df['Score'] > 0.3).sum():,}")
print(f"  >0.5: {(df['Score'] > 0.5).sum():,}")
print(f"  >0.7: {(df['Score'] > 0.7).sum():,}")

# ============================================
# 4. pIC50统计
# ============================================
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    print(f"\npIC50统计:")
    print(f"  平均: {df['DENV_Activity_pIC50 (raw)'].mean():.2f}")
    print(f"  中位数: {df['DENV_Activity_pIC50 (raw)'].median():.2f}")
    print(f"  最大: {df['DENV_Activity_pIC50 (raw)'].max():.2f}")
    print(f"  >7.0: {(df['DENV_Activity_pIC50 (raw)'] > 7.0).sum():,}")
    print(f"  >8.0: {(df['DENV_Activity_pIC50 (raw)'] > 8.0).sum():,}")
    print(f"  >9.0: {(df['DENV_Activity_pIC50 (raw)'] > 9.0).sum():,}")

# ============================================
# 5. Top分子
# ============================================
print("\n" + "=" * 80)
print("Top 20 分子 (按Score排序)")
print("=" * 80)

cols = ['SMILES', 'Score', 'DENV_Activity_pIC50 (raw)']
if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    cols += ['Pyrrolidine_Scaffold (raw)', 'Cyclobutane_Scaffold (raw)']
if 'Drug_Likeness (raw)' in df.columns:
    cols += ['Drug_Likeness (raw)', 'Synthetic_Accessibility (raw)']

top20 = df.nlargest(20, 'Score')[cols].drop_duplicates(subset=['SMILES']).head(20)
print(top20.to_string(index=False))

# ============================================
# 6. 保存文件
# ============================================
# 保存匹配骨架的分子
if 'Pyrrolidine_Scaffold (raw)' in df.columns:
    matched = df[((df['Pyrrolidine_Scaffold (raw)'] > 0) | 
                  (df['Cyclobutane_Scaffold (raw)'] > 0))]
    if len(matched) > 0:
        output_file = f"{run_folder}/matched_scaffolds.csv"
        matched.to_csv(output_file, index=False)
        print(f"\n✅ 匹配骨架分子: {output_file} ({len(matched):,}个)")

# 保存Top分子
top_file = f"{run_folder}/top_molecules.csv"
top20.to_csv(top_file, index=False)
print(f"✅ Top20分子: {top_file}")

# 保存高活性分子
if 'DENV_Activity_pIC50 (raw)' in df.columns:
    high_activity = df[df['DENV_Activity_pIC50 (raw)'] > 8.0].sort_values(
        'DENV_Activity_pIC50 (raw)', ascending=False)
    if len(high_activity) > 0:
        high_file = f"{run_folder}/high_activity_pic50_gt8.csv"
        high_activity.to_csv(high_file, index=False)
        print(f"✅ 高活性分子(>8.0): {high_file} ({len(high_activity):,}个)")

print("\n" + "=" * 80)
print("分析完成！")
print("=" * 80)
EOF

chmod +x analyze_results.py