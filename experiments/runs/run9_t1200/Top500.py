import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

# 读取全量数据
df_all = pd.read_csv("full_evaluation/all_molecules_metrics.csv")

print(f"总分子数: {len(df_all)}")

# 第一轮：高质量
hq = df_all[
    (df_all['QED'] > 0.6) &
    (df_all['RO5_Violations'] <= 1) &
    (df_all['MW'].between(200, 600))
].copy()

print(f"高质量: {len(hq)}")

# 第二轮：严格筛选
best = hq[
    (hq['QED'] > 0.65) &           # 提高标准
    (hq['MW'] < 500) &             # 更严格
    (hq['MW'] > 350) &             # 不要太小
    (hq['LogP'] > 2.5) &           # 不要太亲水
    (hq['LogP'] < 4.5) &           # 不要太亲脂
    (hq['TPSA'] > 60) &            # 最小TPSA
    (hq['TPSA'] < 110) &           # 最大TPSA
    (hq['HBD'] <= 3) &             # 限制供体
    (hq['HBA'] <= 8) &             # 限制受体
    (hq['RotBonds'] <= 10)         # 限制柔性
].copy()

print(f"最佳候选: {len(best)}")

# 按QED排序
best = best.sort_values('QED', ascending=False)

# 保存Top 500
top500 = best.head(500)
top500.to_csv("top500_candidates.csv", index=False)

print(f"\n✓ Top 500已保存")
print("\nTop 10:")
print(top500[['SMILES', 'QED', 'MW', 'LogP', 'TPSA', 'Score']].head(10).to_string())

# 统计
print(f"\nTop 500统计:")
print(f"  QED: {top500['QED'].mean():.3f} (范围 {top500['QED'].min():.3f}-{top500['QED'].max():.3f})")
print(f"  MW:  {top500['MW'].mean():.1f} (范围 {top500['MW'].min():.0f}-{top500['MW'].max():.0f})")
print(f"  LogP: {top500['LogP'].mean():.2f} (范围 {top500['LogP'].min():.2f}-{top500['LogP'].max():.2f})")