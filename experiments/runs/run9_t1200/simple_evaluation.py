"""
LibInvent结果简化评价 - 保证能跑的版本
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, Lipinski
import matplotlib.pyplot as plt
from tqdm import tqdm
import os

# 配置
RESULTS_FILE = "results_1.csv"
OUTPUT_DIR = "simple_evaluation"
SAMPLE_SIZE = 10000  # 只分析1万个，快速运行

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*70)
print("LibInvent 结果快速评价")
print("="*70)

# ============================================================================
# 1. 读取数据
# ============================================================================
print("\n[1/4] 读取数据...")
df = pd.read_csv(RESULTS_FILE)
print(f"✓ 总行数: {len(df)}")

# 提取有效SMILES
valid_smiles = []
valid_scores = []

for idx, row in tqdm(df.iterrows(), total=min(len(df), SAMPLE_SIZE), desc="解析"):
    if idx >= SAMPLE_SIZE:
        break
    
    smi = row['SMILES']
    if pd.isna(smi):
        continue
    
    smi_str = str(smi).strip()
    
    # 跳过片段
    if '[*]' in smi_str or '|' in smi_str:
        continue
    
    # 验证SMILES
    mol = Chem.MolFromSmiles(smi_str)
    if mol:
        valid_smiles.append(smi_str)
        valid_scores.append(row.get('Score', 0))

print(f"✓ 有效SMILES: {len(valid_smiles)}")

# ============================================================================
# 2. 计算基础指标
# ============================================================================
print("\n[2/4] 计算指标...")

results = []

for smi in tqdm(valid_smiles, desc="计算"):
    mol = Chem.MolFromSmiles(smi)
    
    try:
        results.append({
            'SMILES': smi,
            'QED': QED.qed(mol),
            'MW': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'ArRings': Descriptors.NumAromaticRings(mol)
        })
    except Exception as e:
        continue

df_results = pd.DataFrame(results)
df_results['Score'] = valid_scores[:len(df_results)]

print(f"✓ 计算完成: {len(df_results)} 个分子")

# ============================================================================
# 3. 统计分析
# ============================================================================
print("\n[3/4] 统计分析...")

# RO5 violations
df_results['RO5_Viol'] = (
    (df_results['MW'] > 500).astype(int) +
    (df_results['LogP'] > 5).astype(int) +
    (df_results['HBA'] > 10).astype(int) +
    (df_results['HBD'] > 5).astype(int)
)

report = []
report.append("="*70)
report.append("LibInvent 生成结果评价")
report.append("="*70)
report.append("")
report.append(f"分析样本: {len(df_results)} 个分子")
report.append("")

report.append("## 成药性指标")
report.append(f"  QED:     {df_results['QED'].mean():.3f} ± {df_results['QED'].std():.3f}")
report.append(f"  MW:      {df_results['MW'].mean():.1f} ± {df_results['MW'].std():.1f}")
report.append(f"  LogP:    {df_results['LogP'].mean():.2f} ± {df_results['LogP'].std():.2f}")
report.append(f"  TPSA:    {df_results['TPSA'].mean():.1f} ± {df_results['TPSA'].std():.1f}")
report.append("")

report.append("## 质量评估")
qed_good = (df_results['QED'] > 0.5).sum()
qed_exc = (df_results['QED'] > 0.7).sum()
report.append(f"  QED > 0.5:  {qed_good} ({qed_good/len(df_results)*100:.1f}%)")
report.append(f"  QED > 0.7:  {qed_exc} ({qed_exc/len(df_results)*100:.1f}%)")
report.append("")

report.append("## Lipinski Rule of 5")
for i in range(5):
    count = (df_results['RO5_Viol'] == i).sum()
    pct = count / len(df_results) * 100
    report.append(f"  {i} 违反: {count} ({pct:.1f}%)")
report.append("")

# 高质量分子
high_quality = df_results[
    (df_results['QED'] > 0.6) &
    (df_results['RO5_Viol'] <= 1) &
    (df_results['MW'].between(200, 600))
]

report.append("## 高质量分子")
report.append(f"  数量: {len(high_quality)} ({len(high_quality)/len(df_results)*100:.1f}%)")
report.append(f"  定义: QED>0.6, RO5≤1, 200<MW<600")
report.append("")

# 评级
if qed_good/len(df_results) > 0.7:
    rating = "优秀 ⭐⭐⭐"
elif qed_good/len(df_results) > 0.5:
    rating = "良好 ⭐⭐"
else:
    rating = "一般 ⭐"

report.append(f"## 总体评级: {rating}")
report.append("")
report.append("="*70)

report_text = "\n".join(report)
print("\n" + report_text)

# 保存报告
with open(f"{OUTPUT_DIR}/report.txt", 'w', encoding='utf-8') as f:
    f.write(report_text)

# 保存数据
df_results.to_csv(f"{OUTPUT_DIR}/metrics.csv", index=False)
high_quality.to_csv(f"{OUTPUT_DIR}/high_quality.csv", index=False)

print(f"\n✓ 报告已保存: {OUTPUT_DIR}/report.txt")
print(f"✓ 数据已保存: {OUTPUT_DIR}/metrics.csv")
print(f"✓ 高质量分子: {OUTPUT_DIR}/high_quality.csv ({len(high_quality)} 个)")

# ============================================================================
# 4. 可视化
# ============================================================================
print("\n[4/4] 生成图表...")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# 1. QED
axes[0].hist(df_results['QED'], bins=50, color='steelblue', 
             alpha=0.7, edgecolor='black')
axes[0].axvline(0.5, color='orange', linestyle='--', label='Good')
axes[0].axvline(0.7, color='green', linestyle='--', label='Excellent')
axes[0].set_xlabel('QED', fontweight='bold')
axes[0].set_ylabel('Count', fontweight='bold')
axes[0].set_title('Drug-likeness (QED)', fontweight='bold', fontsize=12)
axes[0].legend()
axes[0].grid(alpha=0.3)

# 2. MW
axes[1].hist(df_results['MW'], bins=50, color='coral', 
             alpha=0.7, edgecolor='black')
axes[1].axvline(500, color='red', linestyle='--', label='RO5 limit')
axes[1].set_xlabel('Molecular Weight', fontweight='bold')
axes[1].set_ylabel('Count', fontweight='bold')
axes[1].set_title('Molecular Weight', fontweight='bold', fontsize=12)
axes[1].legend()
axes[1].grid(alpha=0.3)

# 3. LogP
axes[2].hist(df_results['LogP'], bins=50, color='lightgreen', 
             alpha=0.7, edgecolor='black')
axes[2].axvline(5, color='red', linestyle='--', label='RO5 limit')
axes[2].set_xlabel('LogP', fontweight='bold')
axes[2].set_ylabel('Count', fontweight='bold')
axes[2].set_title('Lipophilicity', fontweight='bold', fontsize=12)
axes[2].legend()
axes[2].grid(alpha=0.3)

# 4. TPSA
axes[3].hist(df_results['TPSA'], bins=50, color='plum', 
             alpha=0.7, edgecolor='black')
axes[3].axvline(140, color='red', linestyle='--', label='Good')
axes[3].set_xlabel('TPSA (Ų)', fontweight='bold')
axes[3].set_ylabel('Count', fontweight='bold')
axes[3].set_title('Polar Surface Area', fontweight='bold', fontsize=12)
axes[3].legend()
axes[3].grid(alpha=0.3)

# 5. RO5
ro5_counts = df_results['RO5_Viol'].value_counts().sort_index()
axes[4].bar(ro5_counts.index, ro5_counts.values, 
           color='skyblue', edgecolor='black')
axes[4].set_xlabel('Violations', fontweight='bold')
axes[4].set_ylabel('Count', fontweight='bold')
axes[4].set_title('Lipinski Rule of 5', fontweight='bold', fontsize=12)
axes[4].grid(alpha=0.3, axis='y')

# 6. QED vs MW
scatter = axes[5].scatter(df_results['MW'], df_results['QED'], 
                         c=df_results['Score'], cmap='viridis', 
                         alpha=0.6, s=20)
axes[5].axvline(500, color='red', linestyle='--', alpha=0.5)
axes[5].axhline(0.5, color='orange', linestyle='--', alpha=0.5)
axes[5].set_xlabel('MW', fontweight='bold')
axes[5].set_ylabel('QED', fontweight='bold')
axes[5].set_title('QED vs MW', fontweight='bold', fontsize=12)
plt.colorbar(scatter, ax=axes[5], label='Score')
axes[5].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/plots.png", dpi=300, bbox_inches='tight')
print(f"✓ 图表已保存: {OUTPUT_DIR}/plots.png")

plt.close()

print("\n" + "="*70)
print("✅ 评价完成！")
print("="*70)
print(f"\n查看结果:")
print(f"  - {OUTPUT_DIR}/report.txt")
print(f"  - {OUTPUT_DIR}/plots.png")
print(f"  - {OUTPUT_DIR}/high_quality.csv ({len(high_quality)} 个高质量分子)")
print()