"""
Top 500候选分子完整分析
包含：化学空间可视化、结构图生成、Score=0分析
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from sklearn.decomposition import PCA
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("Top 500 候选分子深度分析")
print("="*80)

# ============================================================================
# 1. 化学空间可视化 (PCA)
# ============================================================================
print("\n[1/4] 化学空间分析...")

top500 = pd.read_csv("top500_candidates.csv")
print(f"✓ 读取 {len(top500)} 个分子")

# 计算Morgan指纹
fps = []
valid_mols = []

for smi in tqdm(top500['SMILES'], desc="计算指纹"):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fps.append(np.array(list(fp)))
        valid_mols.append(mol)

fps_matrix = np.array(fps)
print(f"✓ 有效分子: {len(valid_mols)}")

# PCA降维
pca = PCA(n_components=2)
coords = pca.fit_transform(fps_matrix)

print(f"✓ PCA完成")
print(f"  解释方差: PC1={pca.explained_variance_ratio_[0]*100:.1f}%, PC2={pca.explained_variance_ratio_[1]*100:.1f}%")

# 可视化
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# 左图：按QED着色
ax1 = axes[0]
scatter1 = ax1.scatter(coords[:, 0], coords[:, 1], 
                      c=top500['QED'][:len(coords)], 
                      cmap='RdYlGn', s=50, alpha=0.7, 
                      edgecolors='black', linewidth=0.5)
ax1.set_xlabel('PC1', fontweight='bold', fontsize=12)
ax1.set_ylabel('PC2', fontweight='bold', fontsize=12)
ax1.set_title('Chemical Space (colored by QED)', fontweight='bold', fontsize=13)
ax1.grid(alpha=0.3)
cbar1 = plt.colorbar(scatter1, ax=ax1)
cbar1.set_label('QED', fontweight='bold')

# 右图：按MW着色
ax2 = axes[1]
scatter2 = ax2.scatter(coords[:, 0], coords[:, 1], 
                      c=top500['MW'][:len(coords)], 
                      cmap='viridis', s=50, alpha=0.7, 
                      edgecolors='black', linewidth=0.5)
ax2.set_xlabel('PC1', fontweight='bold', fontsize=12)
ax2.set_ylabel('PC2', fontweight='bold', fontsize=12)
ax2.set_title('Chemical Space (colored by MW)', fontweight='bold', fontsize=13)
ax2.grid(alpha=0.3)
cbar2 = plt.colorbar(scatter2, ax=ax2)
cbar2.set_label('Molecular Weight', fontweight='bold')

plt.tight_layout()
plt.savefig('top500_chemical_space.png', dpi=300, bbox_inches='tight')
print("✓ 化学空间图已保存: top500_chemical_space.png")
plt.close()

# ============================================================================
# 2. 生成Top 50结构图
# ============================================================================
print("\n[2/4] 生成结构图...")

mols_for_grid = []
legends = []

for idx, row in top500.head(50).iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols_for_grid.append(mol)
        legend = f"#{len(mols_for_grid)}\nQED={row['QED']:.3f}\nMW={row['MW']:.0f}\nLogP={row['LogP']:.2f}"
        legends.append(legend)

img = Draw.MolsToGridImage(
    mols_for_grid, 
    molsPerRow=5, 
    subImgSize=(350, 350),
    legends=legends,
    returnPNG=False
)

img.save("top50_structures.png")
print(f"✓ Top 50结构图已保存: top50_structures.png")

# ============================================================================
# 3. Score=0分析
# ============================================================================
print("\n[3/4] Score=0异常分析...")

score_zero = top500[top500['Score'] == 0]
print(f"✓ Score=0的分子: {len(score_zero)} 个")

if len(score_zero) > 0:
    print("\nScore=0但QED高的分子:")
    print(score_zero[['SMILES', 'QED', 'MW', 'LogP', 'TPSA']].head(10).to_string())
    
    # 查看原始数据
    print("\n正在查找原始数据中的信息...")
    try:
        df_full = pd.read_csv("results_1.csv")
        
        # 找出这些分子在原始数据中的行
        for idx, row in score_zero.head(5).iterrows():
            smi = row['SMILES']
            matches = df_full[df_full['SMILES'] == smi]
            
            if len(matches) > 0:
                print(f"\n分子: {smi[:50]}...")
                print(f"  出现次数: {len(matches)}")
                print(f"  Scores: {matches['Score'].values}")
                
                # 检查是否有scoring组件信息
                if 'DENV_Activity (raw)' in matches.columns:
                    print(f"  DENV_Activity (raw): {matches['DENV_Activity (raw)'].values[0]}")
                if 'Stability_Alerts' in matches.columns:
                    print(f"  Stability_Alerts: {matches['Stability_Alerts'].values[0]}")
    except Exception as e:
        print(f"  无法读取原始数据: {e}")

# ============================================================================
# 4. 统计分析
# ============================================================================
print("\n[4/4] 详细统计分析...")

# 计算内部相似度
print("\n计算Top 500内部相似度...")
from rdkit.Chem import DataStructs

similarities = []
n_samples = min(1000, len(fps) * (len(fps) - 1) // 2)

import random
random.seed(42)

for _ in tqdm(range(n_samples), desc="采样相似度对"):
    i, j = random.sample(range(len(fps)), 2)
    fp_i = [fps[i]]
    fp_j = [fps[j]]
    
    # 重新生成RDKit指纹对象
    mol_i = valid_mols[i]
    mol_j = valid_mols[j]
    fp_i_rd = AllChem.GetMorganFingerprintAsBitVect(mol_i, 2, nBits=2048)
    fp_j_rd = AllChem.GetMorganFingerprintAsBitVect(mol_j, 2, nBits=2048)
    
    sim = DataStructs.TanimotoSimilarity(fp_i_rd, fp_j_rd)
    similarities.append(sim)

# 生成统计报告
report = []
report.append("="*80)
report.append("Top 500 候选分子分析报告")
report.append("="*80)
report.append("")

report.append("## 1. 基本统计")
report.append(f"   分子数: {len(top500)}")
report.append(f"   QED范围: {top500['QED'].min():.4f} - {top500['QED'].max():.4f}")
report.append(f"   QED平均: {top500['QED'].mean():.4f} ± {top500['QED'].std():.4f}")
report.append(f"   MW范围: {top500['MW'].min():.1f} - {top500['MW'].max():.1f}")
report.append(f"   MW平均: {top500['MW'].mean():.1f} ± {top500['MW'].std():.1f}")
report.append(f"   LogP范围: {top500['LogP'].min():.2f} - {top500['LogP'].max():.2f}")
report.append(f"   LogP平均: {top500['LogP'].mean():.2f} ± {top500['LogP'].std():.2f}")
report.append("")

report.append("## 2. 多样性")
report.append(f"   内部平均相似度: {np.mean(similarities):.4f}")
report.append(f"   相似度中位数: {np.median(similarities):.4f}")
report.append(f"   相似度范围: {np.min(similarities):.4f} - {np.max(similarities):.4f}")
report.append("   (值越低=越多样)")
report.append("")

report.append("## 3. 成药性")
all_good = (top500['QED'] > 0.9).sum()
report.append(f"   QED > 0.9: {all_good} ({all_good/len(top500)*100:.1f}%)")
mw_good = ((top500['MW'] >= 300) & (top500['MW'] <= 500)).sum()
report.append(f"   MW 300-500: {mw_good} ({mw_good/len(top500)*100:.1f}%)")
logp_good = ((top500['LogP'] >= 2) & (top500['LogP'] <= 4)).sum()
report.append(f"   LogP 2-4: {logp_good} ({logp_good/len(top500)*100:.1f}%)")
report.append("")

report.append("## 4. 异常")
report.append(f"   Score=0: {len(score_zero)} ({len(score_zero)/len(top500)*100:.1f}%)")
if len(score_zero) > 0:
    report.append("   这些分子虽然Score=0，但QED极高")
    report.append("   可能原因:")
    report.append("     - 训练早期生成（某些约束未满足）")
    report.append("     - QSAR预测失败")
    report.append("     - 需要检查原始数据确认")
    report.append("   建议: 保留用于实验验证")
report.append("")

report.append("## 5. 推荐策略")
report.append("   下一步:")
if len(score_zero) > 0:
    report.append(f"   1. 优先验证Score>0的 {len(top500)-len(score_zero)} 个分子")
    report.append(f"   2. Score=0的 {len(score_zero)} 个作为补充（QED极高）")
else:
    report.append("   1. 全部500个都可考虑")

report.append("   2. 如有对接条件，对全部500个做对接筛选")
report.append("   3. 选Top 50采购/合成")
report.append("   4. 实验验证20-30个")
report.append("")
report.append("="*80)

report_text = "\n".join(report)
print("\n" + report_text)

# 保存报告
with open("top500_analysis_report.txt", 'w', encoding='utf-8') as f:
    f.write(report_text)

print("\n✓ 报告已保存: top500_analysis_report.txt")

# ============================================================================
# 5. 生成汇总图
# ============================================================================
print("\n生成汇总图...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. QED分布
axes[0, 0].hist(top500['QED'], bins=30, color='steelblue', alpha=0.7, edgecolor='black')
axes[0, 0].axvline(top500['QED'].mean(), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean={top500["QED"].mean():.3f}')
axes[0, 0].set_xlabel('QED', fontweight='bold')
axes[0, 0].set_ylabel('Count', fontweight='bold')
axes[0, 0].set_title('QED Distribution (Top 500)', fontweight='bold')
axes[0, 0].legend()
axes[0, 0].grid(alpha=0.3)

# 2. MW vs QED
scatter = axes[0, 1].scatter(top500['MW'], top500['QED'], 
                            c=top500['LogP'], cmap='viridis',
                            s=30, alpha=0.6, edgecolors='black')
axes[0, 1].set_xlabel('Molecular Weight', fontweight='bold')
axes[0, 1].set_ylabel('QED', fontweight='bold')
axes[0, 1].set_title('MW vs QED (colored by LogP)', fontweight='bold')
plt.colorbar(scatter, ax=axes[0, 1], label='LogP')
axes[0, 1].grid(alpha=0.3)

# 3. LogP分布
axes[1, 0].hist(top500['LogP'], bins=30, color='lightgreen', alpha=0.7, edgecolor='black')
axes[1, 0].axvline(top500['LogP'].mean(), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean={top500["LogP"].mean():.2f}')
axes[1, 0].axvline(2, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
axes[1, 0].axvline(4, color='orange', linestyle=':', linewidth=1.5, alpha=0.7, 
                   label='Ideal range (2-4)')
axes[1, 0].set_xlabel('LogP', fontweight='bold')
axes[1, 0].set_ylabel('Count', fontweight='bold')
axes[1, 0].set_title('LogP Distribution', fontweight='bold')
axes[1, 0].legend()
axes[1, 0].grid(alpha=0.3)

# 4. 内部相似度
axes[1, 1].hist(similarities, bins=30, color='coral', alpha=0.7, edgecolor='black')
axes[1, 1].axvline(np.mean(similarities), color='red', linestyle='--', 
                   linewidth=2, label=f'Mean={np.mean(similarities):.3f}')
axes[1, 1].set_xlabel('Tanimoto Similarity', fontweight='bold')
axes[1, 1].set_ylabel('Count', fontweight='bold')
axes[1, 1].set_title('Internal Diversity', fontweight='bold')
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('top500_summary.png', dpi=300, bbox_inches='tight')
print("✓ 汇总图已保存: top500_summary.png")

print("\n" + "="*80)
print("✅ 分析完成！")
print("="*80)
print("\n生成的文件:")
print("  1. top500_chemical_space.png - 化学空间PCA图")
print("  2. top50_structures.png - Top 50结构图")
print("  3. top500_summary.png - 统计汇总图")
print("  4. top500_analysis_report.txt - 文本报告")
print()