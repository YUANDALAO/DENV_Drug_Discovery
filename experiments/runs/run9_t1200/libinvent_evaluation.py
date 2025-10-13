"""
LibInvent生成结果综合评价 - 完整版
包含：成药性、合成可行性、多样性、官能团分析
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED, Crippen, Lipinski
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import DataStructs
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from collections import Counter
import warnings
import os
warnings.filterwarnings('ignore')

# ============================================================================
# 配置
# ============================================================================
RESULTS_FILE = "results_1.csv"
OUTPUT_DIR = "libinvent_evaluation"
SAMPLE_SIZE_DRUGLIKENESS = 50000  # 成药性分析样本
SAMPLE_SIZE_DIVERSITY = 5000      # 多样性分析样本

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 核心计算函数
# ============================================================================

def calculate_druglikeness(mol):
    """计算成药性指标"""
    try:
        metrics = {
            'QED': QED.qed(mol),
            'MW': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'ArRings': Descriptors.NumAromaticRings(mol),
            'HeavyAtoms': Lipinski.HeavyAtomCount(mol)
        }
        
        # RO5 violations
        ro5 = 0
        if metrics['MW'] > 500: ro5 += 1
        if metrics['LogP'] > 5: ro5 += 1
        if metrics['HBA'] > 10: ro5 += 1
        if metrics['HBD'] > 5: ro5 += 1
        metrics['RO5_Violations'] = ro5
        
        return metrics
    except:
        return None


def calculate_sa_score(mol):
    """计算合成难度 - 使用简化版本"""
    try:
        # 尝试导入SA Score
        from rdkit.Chem import RDConfig
        import sys
        sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
        import sascorer
        return sascorer.calculateScore(mol)
    except:
        # 如果不可用，返回简化估算
        # 基于环系统复杂度和手性中心数量
        try:
            rings = Descriptors.RingCount(mol)
            stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            rotbonds = Descriptors.NumRotatableBonds(mol)
            
            # 简化SA评分 (1-10)
            sa_simple = 1 + (rings * 0.5) + (stereo * 0.8) + (rotbonds * 0.1)
            return min(10, max(1, sa_simple))
        except:
            return None


def analyze_functional_groups(mol):
    """分析常见官能团"""
    fg_patterns = {
        'Alcohol': '[OH]',
        'Amine': '[NX3;H2,H1;!$(NC=O)]',
        'Amide': 'C(=O)N',
        'Ester': 'C(=O)O[C,c]',
        'Ether': '[OD2]([C,c])[C,c]',
        'Ketone': '[CX3](=O)[C,c]',
        'Halogen': '[F,Cl,Br,I]',
        'Aromatic_N': '[n]',
        'Sulfone': 'S(=O)(=O)',
        'Nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])]'
    }
    
    counts = {}
    for name, smarts in fg_patterns.items():
        try:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                counts[name] = len(matches)
            else:
                counts[name] = 0
        except:
            counts[name] = 0
    
    return counts


def calculate_diversity(smiles_list, max_samples=5000):
    """计算结构多样性指标"""
    
    # 采样
    if len(smiles_list) > max_samples:
        import random
        smiles_sample = random.sample(smiles_list, max_samples)
    else:
        smiles_sample = smiles_list
    
    print(f"\n[多样性分析] 采样 {len(smiles_sample)} 个分子")
    
    # 计算指纹
    fps = []
    scaffolds = []
    
    for smi in tqdm(smiles_sample, desc="计算指纹和骨架"):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            # Morgan指纹
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fps.append(fp)
            
            # Murcko骨架
            try:
                scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                scaffolds.append(scaffold)
            except:
                scaffolds.append(None)
    
    # 计算内部相似度
    print("  计算内部相似度...")
    n_pairs = min(10000, len(fps) * (len(fps) - 1) // 2)
    similarities = []
    
    import random
    for _ in tqdm(range(n_pairs), desc="  采样相似度对"):
        i, j = random.sample(range(len(fps)), 2)
        sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        similarities.append(sim)
    
    # 骨架统计
    scaffold_counts = Counter([s for s in scaffolds if s])
    
    metrics = {
        'n_molecules': len(smiles_sample),
        'mean_similarity': np.mean(similarities),
        'median_similarity': np.median(similarities),
        'unique_scaffolds': len(scaffold_counts),
        'scaffold_diversity': len(scaffold_counts) / len(scaffolds) if scaffolds else 0,
        'top10_coverage': sum(c for _, c in scaffold_counts.most_common(10)) / len(scaffolds) if scaffolds else 0
    }
    
    return metrics, similarities, scaffold_counts


# ============================================================================
# 主评价流程
# ============================================================================

def main():
    
    print("="*80)
    print("LibInvent 生成结果综合评价 - 完整版")
    print("="*80)
    
    # ========================================================================
    # Step 1: 读取数据
    # ========================================================================
    print("\nStep 1/5: 读取生成结果")
    print("-"*80)
    
    df = pd.read_csv(RESULTS_FILE)
    print(f"✓ 文件行数: {len(df)}")
    
    # 提取有效SMILES
    valid_data = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="解析SMILES"):
        smi = row['SMILES']
        
        if pd.isna(smi):
            continue
        
        smi_str = str(smi).strip()
        
        # 跳过片段格式
        if '[*]' in smi_str or '|' in smi_str or len(smi_str) > 500:
            continue
        
        mol = Chem.MolFromSmiles(smi_str)
        if mol:
            valid_data.append({
                'SMILES': smi_str,
                'Score': row.get('Score', 0),
                'mol': mol
            })
    
    print(f"✓ 有效分子: {len(valid_data)} ({len(valid_data)/len(df)*100:.1f}%)")
    
    all_smiles = [d['SMILES'] for d in valid_data]
    unique_smiles = list(set(all_smiles))
    print(f"✓ 独特分子: {len(unique_smiles)} ({len(unique_smiles)/len(all_smiles)*100:.1f}%)")
    
    # ========================================================================
    # Step 2: 成药性分析
    # ========================================================================
    print(f"\nStep 2/5: 成药性分析")
    print("-"*80)
    
    sample_dl = valid_data[:min(SAMPLE_SIZE_DRUGLIKENESS, len(valid_data))]
    print(f"分析样本: {len(sample_dl)} 个分子")
    
    druglikeness_results = []
    sa_scores = []
    
    for item in tqdm(sample_dl, desc="计算指标"):
        mol = item['mol']
        
        dl_metrics = calculate_druglikeness(mol)
        if dl_metrics:
            dl_metrics['SMILES'] = item['SMILES']
            dl_metrics['Score'] = item['Score']
            druglikeness_results.append(dl_metrics)
            
            sa = calculate_sa_score(mol)
            if sa:
                sa_scores.append(sa)
    
    df_dl = pd.DataFrame(druglikeness_results)
    print(f"✓ 成药性计算完成: {len(df_dl)} 个")
    
    # ========================================================================
    # Step 3: 官能团分析
    # ========================================================================
    print(f"\nStep 3/5: 官能团分析")
    print("-"*80)
    
    fg_data = []
    sample_fg = valid_data[:min(10000, len(valid_data))]
    
    for item in tqdm(sample_fg, desc="分析官能团"):
        fg = analyze_functional_groups(item['mol'])
        fg['SMILES'] = item['SMILES']
        fg_data.append(fg)
    
    df_fg = pd.DataFrame(fg_data)
    print(f"✓ 官能团分析完成: {len(df_fg)} 个")
    
    # ========================================================================
    # Step 4: 多样性分析
    # ========================================================================
    print(f"\nStep 4/5: 多样性分析")
    print("-"*80)
    
    diversity_metrics, similarities, scaffold_counts = calculate_diversity(
        all_smiles, 
        max_samples=SAMPLE_SIZE_DIVERSITY
    )
    
    print(f"✓ 多样性计算完成")
    print(f"  平均内部相似度: {diversity_metrics['mean_similarity']:.3f}")
    print(f"  独特骨架数: {diversity_metrics['unique_scaffolds']}")
    
    # ========================================================================
    # Step 5: 生成报告
    # ========================================================================
    print(f"\nStep 5/5: 生成报告和可视化")
    print("-"*80)
    
    # 生成文本报告
    report = generate_text_report(
        df, df_dl, df_fg, sa_scores, 
        diversity_metrics, scaffold_counts
    )
    
    with open(f"{OUTPUT_DIR}/evaluation_report.txt", 'w', encoding='utf-8') as f:
        f.write(report)
    
    print(f"✓ 报告已保存: {OUTPUT_DIR}/evaluation_report.txt")
    
    # 保存数据
    df_dl.to_csv(f"{OUTPUT_DIR}/druglikeness_metrics.csv", index=False)
    df_fg.to_csv(f"{OUTPUT_DIR}/functional_groups.csv", index=False)
    
    # 保存高质量分子
    high_quality = df_dl[
        (df_dl['QED'] > 0.6) &
        (df_dl['RO5_Violations'] <= 1) &
        (df_dl['MW'].between(200, 600))
    ]
    high_quality.to_csv(f"{OUTPUT_DIR}/high_quality_molecules.csv", index=False)
    print(f"✓ 高质量分子: {OUTPUT_DIR}/high_quality_molecules.csv ({len(high_quality)} 个)")
    
    # 生成可视化
    generate_plots(df_dl, df_fg, sa_scores, similarities)
    print(f"✓ 图表已保存: {OUTPUT_DIR}/evaluation_plots.png")
    
    # 打印报告
    print("\n" + "="*80)
    print(report)
    print("="*80)
    
    print(f"\n✅ 评价完成！所有结果保存在: {OUTPUT_DIR}/")


def generate_text_report(df_raw, df_dl, df_fg, sa_scores, div_metrics, scaffolds):
    """生成文本报告"""
    
    report = []
    report.append("="*80)
    report.append("LibInvent 生成结果评价报告")
    report.append("="*80)
    report.append("")
    
    # 1. 总体统计
    report.append("## 1. 总体统计")
    report.append(f"   原始生成数: {len(df_raw)}")
    report.append(f"   有效分子数: {len(df_dl)}")
    report.append(f"   有效率: {len(df_dl)/len(df_raw)*100:.1f}%")
    report.append("")
    
    # 2. 成药性
    report.append("## 2. 成药性评估")
    report.append(f"   平均QED: {df_dl['QED'].mean():.3f} ± {df_dl['QED'].std():.3f}")
    qed_good = (df_dl['QED'] > 0.5).sum()
    qed_exc = (df_dl['QED'] > 0.7).sum()
    report.append(f"   QED > 0.5: {qed_good} ({qed_good/len(df_dl)*100:.1f}%)")
    report.append(f"   QED > 0.7: {qed_exc} ({qed_exc/len(df_dl)*100:.1f}%)")
    report.append("")
    
    report.append("   物化性质:")
    report.append(f"   MW:   {df_dl['MW'].mean():.1f} ± {df_dl['MW'].std():.1f}")
    report.append(f"   LogP: {df_dl['LogP'].mean():.2f} ± {df_dl['LogP'].std():.2f}")
    report.append(f"   TPSA: {df_dl['TPSA'].mean():.1f} ± {df_dl['TPSA'].std():.1f}")
    report.append(f"   HBA:  {df_dl['HBA'].mean():.1f} ± {df_dl['HBA'].std():.1f}")
    report.append(f"   HBD:  {df_dl['HBD'].mean():.1f} ± {df_dl['HBD'].std():.1f}")
    report.append("")
    
    report.append("   Lipinski Rule of 5:")
    for v in range(5):
        cnt = (df_dl['RO5_Violations'] == v).sum()
        pct = cnt / len(df_dl) * 100
        report.append(f"   {v}违反: {cnt} ({pct:.1f}%)")
    report.append("")
    
    # 3. 合成可行性
    if sa_scores:
        report.append("## 3. 合成可行性")
        report.append(f"   平均SA Score: {np.mean(sa_scores):.2f} ± {np.std(sa_scores):.2f}")
        easy = sum(1 for s in sa_scores if s < 3)
        medium = sum(1 for s in sa_scores if 3 <= s < 6)
        hard = sum(1 for s in sa_scores if s >= 6)
        report.append(f"   容易(SA<3): {easy} ({easy/len(sa_scores)*100:.1f}%)")
        report.append(f"   中等(3≤SA<6): {medium} ({medium/len(sa_scores)*100:.1f}%)")
        report.append(f"   困难(SA≥6): {hard} ({hard/len(sa_scores)*100:.1f}%)")
        report.append("")
    
    # 4. 多样性
    report.append("## 4. 结构多样性")
    report.append(f"   内部平均相似度: {div_metrics['mean_similarity']:.3f}")
    report.append(f"   独特骨架数: {div_metrics['unique_scaffolds']}")
    report.append(f"   骨架多样性指数: {div_metrics['scaffold_diversity']:.3f}")
    report.append(f"   Top10骨架覆盖: {div_metrics['top10_coverage']*100:.1f}%")
    report.append("")
    
    report.append("   最常见的10个骨架:")
    for i, (scf, cnt) in enumerate(scaffolds.most_common(10), 1):
        pct = cnt / sum(scaffolds.values()) * 100
        report.append(f"   {i}. {scf[:60]} ({cnt}, {pct:.1f}%)")
    report.append("")
    
    # 5. 官能团
    report.append("## 5. 官能团分布")
    fg_means = df_fg.drop('SMILES', axis=1, errors='ignore').mean()
    fg_top = fg_means.sort_values(ascending=False).head(10)
    for fg, cnt in fg_top.items():
        report.append(f"   {fg}: {cnt:.2f} 个/分子")
    report.append("")
    
    # 6. 质量评估
    report.append("## 6. 综合质量评估")
    high_quality = df_dl[
        (df_dl['QED'] > 0.6) &
        (df_dl['RO5_Violations'] <= 1) &
        (df_dl['MW'].between(200, 600))
    ]
    report.append(f"   高质量分子: {len(high_quality)} ({len(high_quality)/len(df_dl)*100:.1f}%)")
    report.append(f"   (定义: QED>0.6, RO5≤1, 200<MW<600)")
    report.append("")
    
    # 评级
    qed_good_pct = qed_good / len(df_dl) * 100
    if qed_good_pct > 70:
        dl_rating = "优秀 ⭐⭐⭐"
    elif qed_good_pct > 50:
        dl_rating = "良好 ⭐⭐"
    else:
        dl_rating = "一般 ⭐"
    
    if div_metrics['mean_similarity'] < 0.4:
        div_rating = "优秀 ⭐⭐⭐"
    elif div_metrics['mean_similarity'] < 0.6:
        div_rating = "良好 ⭐⭐"
    else:
        div_rating = "一般 ⭐"
    
    report.append("   评级:")
    report.append(f"   成药性: {dl_rating}")
    report.append(f"   多样性: {div_rating}")
    report.append("")
    
    # 7. 建议
    report.append("## 7. 改进建议")
    if qed_good_pct < 50:
        report.append("   ⚠️ 成药性偏低，建议:")
        report.append("      - 增加QED权重到2.0以上")
        report.append("      - 限制MW上限到450")
        report.append("      - 控制LogP在2-4范围")
    
    if div_metrics['mean_similarity'] > 0.6:
        report.append("   ⚠️ 多样性不足，建议:")
        report.append("      - 增加多样性惩罚")
        report.append("      - 扩展R-group库")
        report.append("      - 尝试不同骨架变体")
    
    report.append("")
    report.append("="*80)
    
    return "\n".join(report)


def generate_plots(df_dl, df_fg, sa_scores, similarities):
    """生成可视化图表"""
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    axes = axes.flatten()
    
    # 1. QED
    axes[0].hist(df_dl['QED'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    axes[0].axvline(0.5, color='orange', linestyle='--', linewidth=2, label='Good')
    axes[0].axvline(0.7, color='green', linestyle='--', linewidth=2, label='Excellent')
    axes[0].set_xlabel('QED', fontweight='bold', fontsize=11)
    axes[0].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[0].set_title('Drug-likeness (QED)', fontweight='bold', fontsize=12)
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # 2. MW
    axes[1].hist(df_dl['MW'], bins=50, color='coral', alpha=0.7, edgecolor='black')
    axes[1].axvline(500, color='red', linestyle='--', linewidth=2, label='RO5')
    axes[1].set_xlabel('Molecular Weight', fontweight='bold', fontsize=11)
    axes[1].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[1].set_title('Molecular Weight', fontweight='bold', fontsize=12)
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    # 3. LogP
    axes[2].hist(df_dl['LogP'], bins=50, color='lightgreen', alpha=0.7, edgecolor='black')
    axes[2].axvline(5, color='red', linestyle='--', linewidth=2, label='RO5')
    axes[2].set_xlabel('LogP', fontweight='bold', fontsize=11)
    axes[2].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[2].set_title('Lipophilicity', fontweight='bold', fontsize=12)
    axes[2].legend()
    axes[2].grid(alpha=0.3)
    
    # 4. TPSA
    axes[3].hist(df_dl['TPSA'], bins=50, color='plum', alpha=0.7, edgecolor='black')
    axes[3].axvline(140, color='red', linestyle='--', linewidth=2)
    axes[3].set_xlabel('TPSA (Ų)', fontweight='bold', fontsize=11)
    axes[3].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[3].set_title('Polar Surface Area', fontweight='bold', fontsize=12)
    axes[3].grid(alpha=0.3)
    
    # 5. RO5
    ro5_counts = df_dl['RO5_Violations'].value_counts().sort_index()
    axes[4].bar(ro5_counts.index, ro5_counts.values, color='skyblue', edgecolor='black')
    axes[4].set_xlabel('Violations', fontweight='bold', fontsize=11)
    axes[4].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[4].set_title('Lipinski Rule of 5', fontweight='bold', fontsize=12)
    axes[4].grid(alpha=0.3, axis='y')
    
    # 6. SA Score
    if sa_scores:
        axes[5].hist(sa_scores, bins=50, color='wheat', alpha=0.7, edgecolor='black')
        axes[5].axvline(3, color='green', linestyle='--', linewidth=2, label='Easy')
        axes[5].axvline(6, color='orange', linestyle='--', linewidth=2, label='Hard')
        axes[5].set_xlabel('SA Score', fontweight='bold', fontsize=11)
        axes[5].set_ylabel('Count', fontweight='bold', fontsize=11)
        axes[5].set_title('Synthetic Accessibility', fontweight='bold', fontsize=12)
        axes[5].legend()
        axes[5].grid(alpha=0.3)
    
    # 7. 内部相似度
    axes[6].hist(similarities, bins=50, color='lightcoral', alpha=0.7, edgecolor='black')
    axes[6].axvline(np.mean(similarities), color='red', linestyle='--', linewidth=2,
                   label=f'Mean={np.mean(similarities):.3f}')
    axes[6].set_xlabel('Internal Tanimoto', fontweight='bold', fontsize=11)
    axes[6].set_ylabel('Count', fontweight='bold', fontsize=11)
    axes[6].set_title('Structural Diversity', fontweight='bold', fontsize=12)
    axes[6].legend()
    axes[6].grid(alpha=0.3)
    
    # 8. QED vs MW
    scatter = axes[7].scatter(df_dl['MW'], df_dl['QED'], 
                             c=df_dl['Score'], cmap='viridis', 
                             alpha=0.5, s=15)
    axes[7].axvline(500, color='red', linestyle='--', alpha=0.5, linewidth=2)
    axes[7].axhline(0.5, color='orange', linestyle='--', alpha=0.5, linewidth=2)
    axes[7].set_xlabel('MW', fontweight='bold', fontsize=11)
    axes[7].set_ylabel('QED', fontweight='bold', fontsize=11)
    axes[7].set_title('QED vs MW', fontweight='bold', fontsize=12)
    plt.colorbar(scatter, ax=axes[7], label='Score')
    axes[7].grid(alpha=0.3)
    
    # 9. 官能团
    fg_means = df_fg.drop('SMILES', axis=1, errors='ignore').mean()
    fg_top10 = fg_means.sort_values(ascending=False).head(10)
    axes[8].barh(range(len(fg_top10)), fg_top10.values, color='teal', edgecolor='black')
    axes[8].set_yticks(range(len(fg_top10)))
    axes[8].set_yticklabels(fg_top10.index, fontsize=9)
    axes[8].set_xlabel('Avg Count/Molecule', fontweight='bold', fontsize=11)
    axes[8].set_title('Top 10 Functional Groups', fontweight='bold', fontsize=12)
    axes[8].grid(alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/evaluation_plots.png", dpi=300, bbox_inches='tight')
    plt.close()


# ============================================================================
# 运行
# ============================================================================

if __name__ == "__main__":
    main()