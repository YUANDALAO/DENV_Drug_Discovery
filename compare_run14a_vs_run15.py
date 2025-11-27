#!/usr/bin/env python3
"""
深度对比分析: run14a vs run15
重点分析毒性组件的影响
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.figsize'] = (16, 10)

def load_gold_candidates(run_dir):
    """加载金标准候选物"""
    import glob
    
    # 查找gold文件
    gold_files = glob.glob(f"{run_dir}/candidates_gold*.csv")
    if not gold_files:
        gold_files = glob.glob(f"{run_dir}/*金标准*.csv")
    if not gold_files:
        gold_files = glob.glob(f"{run_dir}/promising*.csv")
    
    if gold_files:
        df = pd.read_csv(gold_files[0])
        print(f"✓ 加载 {run_dir}: {len(df)} 个gold候选物")
        return df
    else:
        print(f"✗ 未找到 {run_dir} 的gold文件")
        return None

def calculate_molecular_features(df):
    """计算分子特征"""
    features = []
    
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            continue
        
        # 基本描述符
        feat = {
            'SMILES': smiles,
            'MW': Descriptors.MolWt(mol),
            'LogP': Crippen.MolLogP(mol),
            'TPSA': Descriptors.TPSA(mol),
            'HBA': Lipinski.NumHAcceptors(mol),
            'HBD': Lipinski.NumHDonors(mol),
            'RotBonds': Descriptors.NumRotatableBonds(mol),
            'AromaticRings': Descriptors.NumAromaticRings(mol),
            'HeavyAtoms': Lipinski.HeavyAtomCount(mol),
            'FractionCsp3': Lipinski.FractionCsp3(mol),
            'NumRings': Descriptors.RingCount(mol),
            'NumAliphaticRings': Descriptors.NumAliphaticRings(mol),
        }
        
        # 从CSV中提取的数据
        if 'DENV_Activity (raw)' in row:
            feat['pIC50'] = row['DENV_Activity (raw)']
        elif 'DENV_Activity' in row:
            feat['pIC50'] = row['DENV_Activity']
        
        if 'QED (raw)' in row:
            feat['QED'] = row['QED (raw)']
        elif 'QED' in row:
            feat['QED'] = row['QED']
        
        if 'SA (raw)' in row:
            feat['SA'] = row['SA (raw)']
        elif 'SA' in row:
            feat['SA'] = row['SA']
        
        if 'Score' in row:
            feat['TotalScore'] = row['Score']
        
        # 提取骨架
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            feat['Scaffold'] = Chem.MolToSmiles(scaffold)
        except:
            feat['Scaffold'] = None
        
        # 检测可能的毒性相关基团
        feat['HasNitro'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[N+](=O)[O-]'))
        feat['HasAzo'] = mol.HasSubstructMatch(Chem.MolFromSmarts('N=N'))
        feat['HasHalogenatedCarbon'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[F,Cl,Br,I][C,c][F,Cl,Br,I]'))
        feat['HasQuinone'] = mol.HasSubstructMatch(Chem.MolFromSmarts('C1=CC(=O)C=CC1=O'))
        feat['HasMichaelAcceptor'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[C]=[C]-[C]=[O,S]'))
        feat['HasPeroxide'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]OO[#6]'))
        
        # 计算毒性风险评分
        toxicity_flags = [
            feat['HasNitro'],
            feat['HasAzo'],
            feat['HasHalogenatedCarbon'],
            feat['HasQuinone'],
            feat['HasMichaelAcceptor'],
            feat['HasPeroxide']
        ]
        feat['ToxicityRiskScore'] = sum(toxicity_flags)
        
        features.append(feat)
    
    return pd.DataFrame(features)

def compare_distributions(df1, df2, feature, label1, label2):
    """对比两组的分布"""
    if feature not in df1.columns or feature not in df2.columns:
        return None
    
    vals1 = df1[feature].dropna()
    vals2 = df2[feature].dropna()
    
    if len(vals1) == 0 or len(vals2) == 0:
        return None
    
    stats = {
        'Feature': feature,
        f'{label1}_Mean': vals1.mean(),
        f'{label1}_Std': vals1.std(),
        f'{label1}_Min': vals1.min(),
        f'{label1}_Max': vals1.max(),
        f'{label2}_Mean': vals2.mean(),
        f'{label2}_Std': vals2.std(),
        f'{label2}_Min': vals2.min(),
        f'{label2}_Max': vals2.max(),
        'Diff_Mean': vals2.mean() - vals1.mean(),
        'Diff_Pct': ((vals2.mean() - vals1.mean()) / vals1.mean() * 100) if vals1.mean() != 0 else 0
    }
    
    return stats

def main():
    print("="*80)
    print("深度对比分析: run14a vs run15")
    print("="*80)
    
    # 加载数据
    df14a = load_gold_candidates("experiments/runs/run14a")
    df15 = load_gold_candidates("experiments/runs/run15")
    
    if df14a is None or df15 is None:
        print("❌ 数据加载失败")
        return
    
    print(f"\n📊 数据量:")
    print(f"  run14a: {len(df14a)} 个gold候选物")
    print(f"  run15:  {len(df15)} 个gold候选物 ({len(df15)-len(df14a):+d})")
    
    # 计算分子特征
    print("\n🧪 计算分子特征...")
    feat14a = calculate_molecular_features(df14a)
    feat15 = calculate_molecular_features(df15)
    
    print(f"  run14a: {len(feat14a)} 个有效分子")
    print(f"  run15:  {len(feat15)} 个有效分子")
    
    # 对比分析
    print("\n" + "="*80)
    print("📈 特征对比分析")
    print("="*80)
    
    features_to_compare = [
        'pIC50', 'QED', 'SA', 'MW', 'LogP', 'TPSA',
        'HBA', 'HBD', 'RotBonds', 'AromaticRings', 
        'HeavyAtoms', 'FractionCsp3', 'NumRings', 'ToxicityRiskScore'
    ]
    
    comparison_results = []
    for feature in features_to_compare:
        stats = compare_distributions(feat14a, feat15, feature, 'run14a', 'run15')
        if stats:
            comparison_results.append(stats)
    
    comp_df = pd.DataFrame(comparison_results)
    
    # 打印关键差异
    print("\n🔍 关键差异 (Mean ± Std):")
    print("-"*80)
    
    for _, row in comp_df.iterrows():
        feature = row['Feature']
        mean14a = row['run14a_Mean']
        std14a = row['run14a_Std']
        mean15 = row['run15_Mean']
        std15 = row['run15_Std']
        diff = row['Diff_Mean']
        diff_pct = row['Diff_Pct']
        
        arrow = "↑" if diff > 0 else "↓" if diff < 0 else "="
        
        print(f"{feature:20s}  run14a: {mean14a:6.2f}±{std14a:5.2f}  "
              f"run15: {mean15:6.2f}±{std15:5.2f}  "
              f"{arrow} {diff:+7.2f} ({diff_pct:+6.1f}%)")
    
    # 毒性特征分析
    print("\n" + "="*80)
    print("☠️  毒性相关特征分析")
    print("="*80)
    
    tox_features = ['HasNitro', 'HasAzo', 'HasHalogenatedCarbon', 
                    'HasQuinone', 'HasMichaelAcceptor', 'HasPeroxide']
    
    print("\n含有潜在毒性基团的分子比例:")
    print("-"*80)
    
    for tox_feat in tox_features:
        if tox_feat in feat14a.columns and tox_feat in feat15.columns:
            pct14a = feat14a[tox_feat].mean() * 100
            pct15 = feat15[tox_feat].mean() * 100
            diff = pct15 - pct14a
            
            arrow = "↓" if diff < 0 else "↑" if diff > 0 else "="
            
            print(f"{tox_feat:25s}  run14a: {pct14a:5.1f}%  "
                  f"run15: {pct15:5.1f}%  {arrow} {diff:+6.1f}%")
    
    # 毒性风险评分
    print("\n毒性风险评分分布:")
    print("-"*80)
    
    risk14a = feat14a['ToxicityRiskScore'].value_counts().sort_index()
    risk15 = feat15['ToxicityRiskScore'].value_counts().sort_index()
    
    all_scores = sorted(set(risk14a.index) | set(risk15.index))
    
    for score in all_scores:
        count14a = risk14a.get(score, 0)
        count15 = risk15.get(score, 0)
        pct14a = count14a / len(feat14a) * 100
        pct15 = count15 / len(feat15) * 100
        
        print(f"风险评分 {score}:  run14a: {count14a:3d} ({pct14a:5.1f}%)  "
              f"run15: {count15:3d} ({pct15:5.1f}%)")
    
    # 骨架多样性
    print("\n" + "="*80)
    print("🧩 骨架多样性分析")
    print("="*80)
    
    scaffolds14a = set(feat14a['Scaffold'].dropna())
    scaffolds15 = set(feat15['Scaffold'].dropna())
    
    common_scaffolds = scaffolds14a & scaffolds15
    unique14a = scaffolds14a - scaffolds15
    unique15 = scaffolds15 - scaffolds14a
    
    print(f"\nrun14a独特骨架: {len(scaffolds14a)} 个")
    print(f"run15独特骨架:  {len(scaffolds15)} 个")
    print(f"共享骨架:       {len(common_scaffolds)} 个")
    print(f"run14a特有:     {len(unique14a)} 个")
    print(f"run15特有:      {len(unique15)} 个")
    
    # 最常见骨架
    print("\nrun14a最常见骨架 (Top 5):")
    scaffold_counts14a = feat14a['Scaffold'].value_counts().head(5)
    for scaffold, count in scaffold_counts14a.items():
        print(f"  {count:3d}x  {scaffold[:60]}...")
    
    print("\nrun15最常见骨架 (Top 5):")
    scaffold_counts15 = feat15['Scaffold'].value_counts().head(5)
    for scaffold, count in scaffold_counts15.items():
        print(f"  {count:3d}x  {scaffold[:60]}...")
    
    # 保存详细对比结果
    print("\n" + "="*80)
    print("💾 保存分析结果")
    print("="*80)
    
    # 保存对比统计
    comp_df.to_csv('run14a_vs_run15_comparison.csv', index=False)
    print("✓ run14a_vs_run15_comparison.csv")
    
    # 保存分子特征
    feat14a.to_csv('run14a_gold_features.csv', index=False)
    feat15.to_csv('run15_gold_features.csv', index=False)
    print("✓ run14a_gold_features.csv")
    print("✓ run15_gold_features.csv")
    
    # 生成可视化
    print("\n📊 生成可视化图表...")
    
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle('run14a vs run15: 金标准分子特征对比', fontsize=16, fontweight='bold')
    
    plot_features = [
        ('pIC50', 'Predicted Activity'),
        ('QED', 'Drug-likeness'),
        ('SA', 'Synthetic Accessibility'),
        ('MW', 'Molecular Weight'),
        ('LogP', 'Lipophilicity'),
        ('TPSA', 'Polar Surface Area'),
        ('RotBonds', 'Rotatable Bonds'),
        ('AromaticRings', 'Aromatic Rings'),
        ('ToxicityRiskScore', 'Toxicity Risk Score')
    ]
    
    for idx, (feature, title) in enumerate(plot_features):
        ax = axes[idx // 3, idx % 3]
        
        if feature in feat14a.columns and feature in feat15.columns:
            data14a = feat14a[feature].dropna()
            data15 = feat15[feature].dropna()
            
            if len(data14a) > 0 and len(data15) > 0:
                ax.hist(data14a, bins=20, alpha=0.6, label='run14a', color='blue', edgecolor='black')
                ax.hist(data15, bins=20, alpha=0.6, label='run15', color='red', edgecolor='black')
                
                ax.axvline(data14a.mean(), color='blue', linestyle='--', linewidth=2, label=f'run14a mean: {data14a.mean():.2f}')
                ax.axvline(data15.mean(), color='red', linestyle='--', linewidth=2, label=f'run15 mean: {data15.mean():.2f}')
                
                ax.set_xlabel(feature)
                ax.set_ylabel('Frequency')
                ax.set_title(title)
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('run14a_vs_run15_distributions.png', dpi=300, bbox_inches='tight')
    print("✓ run14a_vs_run15_distributions.png")
    
    # 生成总结报告
    print("\n" + "="*80)
    print("📋 总结报告")
    print("="*80)
    
    print("\n🎯 主要发现:")
    
    # 活性差异
    pic50_diff = comp_df[comp_df['Feature'] == 'pIC50']['Diff_Mean'].values[0] if 'pIC50' in comp_df['Feature'].values else 0
    if abs(pic50_diff) > 0.05:
        print(f"\n1. 活性变化: run15的平均pIC50 {pic50_diff:+.2f}")
        if pic50_diff > 0:
            print("   ✓ 毒性组件帮助提高了活性预测")
        else:
            print("   ✗ 毒性组件可能过度限制了活性")
    
    # 毒性风险
    avg_risk14a = feat14a['ToxicityRiskScore'].mean()
    avg_risk15 = feat15['ToxicityRiskScore'].mean()
    print(f"\n2. 毒性风险: run14a平均{avg_risk14a:.2f}, run15平均{avg_risk15:.2f}")
    if avg_risk15 < avg_risk14a:
        print(f"   ✓ 毒性组件成功降低了毒性风险 ({(avg_risk14a-avg_risk15)/avg_risk14a*100:.1f}%)")
    else:
        print(f"   ⚠️  毒性风险未明显降低")
    
    # 骨架多样性
    diversity_ratio14a = len(scaffolds14a) / len(feat14a)
    diversity_ratio15 = len(scaffolds15) / len(feat15)
    print(f"\n3. 骨架多样性: run14a={diversity_ratio14a:.2f}, run15={diversity_ratio15:.2f}")
    if diversity_ratio15 > diversity_ratio14a:
        print(f"   ✓ run15的骨架更多样化")
    else:
        print(f"   ⚠️  run15的骨架多样性降低")
    
    # 类药性
    qed_diff = comp_df[comp_df['Feature'] == 'QED']['Diff_Mean'].values[0] if 'QED' in comp_df['Feature'].values else 0
    sa_diff = comp_df[comp_df['Feature'] == 'SA']['Diff_Mean'].values[0] if 'SA' in comp_df['Feature'].values else 0
    print(f"\n4. 类药性: QED变化{qed_diff:+.3f}, SA变化{sa_diff:+.2f}")
    
    print("\n" + "="*80)
    print("✅ 分析完成!")
    print("="*80)
    print("\n生成文件:")
    print("  • run14a_vs_run15_comparison.csv - 详细统计对比")
    print("  • run14a_gold_features.csv - run14a分子特征")
    print("  • run15_gold_features.csv - run15分子特征")
    print("  • run14a_vs_run15_distributions.png - 特征分布图")
    print("="*80)

if __name__ == "__main__":
    main()
