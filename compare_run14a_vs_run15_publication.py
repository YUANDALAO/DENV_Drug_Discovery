#!/usr/bin/env python3
"""
深度对比分析: run14a vs run15
Publication-ready版本 - 顶刊配色 + 修复字体
"""
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# 设置专业出版级别的matplotlib参数
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.figsize'] = (18, 12)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 10
plt.rcParams['ytick.labelsize'] = 10
plt.rcParams['legend.fontsize'] = 9

# Nature/Science期刊配色方案
COLOR_RUN14A = '#0173B2'  # 深蓝色 - run14a (无毒性组件)
COLOR_RUN15 = '#DE8F05'   # 橙色 - run15 (含毒性组件)
COLOR_EDGE = '#333333'    # 深灰色 - 边框

def load_gold_candidates(run_dir):
    """加载金标准候选物"""
    import glob
    
    gold_files = glob.glob(f"{run_dir}/candidates_gold*.csv")
    if not gold_files:
        gold_files = glob.glob(f"{run_dir}/*gold*.csv")
    if not gold_files:
        gold_files = glob.glob(f"{run_dir}/promising*.csv")
    
    if gold_files:
        df = pd.read_csv(gold_files[0])
        print(f"✓ Loaded {run_dir}: {len(df)} gold candidates")
        return df
    else:
        print(f"✗ Not found: {run_dir}")
        return None

def calculate_molecular_features(df):
    """计算分子特征"""
    features = []
    
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            continue
        
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
            'FractionCsp3': Lipinski.FractionCSP3(mol),
            'NumRings': Descriptors.RingCount(mol),
            'NumAliphaticRings': Descriptors.NumAliphaticRings(mol),
        }
        
        for col_name, feat_name in [
            ('DENV_Activity (raw)', 'pIC50'),
            ('DENV_Activity', 'pIC50'),
            ('QED (raw)', 'QED'),
            ('QED', 'QED'),
            ('SA (raw)', 'SA'),
            ('SA', 'SA'),
            ('Score', 'TotalScore'),
        ]:
            if col_name in df.columns and feat_name not in feat:
                feat[feat_name] = row[col_name]
        
        try:
            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
            feat['Scaffold'] = Chem.MolToSmiles(scaffold)
        except:
            feat['Scaffold'] = None
        
        toxicity_substructures = {
            'HasNitro': '[N+](=O)[O-]',
            'HasAzo': 'N=N',
            'HasHalogenatedCarbon': '[F,Cl,Br,I][C,c][F,Cl,Br,I]',
            'HasQuinone': 'C1=CC(=O)C=CC1=O',
            'HasMichaelAcceptor': '[C]=[C]-[C]=[O,S]',
            'HasPeroxide': '[#6]OO[#6]',
            'HasAcylChloride': 'C(=O)Cl',
            'HasSulfonylChloride': 'S(=O)(=O)Cl',
        }
        
        toxicity_flags = []
        for name, smarts in toxicity_substructures.items():
            try:
                pattern = Chem.MolFromSmarts(smarts)
                has_match = mol.HasSubstructMatch(pattern) if pattern else False
                feat[name] = has_match
                toxicity_flags.append(has_match)
            except:
                feat[name] = False
        
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
    print("Publication-Ready Analysis: run14a vs run15")
    print("="*80)
    
    # 加载数据
    df14a = load_gold_candidates("experiments/runs/run14a")
    df15 = load_gold_candidates("experiments/runs/run15")
    
    if df14a is None or df15 is None:
        print("Error: Failed to load data")
        return
    
    print(f"\nData Summary:")
    print(f"  run14a: {len(df14a)} gold candidates")
    print(f"  run15:  {len(df15)} gold candidates ({len(df15)-len(df14a):+d})")
    
    # 计算分子特征
    print("\nCalculating molecular features...")
    feat14a = calculate_molecular_features(df14a)
    feat15 = calculate_molecular_features(df15)
    
    print(f"  run14a: {len(feat14a)} valid molecules")
    print(f"  run15:  {len(feat15)} valid molecules")
    
    # 对比分析
    print("\n" + "="*80)
    print("Feature Comparison")
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
    print("\nKey Differences (Mean ± SD):")
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
    print("Toxicity Features Analysis")
    print("="*80)
    
    tox_features = ['HasNitro', 'HasAzo', 'HasHalogenatedCarbon', 
                    'HasQuinone', 'HasMichaelAcceptor', 'HasPeroxide',
                    'HasAcylChloride', 'HasSulfonylChloride']
    
    print("\nMolecules with Toxicity Substructures:")
    print("-"*80)
    
    for tox_feat in tox_features:
        if tox_feat in feat14a.columns and tox_feat in feat15.columns:
            pct14a = feat14a[tox_feat].mean() * 100
            pct15 = feat15[tox_feat].mean() * 100
            diff = pct15 - pct14a
            
            arrow = "✓" if diff < 0 else "✗" if diff > 0 else "="
            
            print(f"{tox_feat:25s}  run14a: {pct14a:5.1f}%  "
                  f"run15: {pct15:5.1f}%  {arrow} {diff:+6.1f}%")
    
    # 骨架多样性
    print("\n" + "="*80)
    print("Scaffold Diversity")
    print("="*80)
    
    scaffolds14a = set(feat14a['Scaffold'].dropna())
    scaffolds15 = set(feat15['Scaffold'].dropna())
    
    diversity14a = len(scaffolds14a) / len(feat14a)
    diversity15 = len(scaffolds15) / len(feat15)
    
    print(f"\nrun14a: {len(scaffolds14a)} unique scaffolds (diversity ratio: {diversity14a:.3f})")
    print(f"run15:  {len(scaffolds15)} unique scaffolds (diversity ratio: {diversity15:.3f})")
    
    # 保存结果
    print("\n" + "="*80)
    print("Saving Results")
    print("="*80)
    
    comp_df.to_csv('run14a_vs_run15_comparison.csv', index=False)
    feat14a.to_csv('run14a_gold_features.csv', index=False)
    feat15.to_csv('run15_gold_features.csv', index=False)
    print("✓ CSV files saved")
    
    # 生成Publication-ready可视化
    print("\nGenerating publication-ready figures...")
    
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3, 
                          left=0.08, right=0.98, top=0.94, bottom=0.06)
    
    # 标题
    fig.suptitle('Comparison of run14a vs run15: Impact of Toxicity Filters\n' + 
                 'run15 includes: Toxicity_Alerts, PAINS_Filter, Metabolic_Stability',
                 fontsize=14, fontweight='bold', y=0.98)
    
    plot_features = [
        ('pIC50', 'Predicted Activity (pIC50)', 'higher is better'),
        ('QED', 'Drug-likeness (QED)', 'higher is better'),
        ('SA', 'Synthetic Accessibility', 'lower is better'),
        ('MW', 'Molecular Weight (Da)', ''),
        ('LogP', 'Lipophilicity (LogP)', ''),
        ('TPSA', 'Polar Surface Area (Å²)', ''),
        ('RotBonds', 'Rotatable Bonds', ''),
        ('AromaticRings', 'Aromatic Rings', ''),
        ('ToxicityRiskScore', 'Toxicity Risk Score', 'lower is better')
    ]
    
    for idx, (feature, title, note) in enumerate(plot_features):
        ax = fig.add_subplot(gs[idx // 3, idx % 3])
        
        if feature in feat14a.columns and feature in feat15.columns:
            data14a = feat14a[feature].dropna()
            data15 = feat15[feature].dropna()
            
            if len(data14a) > 0 and len(data15) > 0:
                # 绘制直方图
                bins = np.linspace(
                    min(data14a.min(), data15.min()),
                    max(data14a.max(), data15.max()),
                    20
                )
                
                ax.hist(data14a, bins=bins, alpha=0.7, label='run14a (no toxicity filters)', 
                       color=COLOR_RUN14A, edgecolor=COLOR_EDGE, linewidth=0.5)
                ax.hist(data15, bins=bins, alpha=0.7, label='run15 (with toxicity filters)', 
                       color=COLOR_RUN15, edgecolor=COLOR_EDGE, linewidth=0.5)
                
                # 添加均值线
                mean14a = data14a.mean()
                mean15 = data15.mean()
                
                ax.axvline(mean14a, color=COLOR_RUN14A, linestyle='--', linewidth=2, alpha=0.8)
                ax.axvline(mean15, color=COLOR_RUN15, linestyle='--', linewidth=2, alpha=0.8)
                
                # 添加均值标注
                y_pos = ax.get_ylim()[1] * 0.95
                ax.text(mean14a, y_pos, f'{mean14a:.2f}', 
                       ha='center', va='top', fontsize=9, fontweight='bold',
                       color=COLOR_RUN14A,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor=COLOR_RUN14A, alpha=0.8))
                ax.text(mean15, y_pos * 0.82, f'{mean15:.2f}', 
                       ha='center', va='top', fontsize=9, fontweight='bold',
                       color=COLOR_RUN15,
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                edgecolor=COLOR_RUN15, alpha=0.8))
                
                # 设置标签
                ax.set_xlabel(feature, fontsize=11, fontweight='bold')
                ax.set_ylabel('Frequency', fontsize=11)
                ax.set_title(title + (f'\n({note})' if note else ''), 
                           fontsize=11, fontweight='bold', pad=10)
                
                # 图例
                if idx == 0:
                    ax.legend(loc='upper right', frameon=True, fancybox=False, 
                            shadow=False, framealpha=0.9, edgecolor='black')
                
                # 网格
                ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
                ax.set_axisbelow(True)
                
                # 边框
                for spine in ax.spines.values():
                    spine.set_edgecolor(COLOR_EDGE)
                    spine.set_linewidth(1.2)
    
    # 保存高质量图片
    for fmt in ['png', 'pdf']:
        filename = f'run14a_vs_run15_publication.{fmt}'
        plt.savefig(filename, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        print(f"✓ {filename}")
    
    plt.close()
    
    # 生成总结
    print("\n" + "="*80)
    print("Summary")
    print("="*80)
    
    print("\nKey Findings:")
    
    if 'pIC50' in comp_df['Feature'].values:
        pic50_diff = comp_df[comp_df['Feature'] == 'pIC50']['Diff_Mean'].values[0]
        pic50_pct = comp_df[comp_df['Feature'] == 'pIC50']['Diff_Pct'].values[0]
        print(f"\n1. Activity Change:")
        print(f"   run15 pIC50: {pic50_diff:+.3f} ({pic50_pct:+.1f}%)")
        if abs(pic50_diff) < 0.05:
            print("   → Toxicity filters maintained activity levels")
        elif pic50_diff > 0:
            print("   ✓ Toxicity filters slightly improved activity")
        else:
            print("   → Toxicity filters slightly reduced activity")
    
    avg_risk14a = feat14a['ToxicityRiskScore'].mean()
    avg_risk15 = feat15['ToxicityRiskScore'].mean()
    risk_reduction = (avg_risk14a - avg_risk15) / avg_risk14a * 100 if avg_risk14a > 0 else 0
    
    print(f"\n2. Toxicity Risk:")
    print(f"   run14a: {avg_risk14a:.3f}")
    print(f"   run15:  {avg_risk15:.3f}")
    if risk_reduction > 5:
        print(f"   ✓ Significant toxicity reduction ({risk_reduction:.1f}%)")
    elif risk_reduction > 0:
        print(f"   ✓ Modest toxicity reduction ({risk_reduction:.1f}%)")
    else:
        print(f"   ⚠  No significant toxicity improvement")
    
    print(f"\n3. Scaffold Diversity:")
    print(f"   run14a: {diversity14a:.3f}")
    print(f"   run15:  {diversity15:.3f}")
    if diversity15 > diversity14a * 1.05:
        print(f"   ✓ Increased diversity (+{(diversity15/diversity14a-1)*100:.1f}%)")
    elif diversity15 < diversity14a * 0.95:
        print(f"   ⚠  Reduced diversity ({(diversity15/diversity14a-1)*100:.1f}%)")
    else:
        print(f"   → Similar diversity")
    
    print("\n" + "="*80)
    print("✅ Analysis Complete!")
    print("="*80)
    print("\nGenerated Files:")
    print("  • run14a_vs_run15_comparison.csv")
    print("  • run14a_gold_features.csv") 
    print("  • run15_gold_features.csv")
    print("  • run14a_vs_run15_publication.png (300 DPI)")
    print("  • run14a_vs_run15_publication.pdf (vector)")
    print("="*80)

if __name__ == "__main__":
    main()
