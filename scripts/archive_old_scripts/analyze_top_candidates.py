"""
深度分析Top候选物
提供多个筛选标准和可视化
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

def analyze_candidates(results_csv, output_dir):
    """分析Top候选物"""
    
    df = pd.read_csv(results_csv)
    
    # 计算IC50
    if 'DENV_Activity (raw)' in df.columns:
        df['IC50_nM'] = 10 ** (9 - df['DENV_Activity (raw)'])
    
    print("="*60)
    print("多层次候选物筛选")
    print("="*60)
    
    # 不同严格程度的筛选
    criteria_sets = [
        {
            'name': '金标准 (Gold Standard)',
            'pIC50': 8.0,
            'QED': 0.7,
            'SA': 4.0,
            'MW_min': 300,
            'MW_max': 500,
            'LogP_min': 1,
            'LogP_max': 4
        },
        {
            'name': '高标准 (High Quality)',
            'pIC50': 7.5,
            'QED': 0.6,
            'SA': 4.5,
            'MW_min': 250,
            'MW_max': 550,
            'LogP_min': 0.5,
            'LogP_max': 5
        },
        {
            'name': '中标准 (Good Quality)',
            'pIC50': 7.0,
            'QED': 0.5,
            'SA': 5.0,
            'MW_min': 200,
            'MW_max': 600,
            'LogP_min': 0,
            'LogP_max': 6
        }
    ]
    
    for criteria in criteria_sets:
        print(f"\n{'='*60}")
        print(f"{criteria['name']}")
        print(f"{'='*60}")
        
        mask = (
            (df['DENV_Activity (raw)'] >= criteria['pIC50']) &
            (df['QED (raw)'] >= criteria['QED']) &
            (df['SA (raw)'] <= criteria['SA']) &
            (df['MW (raw)'] >= criteria['MW_min']) &
            (df['MW (raw)'] <= criteria['MW_max']) &
            (df['LogP (raw)'] >= criteria['LogP_min']) &
            (df['LogP (raw)'] <= criteria['LogP_max'])
        )
        
        candidates = df[mask].copy()
        candidates = candidates.sort_values('Score', ascending=False)
        
        print(f"满足条件的分子数: {len(candidates)} ({len(candidates)/len(df)*100:.1f}%)")
        
        if len(candidates) > 0:
            print(f"\n统计信息:")
            print(f"  pIC50: {candidates['DENV_Activity (raw)'].min():.2f} - {candidates['DENV_Activity (raw)'].max():.2f}")
            print(f"  IC50: {candidates['IC50_nM'].min():.1f} - {candidates['IC50_nM'].max():.1f} nM")
            print(f"  QED: {candidates['QED (raw)'].min():.2f} - {candidates['QED (raw)'].max():.2f}")
            print(f"  SA Score: {candidates['SA (raw)'].min():.2f} - {candidates['SA (raw)'].max():.2f}")
            print(f"  总分: {candidates['Score'].min():.3f} - {candidates['Score'].max():.3f}")
            
            # 保存
            output_file = output_dir / f"candidates_{criteria['name'].replace(' ', '_').replace('(', '').replace(')', '')}.csv"
            candidates.to_csv(output_file, index=False)
            print(f"\n✓ 已保存: {output_file}")
            
            # 显示Top 5
            print(f"\nTop 5 分子:")
            display_cols = ['SMILES', 'Score', 'DENV_Activity (raw)', 'IC50_nM', 'QED (raw)', 'SA (raw)', 'MW (raw)', 'LogP (raw)']
            for i, row in candidates.head(5).iterrows():
                print(f"\n{i+1}. Score: {row['Score']:.3f}, IC50: {row['IC50_nM']:.1f} nM")
                print(f"   SMILES: {row['SMILES']}")
                print(f"   pIC50={row['DENV_Activity (raw)']:.2f}, QED={row['QED (raw)']:.2f}, SA={row['SA (raw)']:.2f}")
    
    # 按不同属性排序的Top分子
    print(f"\n{'='*60}")
    print("各属性Top 10分子")
    print(f"{'='*60}")
    
    rankings = [
        ('总分', 'Score', False),
        ('活性 (pIC50)', 'DENV_Activity (raw)', False),
        ('类药性 (QED)', 'QED (raw)', False),
        ('可合成性 (SA)', 'SA (raw)', True),
    ]
    
    for name, col, ascending in rankings:
        print(f"\n{name} Top 10:")
        top = df.nlargest(10, col) if not ascending else df.nsmallest(10, col)
        
        output_file = output_dir / f"top10_by_{col.replace(' ', '_').replace('(', '').replace(')', '')}.csv"
        top.to_csv(output_file, index=False)
        
        for idx, (i, row) in enumerate(top.iterrows(), 1):
            print(f"  {idx}. {row['SMILES'][:60]}")
            print(f"     {name}={row[col]:.3f}, Score={row['Score']:.3f}, IC50={row['IC50_nM']:.1f}nM")
    
    # 综合分析
    print(f"\n{'='*60}")
    print("综合统计")
    print(f"{'='*60}")
    
    # IC50分布
    ic50_ranges = [
        (0, 10, '0-10 nM (极高活性)'),
        (10, 50, '10-50 nM (高活性)'),
        (50, 100, '50-100 nM (中等活性)'),
        (100, 1000, '100-1000 nM (低活性)'),
    ]
    
    print("\nIC50分布:")
    for low, high, label in ic50_ranges:
        count = len(df[(df['IC50_nM'] >= low) & (df['IC50_nM'] < high)])
        pct = count / len(df) * 100
        print(f"  {label}: {count} ({pct:.1f}%)")
    
    # 多维度优秀分子（放宽标准）
    print(f"\n多维度优秀分子 (pIC50>7.5 OR QED>0.75 OR SA<3.5):")
    multi_good = df[
        (df['DENV_Activity (raw)'] > 7.5) |
        (df['QED (raw)'] > 0.75) |
        (df['SA (raw)'] < 3.5)
    ]
    print(f"  符合至少一个优秀标准: {len(multi_good)} ({len(multi_good)/len(df)*100:.1f}%)")
    
    # 同时满足2个
    multi_good_2 = df[
        ((df['DENV_Activity (raw)'] > 7.5) & (df['QED (raw)'] > 0.65)) |
        ((df['DENV_Activity (raw)'] > 7.5) & (df['SA (raw)'] < 4.0)) |
        ((df['QED (raw)'] > 0.65) & (df['SA (raw)'] < 4.0))
    ]
    print(f"  符合至少两个标准: {len(multi_good_2)} ({len(multi_good_2)/len(df)*100:.1f}%)")
    
    multi_good_2 = multi_good_2.sort_values('Score', ascending=False)
    output_file = output_dir / "candidates_multi_criteria.csv"
    multi_good_2.to_csv(output_file, index=False)
    print(f"  ✓ 已保存: {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python analyze_top_candidates.py <实验文件夹>")
        sys.exit(1)
    
    run_dir = Path(sys.argv[1])
    results_file = run_dir / "results_1.csv"
    
    if not results_file.exists():
        results_file = run_dir / "rescored_results.csv"
    
    if not results_file.exists():
        print(f"错误: 找不到结果文件")
        sys.exit(1)
    
    print(f"分析文件: {results_file}\n")
    analyze_candidates(results_file, run_dir)
    
    print(f"\n{'='*60}")
    print("分析完成！")
    print(f"{'='*60}")