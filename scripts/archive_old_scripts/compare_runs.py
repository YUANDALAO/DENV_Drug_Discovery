"""
对比两次实验的详细分析
找出性能下降的原因
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def analyze_run(run_dir, run_name):
    """分析单次实验"""
    
    results_file = Path(run_dir) / "results_1.csv"
    if not results_file.exists():
        results_file = Path(run_dir) / "rescored_results.csv"
    
    if not results_file.exists():
        print(f"⚠️  找不到 {run_name} 的结果文件")
        return None
    
    df = pd.read_csv(results_file)
    
    # 计算IC50
    if 'DENV_Activity (raw)' in df.columns:
        df['IC50_nM'] = 10 ** (9 - df['DENV_Activity (raw)'])
        df['pIC50'] = df['DENV_Activity (raw)']
    
    # 提取scaffolds
    scaffolds = []
    for smiles in df['SMILES']:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(scaffold)
                scaffolds.append(scaffold_smiles)
            else:
                scaffolds.append(None)
        except:
            scaffolds.append(None)
    
    df['Scaffold'] = scaffolds
    
    stats = {
        'name': run_name,
        'total_molecules': len(df),
        'unique_scaffolds': df['Scaffold'].nunique(),
        'avg_ic50': df['IC50_nM'].mean() if 'IC50_nM' in df.columns else None,
        'min_ic50': df['IC50_nM'].min() if 'IC50_nM' in df.columns else None,
        'max_ic50': df['IC50_nM'].max() if 'IC50_nM' in df.columns else None,
        'molecules_per_scaffold': len(df) / df['Scaffold'].nunique(),
    }
    
    # 金标准分子
    if all(col in df.columns for col in ['pIC50', 'QED (raw)', 'SA (raw)', 'MW (raw)', 'LogP (raw)']):
        gold = df[
            (df['pIC50'] >= 8.0) &
            (df['QED (raw)'] >= 0.7) &
            (df['SA (raw)'] <= 4.0) &
            (df['MW (raw)'] >= 300) &
            (df['MW (raw)'] <= 500) &
            (df['LogP (raw)'] >= 1) &
            (df['LogP (raw)'] <= 4)
        ]
        stats['gold_standard'] = len(gold)
    
    # 高标准分子
        high = df[
            (df['pIC50'] >= 7.5) &
            (df['QED (raw)'] >= 0.6) &
            (df['SA (raw)'] <= 4.5) &
            (df['MW (raw)'] >= 250) &
            (df['MW (raw)'] <= 550)
        ]
        stats['high_quality'] = len(high)
    
    return df, stats


def compare_runs(run7_dir, run11_dir, output_dir):
    """对比两次实验"""
    
    print("="*60)
    print("对比分析: Run7 vs Run11")
    print("="*60)
    
    df7, stats7 = analyze_run(run7_dir, "Run7 (1500步)")
    df11, stats11 = analyze_run(run11_dir, "Run11 (4000步)")
    
    if df7 is None or df11 is None:
        return
    
    # 打印统计对比
    print("\n基本统计:")
    print(f"{'指标':<30} {'Run7':<20} {'Run11':<20} {'变化':<15}")
    print("-"*85)
    
    metrics = [
        ('总分子数', 'total_molecules', 'int'),
        ('唯一Scaffolds', 'unique_scaffolds', 'int'),
        ('平均IC50 (nM)', 'avg_ic50', 'float'),
        ('最佳IC50 (nM)', 'min_ic50', 'float'),
        ('每个Scaffold分子数', 'molecules_per_scaffold', 'float'),
        ('金标准候选物', 'gold_standard', 'int'),
        ('高质量候选物', 'high_quality', 'int'),
    ]
    
    for name, key, dtype in metrics:
        val7 = stats7.get(key, 'N/A')
        val11 = stats11.get(key, 'N/A')
        
        if val7 != 'N/A' and val11 != 'N/A':
            if dtype == 'int':
                change = f"{val11 - val7:+d}"
                val7_str = f"{val7:d}"
                val11_str = f"{val11:d}"
            else:
                change = f"{val11 - val7:+.2f}"
                val7_str = f"{val7:.2f}"
                val11_str = f"{val11:.2f}"
        else:
            change = "N/A"
            val7_str = str(val7)
            val11_str = str(val11)
        
        print(f"{name:<30} {val7_str:<20} {val11_str:<20} {change:<15}")
    
    # IC50分布对比
    print("\n\nIC50分布对比:")
    print(f"{'IC50范围':<30} {'Run7':<20} {'Run11':<20}")
    print("-"*70)
    
    ranges = [
        (0, 10, '0-10 nM (极高)'),
        (10, 20, '10-20 nM (很高)'),
        (20, 30, '20-30 nM (高)'),
        (30, 50, '30-50 nM (中高)'),
    ]
    
    for low, high, label in ranges:
        if 'IC50_nM' in df7.columns and 'IC50_nM' in df11.columns:
            count7 = len(df7[(df7['IC50_nM'] >= low) & (df7['IC50_nM'] < high)])
            count11 = len(df11[(df11['IC50_nM'] >= low) & (df11['IC50_nM'] < high)])
            pct7 = count7 / len(df7) * 100
            pct11 = count11 / len(df11) * 100
            print(f"{label:<30} {count7:5d} ({pct7:5.1f}%)      {count11:5d} ({pct11:5.1f}%)")
    
    # Scaffold多样性分析
    print("\n\nScaffold多样性分析:")
    
    scaffold_counts7 = df7['Scaffold'].value_counts()
    scaffold_counts11 = df11['Scaffold'].value_counts()
    
    print(f"\nRun7 Top 5 Scaffolds:")
    for i, (scaffold, count) in enumerate(scaffold_counts7.head(5).items(), 1):
        print(f"  {i}. {scaffold[:60]}")
        print(f"     出现次数: {count} ({count/len(df7)*100:.1f}%)")
    
    print(f"\nRun11 Top 5 Scaffolds:")
    for i, (scaffold, count) in enumerate(scaffold_counts11.head(5).items(), 1):
        print(f"  {i}. {scaffold[:60]}")
        print(f"     出现次数: {count} ({count/len(df11)*100:.1f}%)")
    
    # 分析diversity filter影响
    print("\n\n" + "="*60)
    print("Diversity Filter 影响分析")
    print("="*60)
    
    bucket_size = 25
    
    print(f"\n当前设置: bucket_size = {bucket_size}")
    print(f"\nRun7:")
    print(f"  理论最大保留: {stats7['unique_scaffolds']} scaffolds × {bucket_size} = {stats7['unique_scaffolds'] * bucket_size}")
    print(f"  实际保留: {stats7['total_molecules']}")
    print(f"  利用率: {stats7['total_molecules'] / (stats7['unique_scaffolds'] * bucket_size) * 100:.1f}%")
    
    print(f"\nRun11:")
    print(f"  理论最大保留: {stats11['unique_scaffolds']} scaffolds × {bucket_size} = {stats11['unique_scaffolds'] * bucket_size}")
    print(f"  实际保留: {stats11['total_molecules']}")
    print(f"  利用率: {stats11['total_molecules'] / (stats11['unique_scaffolds'] * bucket_size) * 100:.1f}%")
    
    if stats11['total_molecules'] / (stats11['unique_scaffolds'] * bucket_size) > 0.95:
        print("\n⚠️  警告: Run11的diversity filter已接近饱和！")
        print("   建议: 增加 bucket_size 到 50-100")
    
    # 结论和建议
    print("\n\n" + "="*60)
    print("分析结论")
    print("="*60)
    
    if stats11['gold_standard'] < stats7.get('gold_standard', 0):
        print("\n❌ Run11的金标准候选物数量下降")
        print("\n可能原因:")
        print("  1. Diversity filter饱和，后期优质分子被丢弃")
        print("  2. 训练时间过长导致模式崩塌")
        print("  3. Batch size增大导致探索不足")
        
        print("\n建议改进:")
        print("  1. 增加 bucket_size 从 25 到 50-100")
        print("  2. 减少训练步数到 2000-2500")
        print("  3. 保持 batch_size = 64")
        print("  4. 增加 diversity_filter.minscore 到 0.6-0.7")
        print("  5. 使用 curriculum learning (多阶段训练)")


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("用法: python compare_runs.py <run7目录> <run11目录>")
        sys.exit(1)
    
    run7_dir = sys.argv[1]
    run11_dir = sys.argv[2]
    output_dir = Path(run11_dir)
    
    compare_runs(run7_dir, run11_dir, output_dir)
