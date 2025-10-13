"""
LibInvent全量评价 - 分析所有113万分子
包含：并行计算、分批处理、进度保存
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, Lipinski
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import pickle
from multiprocessing import Pool, cpu_count

# ============================================================================
# 配置
# ============================================================================
RESULTS_FILE = "results_1.csv"
OUTPUT_DIR = "full_evaluation"
BATCH_SIZE = 10000  # 分批处理，避免内存爆炸
N_WORKERS = min(8, cpu_count())  # 并行核心数

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 并行计算函数
# ============================================================================

def process_molecule(args):
    """处理单个分子 - 用于并行计算"""
    idx, smi, score = args
    
    if pd.isna(smi):
        return None
    
    smi_str = str(smi).strip()
    
    # 跳过片段
    if '[*]' in smi_str or '|' in smi_str or len(smi_str) > 500:
        return None
    
    mol = Chem.MolFromSmiles(smi_str)
    if not mol:
        return None
    
    try:
        # 计算所有指标
        result = {
            'Index': idx,
            'SMILES': smi_str,
            'Score': score,
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
        if result['MW'] > 500: ro5 += 1
        if result['LogP'] > 5: ro5 += 1
        if result['HBA'] > 10: ro5 += 1
        if result['HBD'] > 5: ro5 += 1
        result['RO5_Violations'] = ro5
        
        return result
    except Exception as e:
        return None


def process_batch_parallel(batch_data):
    """并行处理一批分子"""
    with Pool(N_WORKERS) as pool:
        results = list(tqdm(
            pool.imap(process_molecule, batch_data),
            total=len(batch_data),
            desc="  处理中"
        ))
    
    # 过滤None
    return [r for r in results if r is not None]


# ============================================================================
# 主流程
# ============================================================================

def main():
    
    print("="*80)
    print("LibInvent 全量评价 - 分析全部113万分子")
    print("="*80)
    print(f"并行核心数: {N_WORKERS}")
    print(f"分批大小: {BATCH_SIZE}")
    print()
    
    # 检查是否有中断的进度
    checkpoint_file = f"{OUTPUT_DIR}/checkpoint.pkl"
    if os.path.exists(checkpoint_file):
        print("发现中断的进度，继续处理...")
        with open(checkpoint_file, 'rb') as f:
            checkpoint = pickle.load(f)
        all_results = checkpoint['results']
        start_batch = checkpoint['last_batch'] + 1
    else:
        all_results = []
        start_batch = 0
    
    # ========================================================================
    # Step 1: 读取数据
    # ========================================================================
    print("Step 1/3: 读取数据")
    print("-"*80)
    
    df = pd.read_csv(RESULTS_FILE)
    total_mols = len(df)
    print(f"✓ 总分子数: {total_mols}")
    
    n_batches = (total_mols + BATCH_SIZE - 1) // BATCH_SIZE
    print(f"✓ 分为 {n_batches} 批处理")
    print()
    
    # ========================================================================
    # Step 2: 分批并行处理
    # ========================================================================
    print("Step 2/3: 分批处理所有分子")
    print("-"*80)
    
    for batch_idx in range(start_batch, n_batches):
        start_idx = batch_idx * BATCH_SIZE
        end_idx = min(start_idx + BATCH_SIZE, total_mols)
        
        print(f"\n批次 {batch_idx+1}/{n_batches}: 处理 {start_idx}-{end_idx}")
        
        # 准备批次数据
        batch_df = df.iloc[start_idx:end_idx]
        batch_data = [
            (idx, row['SMILES'], row.get('Score', 0))
            for idx, row in batch_df.iterrows()
        ]
        
        # 并行处理
        batch_results = process_batch_parallel(batch_data)
        all_results.extend(batch_results)
        
        print(f"  ✓ 本批有效: {len(batch_results)}/{len(batch_data)}")
        print(f"  ✓ 累计有效: {len(all_results)}")
        
        # 每10批保存一次checkpoint
        if (batch_idx + 1) % 10 == 0:
            checkpoint = {
                'results': all_results,
                'last_batch': batch_idx
            }
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(checkpoint, f)
            print(f"  ✓ 进度已保存")
    
    # 删除checkpoint
    if os.path.exists(checkpoint_file):
        os.remove(checkpoint_file)
    
    # ========================================================================
    # Step 3: 汇总分析
    # ========================================================================
    print("\nStep 3/3: 生成报告")
    print("-"*80)
    
    df_results = pd.DataFrame(all_results)
    
    print(f"✓ 处理完成")
    print(f"  总分子数: {total_mols}")
    print(f"  有效分子: {len(df_results)} ({len(df_results)/total_mols*100:.1f}%)")
    print(f"  独特SMILES: {df_results['SMILES'].nunique()}")
    
    # 保存完整数据
    output_file = f"{OUTPUT_DIR}/all_molecules_metrics.csv"
    df_results.to_csv(output_file, index=False)
    print(f"\n✓ 完整数据已保存: {output_file}")
    print(f"  文件大小: {os.path.getsize(output_file) / 1024 / 1024:.1f} MB")
    
    # 生成统计报告
    generate_full_report(df_results, total_mols)
    
    # 生成可视化
    generate_plots(df_results)
    
    # 保存高质量分子
    save_high_quality(df_results)
    
    print("\n" + "="*80)
    print("✅ 全量评价完成！")
    print("="*80)


def generate_full_report(df_results, total_mols):
    """生成统计报告"""
    
    report = []
    report.append("="*80)
    report.append("LibInvent 全量评价报告")
    report.append("="*80)
    report.append("")
    
    # 1. 总体统计
    report.append("## 1. 总体统计 (全量分析)")
    report.append(f"   原始生成数: {total_mols}")
    report.append(f"   有效分子数: {len(df_results)}")
    report.append(f"   有效率: {len(df_results)/total_mols*100:.2f}%")
    report.append(f"   独特SMILES: {df_results['SMILES'].nunique()}")
    report.append(f"   独特率: {df_results['SMILES'].nunique()/len(df_results)*100:.2f}%")
    report.append("")
    
    # 2. 成药性
    report.append("## 2. 成药性评估 (全量)")
    report.append(f"   平均QED: {df_results['QED'].mean():.4f} ± {df_results['QED'].std():.4f}")
    
    qed_good = (df_results['QED'] > 0.5).sum()
    qed_exc = (df_results['QED'] > 0.7).sum()
    report.append(f"   QED > 0.5: {qed_good} ({qed_good/len(df_results)*100:.2f}%)")
    report.append(f"   QED > 0.7: {qed_exc} ({qed_exc/len(df_results)*100:.2f}%)")
    report.append("")
    
    report.append("   物化性质:")
    report.append(f"   MW:   {df_results['MW'].mean():.2f} ± {df_results['MW'].std():.2f}")
    report.append(f"   LogP: {df_results['LogP'].mean():.3f} ± {df_results['LogP'].std():.3f}")
    report.append(f"   TPSA: {df_results['TPSA'].mean():.2f} ± {df_results['TPSA'].std():.2f}")
    report.append("")
    
    report.append("   Lipinski Rule of 5:")
    for v in range(5):
        cnt = (df_results['RO5_Violations'] == v).sum()
        pct = cnt / len(df_results) * 100
        report.append(f"   {v}违反: {cnt} ({pct:.2f}%)")
    report.append("")
    
    # 3. 高质量分子
    high_quality = df_results[
        (df_results['QED'] > 0.6) &
        (df_results['RO5_Violations'] <= 1) &
        (df_results['MW'].between(200, 600))
    ]
    
    report.append("## 3. 高质量分子 (全量)")
    report.append(f"   数量: {len(high_quality)} ({len(high_quality)/len(df_results)*100:.2f}%)")
    report.append(f"   定义: QED>0.6, RO5≤1, 200<MW<600")
    report.append("")
    
    # 4. 分位数分析
    report.append("## 4. 分位数分析")
    for metric in ['QED', 'MW', 'LogP', 'TPSA']:
        quantiles = df_results[metric].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
        report.append(f"   {metric}:")
        report.append(f"     10%: {quantiles[0.1]:.2f}")
        report.append(f"     25%: {quantiles[0.25]:.2f}")
        report.append(f"     50%: {quantiles[0.5]:.2f}")
        report.append(f"     75%: {quantiles[0.75]:.2f}")
        report.append(f"     90%: {quantiles[0.9]:.2f}")
    report.append("")
    
    # 5. 评级
    qed_good_pct = qed_good / len(df_results) * 100
    if qed_good_pct > 70:
        rating = "优秀 ⭐⭐⭐"
    elif qed_good_pct > 50:
        rating = "良好 ⭐⭐"
    else:
        rating = "一般 ⭐"
    
    report.append(f"## 5. 总体评级: {rating}")
    report.append(f"   基于全部 {len(df_results)} 个有效分子")
    report.append("")
    report.append("="*80)
    
    # 保存报告
    report_text = "\n".join(report)
    with open(f"{OUTPUT_DIR}/full_report.txt", 'w', encoding='utf-8') as f:
        f.write(report_text)
    
    print("\n" + report_text)


def generate_plots(df_results):
    """生成可视化"""
    
    # 采样用于可视化（太多点会很慢）
    df_plot = df_results.sample(n=min(50000, len(df_results)), random_state=42)
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # 1. QED
    axes[0].hist(df_plot['QED'], bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    axes[0].axvline(0.5, color='orange', linestyle='--', linewidth=2)
    axes[0].axvline(0.7, color='green', linestyle='--', linewidth=2)
    axes[0].set_xlabel('QED', fontweight='bold')
    axes[0].set_ylabel('Count', fontweight='bold')
    axes[0].set_title('QED Distribution (Full Dataset)', fontweight='bold')
    axes[0].grid(alpha=0.3)
    
    # 2. MW
    axes[1].hist(df_plot['MW'], bins=50, color='coral', alpha=0.7, edgecolor='black')
    axes[1].axvline(500, color='red', linestyle='--', linewidth=2)
    axes[1].set_xlabel('MW', fontweight='bold')
    axes[1].set_ylabel('Count', fontweight='bold')
    axes[1].set_title('Molecular Weight', fontweight='bold')
    axes[1].grid(alpha=0.3)
    
    # 3. LogP
    axes[2].hist(df_plot['LogP'], bins=50, color='lightgreen', alpha=0.7, edgecolor='black')
    axes[2].axvline(5, color='red', linestyle='--', linewidth=2)
    axes[2].set_xlabel('LogP', fontweight='bold')
    axes[2].set_ylabel('Count', fontweight='bold')
    axes[2].set_title('Lipophilicity', fontweight='bold')
    axes[2].grid(alpha=0.3)
    
    # 4. TPSA
    axes[3].hist(df_plot['TPSA'], bins=50, color='plum', alpha=0.7, edgecolor='black')
    axes[3].set_xlabel('TPSA', fontweight='bold')
    axes[3].set_ylabel('Count', fontweight='bold')
    axes[3].set_title('Polar Surface Area', fontweight='bold')
    axes[3].grid(alpha=0.3)
    
    # 5. RO5
    ro5_counts = df_results['RO5_Violations'].value_counts().sort_index()
    axes[4].bar(ro5_counts.index, ro5_counts.values, color='skyblue', edgecolor='black')
    axes[4].set_xlabel('Violations', fontweight='bold')
    axes[4].set_ylabel('Count', fontweight='bold')
    axes[4].set_title('Lipinski RO5 (Full)', fontweight='bold')
    axes[4].grid(alpha=0.3, axis='y')
    
    # 6. QED vs MW
    scatter = axes[5].scatter(df_plot['MW'], df_plot['QED'], 
                             c=df_plot['Score'], cmap='viridis', 
                             alpha=0.3, s=5)
    axes[5].axvline(500, color='red', linestyle='--', alpha=0.5)
    axes[5].axhline(0.5, color='orange', linestyle='--', alpha=0.5)
    axes[5].set_xlabel('MW', fontweight='bold')
    axes[5].set_ylabel('QED', fontweight='bold')
    axes[5].set_title('QED vs MW', fontweight='bold')
    plt.colorbar(scatter, ax=axes[5], label='Score')
    axes[5].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/full_evaluation_plots.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ 图表已保存: {OUTPUT_DIR}/full_evaluation_plots.png")


def save_high_quality(df_results):
    """保存高质量分子"""
    
    high_quality = df_results[
        (df_results['QED'] > 0.6) &
        (df_results['RO5_Violations'] <= 1) &
        (df_results['MW'].between(200, 600))
    ].copy()
    
    high_quality = high_quality.sort_values('QED', ascending=False)
    
    hq_file = f"{OUTPUT_DIR}/high_quality_molecules_full.csv"
    high_quality.to_csv(hq_file, index=False)
    
    print(f"✓ 高质量分子: {hq_file} ({len(high_quality)} 个)")
    
    # 保存Top 1000
    top1000 = high_quality.head(1000)
    top1000.to_csv(f"{OUTPUT_DIR}/top1000_molecules.csv", index=False)
    print(f"✓ Top 1000: {OUTPUT_DIR}/top1000_molecules.csv")


# ============================================================================
# 运行
# ============================================================================

if __name__ == "__main__":
    main()