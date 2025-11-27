#!/usr/bin/env python3
"""
比较所有REINVENT4运行结果并生成PDF报告
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from datetime import datetime
from pathlib import Path

# 设置中文字体（如果需要）
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def find_all_runs(base_dir="experiments/runs"):
    """查找所有run目录"""
    runs = {}
    for run_dir in glob.glob(os.path.join(base_dir, "*")):
        if os.path.isdir(run_dir):
            run_name = os.path.basename(run_dir)
            runs[run_name] = run_dir
    return runs

def load_run_results(run_dir):
    """加载单个run的结果"""
    results = {
        'run_dir': run_dir,
        'config': None,
        'results_csv': None,
        'candidates_gold': None,
        'candidates_high': None,
        'candidates_good': None,
        'promising': None,
        'training_log': None
    }
    
    # 查找CSV文件
    csv_files = {
        'results_1.csv': 'results_csv',
        'candidates_gold.csv': 'candidates_gold',
        'candidates_high.csv': 'candidates_high',
        'candidates_good.csv': 'candidates_good',
        'promising_candidates.csv': 'promising'
    }
    
    for csv_name, key in csv_files.items():
        csv_path = os.path.join(run_dir, csv_name)
        if os.path.exists(csv_path):
            try:
                df = pd.read_csv(csv_path)
                results[key] = df
            except Exception as e:
                print(f"  ⚠️  无法读取 {csv_name}: {e}")
    
    # 查找配置文件
    config_path = os.path.join(run_dir, 'config.toml')
    if os.path.exists(config_path):
        results['config'] = config_path
    
    # 查找训练日志
    log_path = os.path.join(run_dir, 'training.log')
    if os.path.exists(log_path):
        results['training_log'] = log_path
    
    return results

def get_run_statistics(results):
    """获取run的统计信息"""
    stats = {}
    
    if results['results_csv'] is not None:
        df = results['results_csv']
        stats['total_molecules'] = len(df)
        if 'total_score' in df.columns:
            stats['mean_score'] = df['total_score'].mean()
            stats['max_score'] = df['total_score'].max()
            stats['min_score'] = df['total_score'].min()
    
    if results['candidates_gold'] is not None:
        stats['gold_count'] = len(results['candidates_gold'])
    else:
        stats['gold_count'] = 0
    
    if results['candidates_high'] is not None:
        stats['high_count'] = len(results['candidates_high'])
    else:
        stats['high_count'] = 0
    
    if results['candidates_good'] is not None:
        stats['good_count'] = len(results['candidates_good'])
    else:
        stats['good_count'] = 0
    
    if results['promising'] is not None:
        stats['promising_count'] = len(results['promising'])
    else:
        stats['promising_count'] = 0
    
    return stats

def create_comparison_report(all_results, output_pdf="runs_comparison_report.pdf"):
    """创建比较报告PDF"""
    
    with PdfPages(output_pdf) as pdf:
        # 第1页：概览
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('REINVENT4 Runs Comparison Report', fontsize=16, fontweight='bold')
        
        ax = fig.add_subplot(111)
        ax.axis('off')
        
        report_text = f"""
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
Total Runs: {len(all_results)}

═══════════════════════════════════════════════════════════════
SUMMARY OF ALL RUNS
═══════════════════════════════════════════════════════════════
"""
        
        # 收集统计数据
        comparison_data = []
        for run_name, results in sorted(all_results.items()):
            stats = get_run_statistics(results)
            comparison_data.append({
                'Run': run_name,
                **stats
            })
            
            report_text += f"\n{run_name}:\n"
            report_text += f"  Gold candidates: {stats.get('gold_count', 0)}\n"
            report_text += f"  High candidates: {stats.get('high_count', 0)}\n"
            report_text += f"  Good candidates: {stats.get('good_count', 0)}\n"
            report_text += f"  Promising: {stats.get('promising_count', 0)}\n"
            if 'total_molecules' in stats:
                report_text += f"  Total molecules: {stats['total_molecules']}\n"
            if 'mean_score' in stats:
                report_text += f"  Mean score: {stats['mean_score']:.4f}\n"
                report_text += f"  Max score: {stats['max_score']:.4f}\n"
        
        ax.text(0.05, 0.95, report_text, transform=ax.transAxes,
                fontfamily='monospace', fontsize=9, verticalalignment='top')
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
        
        # 第2页：候选分子数量对比图
        if comparison_data:
            comp_df = pd.DataFrame(comparison_data)
            
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
            fig.suptitle('Candidates Comparison Across Runs', fontsize=14, fontweight='bold')
            
            # 图1：Gold candidates
            if 'gold_count' in comp_df.columns:
                ax = axes[0, 0]
                comp_df.plot(x='Run', y='gold_count', kind='bar', ax=ax, 
                            color='gold', legend=False)
                ax.set_title('Gold Candidates')
                ax.set_ylabel('Count')
                ax.tick_params(axis='x', rotation=45)
            
            # 图2：High candidates
            if 'high_count' in comp_df.columns:
                ax = axes[0, 1]
                comp_df.plot(x='Run', y='high_count', kind='bar', ax=ax,
                            color='orange', legend=False)
                ax.set_title('High Quality Candidates')
                ax.set_ylabel('Count')
                ax.tick_params(axis='x', rotation=45)
            
            # 图3：Good candidates
            if 'good_count' in comp_df.columns:
                ax = axes[1, 0]
                comp_df.plot(x='Run', y='good_count', kind='bar', ax=ax,
                            color='lightblue', legend=False)
                ax.set_title('Good Candidates')
                ax.set_ylabel('Count')
                ax.tick_params(axis='x', rotation=45)
            
            # 图4：总体对比
            ax = axes[1, 1]
            metrics = ['gold_count', 'high_count', 'good_count']
            existing_metrics = [m for m in metrics if m in comp_df.columns]
            if existing_metrics:
                comp_df[['Run'] + existing_metrics].set_index('Run').plot(
                    kind='bar', ax=ax, stacked=False)
                ax.set_title('Overall Comparison')
                ax.set_ylabel('Count')
                ax.legend(['Gold', 'High', 'Good'])
                ax.tick_params(axis='x', rotation=45)
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        # 第3页+：每个run的详细信息
        for run_name, results in sorted(all_results.items()):
            fig = plt.figure(figsize=(11, 8.5))
            fig.suptitle(f'Detailed Results: {run_name}', fontsize=14, fontweight='bold')
            
            ax = fig.add_subplot(111)
            ax.axis('off')
            
            detail_text = f"""
Run: {run_name}
Directory: {results['run_dir']}

═══════════════════════════════════════════════════════════════
STATISTICS
═══════════════════════════════════════════════════════════════
"""
            stats = get_run_statistics(results)
            for key, value in stats.items():
                detail_text += f"{key}: {value}\n"
            
            # 显示top候选分子
            if results['candidates_gold'] is not None:
                detail_text += f"\n{'='*63}\nTOP GOLD CANDIDATES\n{'='*63}\n"
                df = results['candidates_gold'].head(10)
                for i, row in df.iterrows():
                    detail_text += f"\n#{i+1}\n"
                    detail_text += f"SMILES: {row.get('SMILES', 'N/A')[:60]}\n"
                    if 'Score' in row:
                        detail_text += f"Score: {row['Score']:.4f}\n"
                    if 'QSAR_Score' in row:
                        detail_text += f"QSAR: {row['QSAR_Score']:.4f}\n"
            
            ax.text(0.05, 0.95, detail_text, transform=ax.transAxes,
                    fontfamily='monospace', fontsize=8, verticalalignment='top')
            
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
        
        # 元数据页
        d = pdf.infodict()
        d['Title'] = 'REINVENT4 Runs Comparison Report'
        d['Author'] = 'REINVENT4 Analysis Script'
        d['Subject'] = 'Comparison of multiple REINVENT4 runs'
        d['Keywords'] = 'REINVENT4, Drug Discovery, DENV'
        d['CreationDate'] = datetime.now()
    
    print(f"\n✅ 报告已生成: {output_pdf}")

def main():
    print("="*70)
    print("REINVENT4 Runs Comparison Tool")
    print("="*70)
    
    # 查找所有runs
    print("\n📁 搜索所有运行目录...")
    runs = find_all_runs()
    print(f"   找到 {len(runs)} 个运行目录")
    
    # 加载所有结果
    print("\n📊 加载运行结果...")
    all_results = {}
    for run_name, run_dir in sorted(runs.items()):
        print(f"\n  处理: {run_name}")
        results = load_run_results(run_dir)
        all_results[run_name] = results
        
        # 显示简要信息
        stats = get_run_statistics(results)
        print(f"    Gold: {stats.get('gold_count', 0)}, "
              f"High: {stats.get('high_count', 0)}, "
              f"Good: {stats.get('good_count', 0)}")
    
    # 生成报告
    print("\n📄 生成PDF报告...")
    create_comparison_report(all_results)
    
    # 生成CSV摘要
    print("\n📋 生成CSV摘要...")
    summary_data = []
    for run_name, results in sorted(all_results.items()):
        stats = get_run_statistics(results)
        summary_data.append({
            'Run_Name': run_name,
            'Run_Directory': results['run_dir'],
            **stats
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = "runs_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"✅ 摘要已保存: {summary_csv}")
    
    print("\n" + "="*70)
    print("✅ 完成！")
    print("="*70)
    print(f"\n生成的文件:")
    print(f"  1. runs_comparison_report.pdf  - 详细PDF报告")
    print(f"  2. runs_summary.csv            - CSV摘要表格")
    print()

if __name__ == "__main__":
    main()
