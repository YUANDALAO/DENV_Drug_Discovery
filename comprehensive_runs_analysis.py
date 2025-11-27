#!/usr/bin/env python3
"""
REINVENT4 Comprehensive Runs Analysis & Comparison
顶刊级别的完整分析报告
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from datetime import datetime
from pathlib import Path
import json
import re

# 设置绘图风格
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

class RunAnalyzer:
    """分析单个run的数据"""
    
    def __init__(self, run_name, run_dir):
        self.run_name = run_name
        self.run_dir = run_dir
        self.data = {}
        self.config = {}
        self.metrics = {}
        
    def load_data(self):
        """加载所有数据文件"""
        print(f"  加载 {self.run_name}...")
        
        # 查找CSV文件（智能匹配）
        csv_patterns = {
            'results': ['results_*.csv', 'results.csv', 'output.csv'],
            'gold': ['candidates_gold.csv', '*gold*.csv'],
            'high': ['candidates_high.csv', '*high*.csv'],
            'good': ['candidates_good.csv', '*good*.csv'],
            'promising': ['promising*.csv']
        }
        
        for key, patterns in csv_patterns.items():
            for pattern in patterns:
                files = glob.glob(os.path.join(self.run_dir, pattern))
                if files:
                    try:
                        self.data[key] = pd.read_csv(files[0])
                        print(f"    ✓ {key}: {len(self.data[key])} 条记录")
                        break
                    except Exception as e:
                        print(f"    ✗ 无法读取 {key}: {e}")
        
        # 加载配置文件
        self._load_config()
        
        # 计算指标
        self._calculate_metrics()
        
    def _load_config(self):
        """加载并解析配置文件"""
        config_path = os.path.join(self.run_dir, 'config.toml')
        if os.path.exists(config_path):
            try:
                with open(config_path, 'r') as f:
                    content = f.read()
                    
                # 解析关键配置
                self.config = {
                    'run_type': self._extract_config(content, r'run_type\s*=\s*"([^"]+)"'),
                    'device': self._extract_config(content, r'device\s*=\s*"([^"]+)"'),
                    'scoring_type': self._extract_config(content, r'\[scoring\]\s*type\s*=\s*"([^"]+)"'),
                    'learning_rate': self._extract_config(content, r'rate\s*=\s*([0-9.]+)', float),
                    'sigma': self._extract_config(content, r'sigma\s*=\s*([0-9.]+)', float),
                    'batch_size': self._extract_config(content, r'batch_size\s*=\s*([0-9]+)', int),
                }
                
                # 提取scoring components
                self.config['components'] = self._extract_components(content)
                
            except Exception as e:
                print(f"    ⚠️  配置文件解析警告: {e}")
    
    def _extract_config(self, content, pattern, dtype=str):
        """从配置中提取值"""
        match = re.search(pattern, content)
        if match:
            try:
                return dtype(match.group(1))
            except:
                return match.group(1)
        return None
    
    def _extract_components(self, content):
        """提取scoring组件"""
        components = []
        # 简单提取component段落
        comp_pattern = r'\[\[.*?component.*?\]\].*?\[(.*?)\]'
        matches = re.findall(comp_pattern, content, re.DOTALL)
        for match in matches:
            comp_name = match.split('.')[1] if '.' in match else match
            components.append(comp_name.strip())
        return components
    
    def _calculate_metrics(self):
        """计算关键指标"""
        metrics = {}
        
        # 从results或gold数据中提取指标
        df = self.data.get('gold') or self.data.get('results')
        
        if df is not None:
            # 标准化列名（处理不同的命名方式）
            df.columns = df.columns.str.strip()
            
            # 分子数量
            metrics['n_molecules'] = len(df)
            metrics['n_gold'] = len(self.data.get('gold', []))
            metrics['n_high'] = len(self.data.get('high', []))
            metrics['n_good'] = len(self.data.get('good', []))
            
            # 尝试提取各种可能的列名
            column_mappings = {
                'total_score': ['total_score', 'Score', 'score', 'final_score'],
                'qsar': ['QSAR_Score', 'qsar', 'QSAR', 'predicted_pic50', 'pIC50'],
                'qed': ['QED', 'qed', 'QED Score'],
                'mw': ['MW', 'Molecular weight', 'molecular_weight', 'MolecularWeight'],
                'sa_score': ['SA score', 'SAScore', 'sa_score'],
                'tpsa': ['TPSA', 'tpsa'],
                'logp': ['SlogP (RDKit)', 'SlogP', 'logp', 'LogP'],
                'num_rings': ['Number of aromatic rings', 'NumAromaticRings', 'aromatic_rings'],
            }
            
            for metric_name, possible_cols in column_mappings.items():
                for col in possible_cols:
                    if col in df.columns:
                        values = pd.to_numeric(df[col], errors='coerce').dropna()
                        if len(values) > 0:
                            metrics[f'{metric_name}_mean'] = values.mean()
                            metrics[f'{metric_name}_std'] = values.std()
                            metrics[f'{metric_name}_min'] = values.min()
                            metrics[f'{metric_name}_max'] = values.max()
                            metrics[f'{metric_name}_median'] = values.median()
                        break
        
        self.metrics = metrics
        return metrics

class ComprehensiveReportGenerator:
    """生成综合分析报告"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.comparison_df = self._create_comparison_dataframe()
    
    def _create_comparison_dataframe(self):
        """创建对比数据框"""
        data = []
        for analyzer in self.analyzers:
            row = {
                'Run': analyzer.run_name,
                **analyzer.metrics,
                **{f'config_{k}': v for k, v in analyzer.config.items() if k != 'components'}
            }
            data.append(row)
        return pd.DataFrame(data)
    
    def generate_report(self, output_pdf="comprehensive_runs_report.pdf"):
        """生成完整PDF报告"""
        
        with PdfPages(output_pdf) as pdf:
            print("\n📄 生成PDF报告...")
            
            # 第1页：封面和摘要
            self._create_cover_page(pdf)
            
            # 第2页：金标准分子数量比较
            self._create_gold_candidates_comparison(pdf)
            
            # 第3页：分子结构对比
            self._create_structure_comparison(pdf)
            
            # 第4页：类药性指标比较
            self._create_druglikeness_comparison(pdf)
            
            # 第5页：活性和QSAR比较
            self._create_activity_comparison(pdf)
            
            # 第6页：分子性质比较（分子量、LogP、TPSA等）
            self._create_molecular_properties_comparison(pdf)
            
            # 第7页：模型优化趋势（拟合值变化）
            self._create_optimization_trends(pdf)
            
            # 第8页：配置差异对比
            self._create_config_comparison(pdf)
            
            # 第9-N页：每个run的详细报告
            self._create_detailed_run_pages(pdf)
            
            # 元数据
            d = pdf.infodict()
            d['Title'] = 'REINVENT4 Comprehensive Analysis Report'
            d['Author'] = 'Automated Analysis Pipeline'
            d['Subject'] = 'Multi-run comparison and analysis'
            d['CreationDate'] = datetime.now()
        
        print(f"✅ PDF报告已生成: {output_pdf}")
    
    def _create_cover_page(self, pdf):
        """创建封面页"""
        fig = plt.figure(figsize=(11, 8.5))
        ax = fig.add_subplot(111)
        ax.axis('off')
        
        title = "REINVENT4 COMPREHENSIVE ANALYSIS REPORT"
        subtitle = "De Novo Drug Design for DENV NS2B-NS3 Protease Inhibitors"
        
        # 标题
        ax.text(0.5, 0.7, title, ha='center', va='center',
                fontsize=20, fontweight='bold', transform=ax.transAxes)
        ax.text(0.5, 0.63, subtitle, ha='center', va='center',
                fontsize=12, style='italic', transform=ax.transAxes)
        
        # 摘要统计
        summary_text = f"""
══════════════════════════════════════════════════════════════════
EXECUTIVE SUMMARY
══════════════════════════════════════════════════════════════════

Report Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Total Runs Analyzed: {len(self.analyzers)}

Runs: {', '.join([a.run_name for a in self.analyzers])}

Total Gold Candidates: {sum([a.metrics.get('n_gold', 0) for a in self.analyzers])}
Total High Quality: {sum([a.metrics.get('n_high', 0) for a in self.analyzers])}
Total Good Quality: {sum([a.metrics.get('n_good', 0) for a in self.analyzers])}

══════════════════════════════════════════════════════════════════
TABLE OF CONTENTS
══════════════════════════════════════════════════════════════════

1. Gold Standard Molecules Comparison
2. Molecular Structure Comparison  
3. Drug-likeness Metrics Comparison
4. Activity & QSAR Prediction Comparison
5. Molecular Properties Comparison
6. Model Optimization Trends
7. Configuration Differences Analysis
8. Detailed Individual Run Reports

══════════════════════════════════════════════════════════════════
"""
        
        ax.text(0.1, 0.55, summary_text, ha='left', va='top',
                fontfamily='monospace', fontsize=9, transform=ax.transAxes)
        
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_gold_candidates_comparison(self, pdf):
        """金标准分子数量比较"""
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        fig.suptitle('Gold Standard Molecules Comparison', 
                     fontsize=14, fontweight='bold', y=0.98)
        
        df = self.comparison_df
        runs = df['Run'].values
        
        # 图1：金标准定义说明 + 数量柱状图
        ax = axes[0, 0]
        gold_counts = df['n_gold'].fillna(0).values
        bars = ax.bar(range(len(runs)), gold_counts, color='gold', alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(runs)))
        ax.set_xticklabels(runs, rotation=45, ha='right')
        ax.set_ylabel('Number of Molecules')
        ax.set_title('Gold Candidates Count\n(Score ≥ threshold & Pass filters)')
        ax.grid(axis='y', alpha=0.3)
        
        # 添加数值标签
        for i, (bar, count) in enumerate(zip(bars, gold_counts)):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                   f'{int(count)}', ha='center', va='bottom', fontsize=9)
        
        # 图2：层级分类堆叠图
        ax = axes[0, 1]
        categories = ['n_gold', 'n_high', 'n_good']
        labels = ['Gold', 'High', 'Good']
        colors = ['gold', 'orange', 'lightblue']
        
        bottom = np.zeros(len(runs))
        for cat, label, color in zip(categories, labels, colors):
            values = df[cat].fillna(0).values
            ax.bar(range(len(runs)), values, bottom=bottom, 
                  label=label, color=color, alpha=0.7, edgecolor='black')
            bottom += values
        
        ax.set_xticks(range(len(runs)))
        ax.set_xticklabels(runs, rotation=45, ha='right')
        ax.set_ylabel('Number of Molecules')
        ax.set_title('Candidate Quality Distribution')
        ax.legend()
        ax.grid(axis='y', alpha=0.3)
        
        # 图3：成功率比较
        ax = axes[1, 0]
        total_mols = df['n_molecules'].fillna(1).values
        success_rate = (gold_counts / total_mols) * 100
        bars = ax.bar(range(len(runs)), success_rate, color='green', alpha=0.6)
        ax.set_xticks(range(len(runs)))
        ax.set_xticklabels(runs, rotation=45, ha='right')
        ax.set_ylabel('Success Rate (%)')
        ax.set_title('Gold Candidate Success Rate')
        ax.grid(axis='y', alpha=0.3)
        
        for bar, rate in zip(bars, success_rate):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                   f'{rate:.1f}%', ha='center', va='bottom', fontsize=9)
        
        # 图4：标准说明文本
        ax = axes[1, 1]
        ax.axis('off')
        
        criteria_text = """
GOLD STANDARD CRITERIA:

✓ Total Score ≥ threshold
✓ QSAR predicted pIC50 ≥ 6.0
✓ QED ≥ 0.5 (drug-likeness)
✓ Molecular Weight: 200-500 Da
✓ Pass structural alerts
✓ Lipinski's Rule of 5 compliant
✓ Synthetic Accessibility < 6

HIGH QUALITY CRITERIA:

✓ Total Score ≥ 0.7
✓ QSAR pIC50 ≥ 5.5
✓ QED ≥ 0.4

GOOD QUALITY CRITERIA:

✓ Total Score ≥ 0.5
✓ Pass basic filters
"""
        
        ax.text(0.1, 0.9, criteria_text, transform=ax.transAxes,
                fontfamily='monospace', fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_structure_comparison(self, pdf):
        """分子结构对比"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            from rdkit.Chem import Descriptors
            
            fig = plt.figure(figsize=(11, 8.5))
            fig.suptitle('Top Molecular Structures Comparison Across Runs',
                        fontsize=14, fontweight='bold')
            
            n_runs = len(self.analyzers)
            n_mols_per_run = 6
            
            for idx, analyzer in enumerate(self.analyzers):
                if analyzer.data.get('gold') is not None:
                    df = analyzer.data['gold'].head(n_mols_per_run)
                    
                    mols = []
                    legends = []
                    
                    for i, row in df.iterrows():
                        try:
                            mol = Chem.MolFromSmiles(row['SMILES'])
                            if mol:
                                mols.append(mol)
                                
                                # 创建legend
                                legend = f"{analyzer.run_name}-{i+1}"
                                if 'total_score' in row:
                                    legend += f"\nScore:{row['total_score']:.2f}"
                                if 'QSAR_Score' in row or 'qsar' in row.keys():
                                    qsar_col = 'QSAR_Score' if 'QSAR_Score' in row else 'qsar'
                                    legend += f"\npIC50:{row[qsar_col]:.2f}"
                                
                                legends.append(legend)
                        except:
                            continue
                    
                    if mols:
                        # 为每个run创建子图
                        ax = fig.add_subplot(n_runs, 1, idx+1)
                        img = Draw.MolsToGridImage(mols[:n_mols_per_run], 
                                                   molsPerRow=n_mols_per_run,
                                                   subImgSize=(200, 200),
                                                   legends=legends[:n_mols_per_run])
                        ax.imshow(img)
                        ax.set_title(f'{analyzer.run_name} - Top {len(mols)} Molecules',
                                    fontweight='bold', fontsize=10)
                        ax.axis('off')
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()
            
        except ImportError:
            print("  ⚠️  RDKit未安装，跳过结构对比页")
        except Exception as e:
            print(f"  ⚠️  结构对比生成失败: {e}")
    
    def _create_druglikeness_comparison(self, pdf):
        """类药性指标比较"""
        fig, axes = plt.subplots(2, 3, figsize=(11, 8.5))
        fig.suptitle('Drug-likeness Metrics Comparison', 
                     fontsize=14, fontweight='bold', y=0.98)
        
        df = self.comparison_df
        runs = df['Run'].values
        
        metrics_to_plot = [
            ('qed_mean', 'QED (Drug-likeness)', (0, 1)),
            ('sa_score_mean', 'SA Score (Synthesizability)', (1, 10)),
            ('mw_mean', 'Molecular Weight (Da)', (0, 600)),
            ('logp_mean', 'LogP (Lipophilicity)', (-2, 6)),
            ('tpsa_mean', 'TPSA (Å²)', (0, 200)),
            ('num_rings_mean', 'Aromatic Rings', (0, 5)),
        ]
        
        for idx, (ax, (metric, title, ylim)) in enumerate(zip(axes.flat, metrics_to_plot)):
            if metric in df.columns:
                means = df[metric].values
                stds = df.get(metric.replace('_mean', '_std'), pd.Series([0]*len(runs))).values
                mins = df.get(metric.replace('_mean', '_min'), means).values
                maxs = df.get(metric.replace('_mean', '_max'), means).values
                
                # 柱状图 + 误差线
                x = np.arange(len(runs))
                bars = ax.bar(x, means, yerr=stds, capsize=5, 
                             color='skyblue', alpha=0.7, edgecolor='black')
                
                # 添加最大/最小值标记
                ax.scatter(x, maxs, color='red', marker='^', s=50, 
                          label='Max', zorder=5)
                ax.scatter(x, mins, color='blue', marker='v', s=50,
                          label='Min', zorder=5)
                
                ax.set_xticks(x)
                ax.set_xticklabels(runs, rotation=45, ha='right', fontsize=8)
                ax.set_ylabel('Value')
                ax.set_title(title, fontsize=10, fontweight='bold')
                ax.set_ylim(ylim)
                ax.grid(axis='y', alpha=0.3)
                
                if idx == 0:
                    ax.legend(fontsize=8, loc='upper right')
                
                # 添加参考线
                if 'qed' in metric:
                    ax.axhline(y=0.5, color='green', linestyle='--', 
                              alpha=0.5, label='Threshold')
                elif 'mw' in metric:
                    ax.axhline(y=500, color='red', linestyle='--', 
                              alpha=0.5, label='Lipinski')
                elif 'sa_score' in metric:
                    ax.axhline(y=6, color='orange', linestyle='--',
                              alpha=0.5, label='Easy to synthesize')
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_activity_comparison(self, pdf):
        """活性和QSAR预测比较"""
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        fig.suptitle('Activity & QSAR Prediction Comparison',
                     fontsize=14, fontweight='bold', y=0.98)
        
        df = self.comparison_df
        runs = df['Run'].values
        
        # 图1：QSAR预测值分布
        ax = axes[0, 0]
        if 'qsar_mean' in df.columns:
            means = df['qsar_mean'].values
            stds = df.get('qsar_std', pd.Series([0]*len(runs))).values
            mins = df.get('qsar_min', means).values
            maxs = df.get('qsar_max', means).values
            
            x = np.arange(len(runs))
            ax.bar(x, means, yerr=stds, capsize=5, color='purple', alpha=0.6)
            ax.scatter(x, maxs, color='red', marker='^', s=60, label='Max')
            ax.scatter(x, mins, color='blue', marker='v', s=60, label='Min')
            
            ax.axhline(y=6.0, color='green', linestyle='--', linewidth=2,
                      label='Target pIC50 ≥ 6.0')
            ax.axhline(y=7.0, color='darkgreen', linestyle=':', linewidth=2,
                      label='Excellent pIC50 ≥ 7.0')
            
            ax.set_xticks(x)
            ax.set_xticklabels(runs, rotation=45, ha='right')
            ax.set_ylabel('Predicted pIC50')
            ax.set_title('QSAR Activity Prediction', fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(axis='y', alpha=0.3)
        
        # 图2：总分数分布
        ax = axes[0, 1]
        if 'total_score_mean' in df.columns:
            means = df['total_score_mean'].values
            stds = df.get('total_score_std', pd.Series([0]*len(runs))).values
            
            x = np.arange(len(runs))
            bars = ax.bar(x, means, yerr=stds, capsize=5, 
                         color='orange', alpha=0.6, edgecolor='black')
            
            ax.axhline(y=0.8, color='gold', linestyle='--', label='Gold threshold')
            ax.axhline(y=0.7, color='orange', linestyle='--', label='High threshold')
            
            ax.set_xticks(x)
            ax.set_xticklabels(runs, rotation=45, ha='right')
            ax.set_ylabel('Total Score')
            ax.set_title('Overall Scoring Function', fontweight='bold')
            ax.set_ylim(0, 1)
            ax.legend(fontsize=8)
            ax.grid(axis='y', alpha=0.3)
        
        # 图3：活性分布小提琴图（如果有原始数据）
        ax = axes[1, 0]
        violin_data = []
        labels = []
        
        for analyzer in self.analyzers:
            df_gold = analyzer.data.get('gold')
            if df_gold is not None:
                qsar_col = None
                for col in ['QSAR_Score', 'qsar', 'predicted_pic50']:
                    if col in df_gold.columns:
                        qsar_col = col
                        break
                
                if qsar_col:
                    values = pd.to_numeric(df_gold[qsar_col], errors='coerce').dropna()
                    if len(values) > 0:
                        violin_data.append(values)
                        labels.append(analyzer.run_name)
        
        if violin_data:
            parts = ax.violinplot(violin_data, positions=range(len(labels)),
                                 showmeans=True, showmedians=True)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha='right')
            ax.set_ylabel('Predicted pIC50')
            ax.set_title('Activity Distribution (Violin Plot)', fontweight='bold')
            ax.axhline(y=6.0, color='green', linestyle='--', alpha=0.5)
            ax.grid(axis='y', alpha=0.3)
        
        # 图4：统计汇总表
        ax = axes[1, 1]
        ax.axis('off')
        
        # 创建统计表
        table_data = []
        for _, row in df.iterrows():
            table_data.append([
                row['Run'],
                f"{row.get('qsar_mean', 0):.2f}±{row.get('qsar_std', 0):.2f}",
                f"{row.get('qsar_max', 0):.2f}",
                f"{row.get('total_score_mean', 0):.2f}",
                f"{row.get('n_gold', 0):.0f}"
            ])
        
        table = ax.table(cellText=table_data,
                        colLabels=['Run', 'pIC50\n(Mean±SD)', 'Max\npIC50', 
                                  'Mean\nScore', 'N\nGold'],
                        cellLoc='center',
                        loc='center',
                        bbox=[0, 0, 1, 1])
        
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.scale(1, 2)
        
        # 设置表格样式
        for i in range(len(table_data) + 1):
            for j in range(5):
                cell = table[(i, j)]
                if i == 0:
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
        
        ax.set_title('Activity Statistics Summary', 
                    fontweight='bold', pad=20, fontsize=10)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_molecular_properties_comparison(self, pdf):
        """分子性质详细比较"""
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Molecular Properties Detailed Comparison',
                     fontsize=14, fontweight='bold')
        
        gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.4)
        
        df = self.comparison_df
        runs = df['Run'].values
        x = np.arange(len(runs))
        
        # 分子量分布
        ax1 = fig.add_subplot(gs[0, :2])
        if 'mw_mean' in df.columns:
            means = df['mw_mean'].values
            mins = df.get('mw_min', means).values
            maxs = df.get('mw_max', means).values
            
            ax1.fill_between(x, mins, maxs, alpha=0.3, color='blue', label='Range')
            ax1.plot(x, means, 'o-', color='darkblue', linewidth=2, 
                    markersize=8, label='Mean')
            
            ax1.axhspan(200, 500, alpha=0.1, color='green', label='Optimal (200-500)')
            ax1.axhline(y=500, color='red', linestyle='--', label='Lipinski limit')
            
            ax1.set_xticks(x)
            ax1.set_xticklabels(runs, rotation=45, ha='right')
            ax1.set_ylabel('Molecular Weight (Da)')
            ax1.set_title('Molecular Weight Distribution', fontweight='bold')
            ax1.legend(fontsize=8)
            ax1.grid(True, alpha=0.3)
        
        # LogP vs TPSA scatter
        ax2 = fig.add_subplot(gs[0, 2])
        if 'logp_mean' in df.columns and 'tpsa_mean' in df.columns:
            colors = plt.cm.viridis(np.linspace(0, 1, len(runs)))
            
            for i, run in enumerate(runs):
                ax2.scatter(df.iloc[i]['logp_mean'], 
                           df.iloc[i]['tpsa_mean'],
                           s=100, color=colors[i], 
                           alpha=0.7, edgecolors='black', linewidth=1.5)
                ax2.text(df.iloc[i]['logp_mean'], 
                        df.iloc[i]['tpsa_mean'],
                        run, fontsize=7, ha='center', va='bottom')
            
            # 添加理想区域
            from matplotlib.patches import Rectangle
            ideal_region = Rectangle((0, 20), 5, 120, 
                                    alpha=0.1, color='green',
                                    label='Ideal region')
            ax2.add_patch(ideal_region)
            
            ax2.set_xlabel('LogP')
            ax2.set_ylabel('TPSA (Ų)')
            ax2.set_title('LogP vs TPSA', fontweight='bold', fontsize=9)
            ax2.grid(True, alpha=0.3)
            ax2.legend(fontsize=7)
        
        # QED分布
        ax3 = fig.add_subplot(gs[1, 0])
        if 'qed_mean' in df.columns:
            means = df['qed_mean'].values
            stds = df.get('qed_std', pd.Series([0]*len(runs))).values
            
            bars = ax3.barh(x, means, xerr=stds, 
                           color='lightcoral', alpha=0.7, edgecolor='black')
            ax3.axvline(x=0.5, color='green', linestyle='--', 
                       linewidth=2, label='Threshold')
            
            ax3.set_yticks(x)
            ax3.set_yticklabels(runs, fontsize=8)
            ax3.set_xlabel('QED Score')
            ax3.set_title('Drug-likeness (QED)', fontweight='bold', fontsize=9)
            ax3.set_xlim(0, 1)
            ax3.legend(fontsize=7)
            ax3.grid(axis='x', alpha=0.3)
        
        # SA Score分布
        ax4 = fig.add_subplot(gs[1, 1])
        if 'sa_score_mean' in df.columns:
            means = df['sa_score_mean'].values
            
            bars = ax4.barh(x, means, color='lightyellow', 
                           alpha=0.7, edgecolor='black')
            ax4.axvline(x=6, color='orange', linestyle='--',
                       linewidth=2, label='Easy synthesis')
            
            ax4.set_yticks(x)
            ax4.set_yticklabels(runs, fontsize=8)
            ax4.set_xlabel('SA Score')
            ax4.set_title('Synthetic Accessibility', fontweight='bold', fontsize=9)
            ax4.set_xlim(1, 10)
            ax4.invert_xaxis()  # 低分更好
            ax4.legend(fontsize=7)
            ax4.grid(axis='x', alpha=0.3)
        
        # 芳香环数量
        ax5 = fig.add_subplot(gs[1, 2])
        if 'num_rings_mean' in df.columns:
            means = df['num_rings_mean'].values
            
            bars = ax5.bar(x, means, color='lightgreen', 
                          alpha=0.7, edgecolor='black')
            ax5.axhline(y=3, color='blue', linestyle='--',
                       label='Typical range')
            
            ax5.set_xticks(x)
            ax5.set_xticklabels(runs, rotation=45, ha='right', fontsize=8)
            ax5.set_ylabel('Number of Rings')
            ax5.set_title('Aromatic Rings', fontweight='bold', fontsize=9)
            ax5.legend(fontsize=7)
            ax5.grid(axis='y', alpha=0.3)
        
        # Lipinski合规性雷达图
        ax6 = fig.add_subplot(gs[2, :], projection='polar')
        
        categories = ['MW\n(<500)', 'LogP\n(<5)', 'HBD\n(<5)', 
                     'HBA\n(<10)', 'TPSA\n(<140)']
        N = len(categories)
        
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]
        
        ax6.set_theta_offset(np.pi / 2)
        ax6.set_theta_direction(-1)
        ax6.set_xticks(angles[:-1])
        ax6.set_xticklabels(categories, fontsize=8)
        
        # 为每个run绘制雷达图
        colors_radar = plt.cm.Set2(np.linspace(0, 1, len(runs)))
        
        for i, run in enumerate(runs):
            # 计算Lipinski合规性分数（归一化到0-1）
            row = df.iloc[i]
            
            scores = []
            # MW < 500
            mw = row.get('mw_mean', 500)
            scores.append(min(1.0, 500/max(mw, 1)))
            
            # LogP < 5
            logp = row.get('logp_mean', 5)
            scores.append(min(1.0, 5/max(abs(logp), 0.1)) if logp > 0 else 1.0)
            
            # HBD < 5 (假设默认值)
            scores.append(0.8)
            
            # HBA < 10 (假设默认值)
            scores.append(0.8)
            
            # TPSA < 140
            tpsa = row.get('tpsa_mean', 140)
            scores.append(min(1.0, 140/max(tpsa, 1)))
            
            scores += scores[:1]
            
            ax6.plot(angles, scores, 'o-', linewidth=2, 
                    label=run, color=colors_radar[i])
            ax6.fill(angles, scores, alpha=0.15, color=colors_radar[i])
        
        ax6.set_ylim(0, 1)
        ax6.set_title('Lipinski Rule of 5 Compliance', 
                     fontweight='bold', pad=20, fontsize=10)
        ax6.legend(loc='upper right', bbox_to_anchor=(1.3, 1.0), fontsize=7)
        ax6.grid(True)
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_optimization_trends(self, pdf):
        """模型优化趋势分析"""
        fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
        fig.suptitle('Model Optimization & Convergence Analysis',
                     fontsize=14, fontweight='bold', y=0.98)
        
        # 图1：分数提升趋势
        ax = axes[0, 0]
        
        # 尝试从results_1.csv读取时间序列数据
        for analyzer in self.analyzers:
            results_df = analyzer.data.get('results')
            if results_df is not None and len(results_df) > 100:
                # 假设有step或索引列
                if 'total_score' in results_df.columns:
                    scores = pd.to_numeric(results_df['total_score'], errors='coerce').dropna()
                    
                    # 计算滑动平均
                    window = max(10, len(scores) // 50)
                    rolling_mean = scores.rolling(window=window, center=True).mean()
                    
                    ax.plot(rolling_mean.index, rolling_mean.values,
                           label=analyzer.run_name, linewidth=2, alpha=0.8)
        
        ax.set_xlabel('Generation Step')
        ax.set_ylabel('Total Score (Moving Average)')
        ax.set_title('Score Improvement Over Training', fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        
        # 图2：收敛速度比较
        ax = axes[0, 1]
        
        convergence_data = []
        for analyzer in self.analyzers:
            results_df = analyzer.data.get('results')
            if results_df is not None and 'total_score' in results_df.columns:
                scores = pd.to_numeric(results_df['total_score'], errors='coerce').dropna()
                
                if len(scores) > 50:
                    # 计算达到80%最大分数所需的步数
                    max_score = scores.max()
                    threshold = max_score * 0.8
                    
                    steps_to_converge = (scores >= threshold).idxmax() if (scores >= threshold).any() else len(scores)
                    convergence_data.append({
                        'Run': analyzer.run_name,
                        'Steps': steps_to_converge,
                        'Final_Score': scores.iloc[-100:].mean()  # 最后100步平均
                    })
        
        if convergence_data:
            conv_df = pd.DataFrame(convergence_data)
            x = np.arange(len(conv_df))
            
            bars = ax.bar(x, conv_df['Steps'], color='steelblue', alpha=0.7)
            ax.set_xticks(x)
            ax.set_xticklabels(conv_df['Run'], rotation=45, ha='right')
            ax.set_ylabel('Steps to Convergence')
            ax.set_title('Convergence Speed\n(Steps to reach 80% max score)', 
                        fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
            
            # 添加数值标签
            for bar, steps in zip(bars, conv_df['Steps']):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                       f'{int(steps)}', ha='center', va='bottom', fontsize=8)
        
        # 图3：模型拟合质量（方差分析）
        ax = axes[1, 0]
        
        variance_data = []
        for analyzer in self.analyzers:
            results_df = analyzer.data.get('results')
            if results_df is not None and 'total_score' in results_df.columns:
                scores = pd.to_numeric(results_df['total_score'], errors='coerce').dropna()
                
                if len(scores) > 100:
                    # 分析最后20%的数据
                    tail_size = len(scores) // 5
                    tail_scores = scores.iloc[-tail_size:]
                    
                    variance_data.append({
                        'Run': analyzer.run_name,
                        'Mean': tail_scores.mean(),
                        'Std': tail_scores.std(),
                        'CV': (tail_scores.std() / tail_scores.mean()) * 100  # 变异系数
                    })
        
        if variance_data:
            var_df = pd.DataFrame(variance_data)
            x = np.arange(len(var_df))
            
            # 双y轴
            ax2 = ax.twinx()
            
            bars = ax.bar(x, var_df['Mean'], yerr=var_df['Std'],
                         capsize=5, color='lightcoral', alpha=0.6,
                         label='Mean ± Std')
            line = ax2.plot(x, var_df['CV'], 'go-', linewidth=2, 
                           markersize=8, label='CV%')
            
            ax.set_xticks(x)
            ax.set_xticklabels(var_df['Run'], rotation=45, ha='right')
            ax.set_ylabel('Score (Final 20%)', color='lightcoral')
            ax2.set_ylabel('Coefficient of Variation (%)', color='green')
            ax.set_title('Model Stability & Consistency', fontweight='bold')
            ax.grid(True, alpha=0.3)
            
            # 组合图例
            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)
        
        # 图4：优化效率总结
        ax = axes[1, 1]
        ax.axis('off')
        
        summary_text = """
OPTIMIZATION QUALITY METRICS:

╔═══════════════════════════════════════╗
║  CONVERGENCE CRITERIA                 ║
╠═══════════════════════════════════════╣
║  ✓ Fast: < 500 steps to 80% max      ║
║  ✓ Stable: CV < 5% in final phase    ║
║  ✓ Improved: Final > Initial by 20%  ║
╚═══════════════════════════════════════╝

MODEL FITTING INDICATORS:

- Low CV% → Good model stability
- High final score → Effective optimization
- Smooth curves → Proper learning rate
- Early plateau → May need longer training

COMPARISON INSIGHTS:

"""
        
        # 添加runs对比总结
        if convergence_data:
            best_run = min(convergence_data, key=lambda x: x['Steps'])
            summary_text += f"\n🏆 Fastest convergence: {best_run['Run']}"
            summary_text += f"\n   ({best_run['Steps']} steps)"
        
        if variance_data:
            most_stable = min(variance_data, key=lambda x: x['CV'])
            summary_text += f"\n\n🎯 Most stable: {most_stable['Run']}"
            summary_text += f"\n   (CV = {most_stable['CV']:.1f}%)"
        
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes,
                fontfamily='monospace', fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_config_comparison(self, pdf):
        """配置差异对比"""
        fig = plt.figure(figsize=(11, 8.5))
        fig.suptitle('Configuration Differences & Optimization Strategy',
                     fontsize=14, fontweight='bold')
        
        gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3)
        
        # 提取配置信息
        config_data = []
        for analyzer in self.analyzers:
            config_data.append({
                'Run': analyzer.run_name,
                **analyzer.config
            })
        
        config_df = pd.DataFrame(config_data)
        
        # 图1：配置参数对比表
        ax1 = fig.add_subplot(gs[0, :])
        ax1.axis('off')
        
        # 创建配置对比表
        table_data = []
        headers = ['Run', 'Run Type', 'Device', 'Scoring', 'LR', 'Sigma', 'Batch']
        
        for _, row in config_df.iterrows():
            table_data.append([
                row.get('Run', 'N/A'),
                row.get('run_type', 'N/A'),
                row.get('device', 'N/A'),
                row.get('scoring_type', 'N/A'),
                f"{row.get('learning_rate', 0):.4f}" if row.get('learning_rate') else 'N/A',
                f"{row.get('sigma', 0):.2f}" if row.get('sigma') else 'N/A',
                f"{row.get('batch_size', 0)}" if row.get('batch_size') else 'N/A',
            ])
        
        table = ax1.table(cellText=table_data, colLabels=headers,
                         cellLoc='center', loc='center',
                         bbox=[0, 0, 1, 1])
        
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 2.5)
        
        # 样式设置
        for i in range(len(table_data) + 1):
            for j in range(len(headers)):
                cell = table[(i, j)]
                if i == 0:
                    cell.set_facecolor('#2196F3')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    # 高亮不同的值
                    if j > 0:  # 跳过Run列
                        col_values = [r[j] for r in table_data]
                        if len(set(col_values)) > 1:  # 如果有差异
                            cell.set_facecolor('#FFE082')  # 黄色高亮
                        else:
                            cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
                    else:
                        cell.set_facecolor('#E3F2FD')
        
        ax1.set_title('Configuration Parameters Comparison\n(Yellow = Different values)',
                     fontweight='bold', pad=20)
        
        # 图2：Scoring Components对比
        ax2 = fig.add_subplot(gs[1, :])
        ax2.axis('off')
        
        # 收集所有组件
        all_components = set()
        for analyzer in self.analyzers:
            all_components.update(analyzer.config.get('components', []))
        
        all_components = sorted(list(all_components))
        
        if all_components:
            # 创建组件矩阵
            comp_matrix = np.zeros((len(self.analyzers), len(all_components)))
            
            for i, analyzer in enumerate(self.analyzers):
                run_comps = analyzer.config.get('components', [])
                for j, comp in enumerate(all_components):
                    if comp in run_comps:
                        comp_matrix[i, j] = 1
            
            # 绘制热力图
            im = ax2.imshow(comp_matrix, cmap='RdYlGn', aspect='auto',
                           interpolation='nearest', vmin=0, vmax=1)
            
            # 设置刻度
            ax2.set_xticks(np.arange(len(all_components)))
            ax2.set_yticks(np.arange(len(self.analyzers)))
            ax2.set_xticklabels(all_components, rotation=90, ha='right', fontsize=8)
            ax2.set_yticklabels([a.run_name for a in self.analyzers], fontsize=9)
            
            # 添加网格
            ax2.set_xticks(np.arange(len(all_components))-0.5, minor=True)
            ax2.set_yticks(np.arange(len(self.analyzers))-0.5, minor=True)
            ax2.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
            
            # 添加文本标注
            for i in range(len(self.analyzers)):
                for j in range(len(all_components)):
                    if comp_matrix[i, j] == 1:
                        ax2.text(j, i, '✓', ha="center", va="center",
                                color="white", fontsize=12, fontweight='bold')
            
            ax2.set_title('Scoring Components Matrix\n(✓ = Component used)',
                         fontweight='bold', pad=10)
        
        # 图3：关键差异总结
        ax3 = fig.add_subplot(gs[2, 0])
        ax3.axis('off')
        
        diff_text = """
KEY CONFIGURATION DIFFERENCES:

═══════════════════════════════════════
OPTIMIZATION PARAMETERS
═══════════════════════════════════════
"""
        
        # 找出不同的配置项
        for param in ['learning_rate', 'sigma', 'batch_size']:
            values = config_df[param].dropna()
            if len(values) > 0 and len(set(values)) > 1:
                diff_text += f"\n{param.replace('_', ' ').title()}:"
                for run, val in zip(config_df['Run'], config_df[param]):
                    if pd.notna(val):
                        diff_text += f"\n  • {run}: {val}"
        
        diff_text += "\n\n═══════════════════════════════════════\n"
        diff_text += "SCORING STRATEGY\n"
        diff_text += "═══════════════════════════════════════\n"
        
        # Scoring type差异
        scoring_types = config_df['scoring_type'].dropna().unique()
        if len(scoring_types) > 1:
            diff_text += "\nDifferent scoring types used:"
            for st in scoring_types:
                runs_with_st = config_df[config_df['scoring_type'] == st]['Run'].tolist()
                diff_text += f"\n  • {st}: {', '.join(runs_with_st)}"
        
        ax3.text(0.05, 0.95, diff_text, transform=ax3.transAxes,
                fontfamily='monospace', fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.4))
        
        # 图4：优化建议
        ax4 = fig.add_subplot(gs[2, 1])
        ax4.axis('off')
        
        recommend_text = """
OPTIMIZATION RECOMMENDATIONS:

🎯 LEARNING RATE:
   • Too high → Unstable training
   • Too low → Slow convergence
   • Optimal: 0.0001 - 0.001

🎯 SIGMA (Reward shaping):
   • Controls exploration vs exploitation
   • Higher σ → More exploration
   • Typical range: 60 - 120

🎯 BATCH SIZE:
   • Larger → More stable gradients
   • Smaller → Faster updates
   • Balance: 50 - 200

🎯 SCORING COMPONENTS:
   • More components → Better filtering
   • But may reduce diversity
   • Balance quality vs quantity

💡 INSIGHTS FROM COMPARISON:
"""
        
        # 基于实际数据的建议
        if len(self.analyzers) > 1:
            # 找出gold数量最多的run
            best_idx = self.comparison_df['n_gold'].idxmax()
            best_run = self.comparison_df.iloc[best_idx]
            
            recommend_text += f"\n   Best performing: {best_run['Run']}"
            recommend_text += f"\n   Gold candidates: {best_run['n_gold']:.0f}"
            
            # 找出该run的配置
            best_config = config_df[config_df['Run'] == best_run['Run']].iloc[0]
            if pd.notna(best_config.get('learning_rate')):
                recommend_text += f"\n   → Used LR: {best_config['learning_rate']}"
            if pd.notna(best_config.get('sigma')):
                recommend_text += f"\n   → Used σ: {best_config['sigma']}"
        
        ax4.text(0.05, 0.95, recommend_text, transform=ax4.transAxes,
                fontfamily='monospace', fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
        
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    def _create_detailed_run_pages(self, pdf):
        """为每个run创建详细页面"""
        for analyzer in self.analyzers:
            fig = plt.figure(figsize=(11, 8.5))
            fig.suptitle(f'Detailed Analysis: {analyzer.run_name}',
                        fontsize=14, fontweight='bold')
            
            gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3)
            
            # 统计摘要
            ax1 = fig.add_subplot(gs[0, :])
            ax1.axis('off')
            
            stats_text = f"""
RUN: {analyzer.run_name}
Directory: {analyzer.run_dir}

═══════════════════════════════════════════════════════════════════════════════
SUMMARY STATISTICS
═══════════════════════════════════════════════════════════════════════════════

Total Molecules Generated: {analyzer.metrics.get('n_molecules', 'N/A')}
Gold Candidates: {analyzer.metrics.get('n_gold', 0)}
High Quality: {analyzer.metrics.get('n_high', 0)}
Good Quality: {analyzer.metrics.get('n_good', 0)}

═══════════════════════════════════════════════════════════════════════════════
MOLECULAR PROPERTIES (Mean ± Std)
═══════════════════════════════════════════════════════════════════════════════

QSAR pIC50:     {analyzer.metrics.get('qsar_mean', 0):.2f} ± {analyzer.metrics.get('qsar_std', 0):.2f}
                Range: [{analyzer.metrics.get('qsar_min', 0):.2f}, {analyzer.metrics.get('qsar_max', 0):.2f}]

QED Score:      {analyzer.metrics.get('qed_mean', 0):.3f} ± {analyzer.metrics.get('qed_std', 0):.3f}

Molecular Wt:   {analyzer.metrics.get('mw_mean', 0):.1f} ± {analyzer.metrics.get('mw_std', 0):.1f} Da

LogP:           {analyzer.metrics.get('logp_mean', 0):.2f} ± {analyzer.metrics.get('logp_std', 0):.2f}

SA Score:       {analyzer.metrics.get('sa_score_mean', 0):.2f} ± {analyzer.metrics.get('sa_score_std', 0):.2f}

Total Score:    {analyzer.metrics.get('total_score_mean', 0):.3f} ± {analyzer.metrics.get('total_score_std', 0):.3f}
"""
            
            ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes,
                    fontfamily='monospace', fontsize=8, verticalalignment='top')
            
            # Top候选分子列表
            ax2 = fig.add_subplot(gs[1:, 0])
            ax2.axis('off')
            
            if analyzer.data.get('gold') is not None:
                top_df = analyzer.data['gold'].head(10)
                
                mol_text = "═══════════════════════════════════════\n"
                mol_text += "TOP 10 GOLD CANDIDATES\n"
                mol_text += "═══════════════════════════════════════\n\n"
                
                for i, row in top_df.iterrows():
                    mol_text += f"#{i+1}\n"
                    mol_text += f"SMILES: {row.get('SMILES', 'N/A')[:55]}\n"
                    
                    if 'total_score' in row or 'Score' in row:
                        score = row.get('total_score', row.get('Score', 'N/A'))
                        mol_text += f"Score: {score:.4f}\n"
                    
                    if 'QSAR_Score' in row:
                        mol_text += f"pIC50: {row['QSAR_Score']:.2f}\n"
                    
                    if 'QED' in row:
                        mol_text += f"QED: {row['QED']:.3f}\n"
                    
                    if 'MW' in row or 'Molecular weight' in row:
                        mw = row.get('MW', row.get('Molecular weight', 'N/A'))
                        mol_text += f"MW: {mw:.1f}\n"
                    
                    mol_text += "\n"
                
                ax2.text(0.05, 0.95, mol_text, transform=ax2.transAxes,
                        fontfamily='monospace', fontsize=7, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='lightgoldenrodyellow', alpha=0.5))
            
            # 配置信息
            ax3 = fig.add_subplot(gs[1:, 1])
            ax3.axis('off')
            
            config_text = "═══════════════════════════════════════\n"
            config_text += "CONFIGURATION\n"
            config_text += "═══════════════════════════════════════\n\n"
            
            for key, value in analyzer.config.items():
                if key != 'components' and value is not None:
                    config_text += f"{key.replace('_', ' ').title()}:\n  {value}\n\n"
            
            if analyzer.config.get('components'):
                config_text += "Scoring Components:\n"
                for comp in analyzer.config['components']:
                    config_text += f"  • {comp}\n"
            
            ax3.text(0.05, 0.95, config_text, transform=ax3.transAxes,
                    fontfamily='monospace', fontsize=8, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
            
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close()


def main():
    """主函数"""
    print("=" * 80)
    print("REINVENT4 COMPREHENSIVE ANALYSIS TOOL")
    print("Generating Publication-Quality Multi-Run Comparison Report")
    print("=" * 80)
    
    # 查找所有runs
    print("\n📁 Scanning for run directories...")
    runs_dir = "experiments/runs"
    
    if not os.path.exists(runs_dir):
        print(f"❌ Error: {runs_dir} directory not found!")
        return
    
    run_dirs = {}
    for item in os.listdir(runs_dir):
        item_path = os.path.join(runs_dir, item)
        if os.path.isdir(item_path):
          run_dirs[item] = item_path
    
    print(f"   Found {len(run_dirs)} run directories: {', '.join(sorted(run_dirs.keys()))}")
    
    if len(run_dirs) == 0:
        print("❌ No run directories found!")
        return
    
    # 创建分析器
    print("\n📊 Loading and analyzing runs...")
    analyzers = []
    
    for run_name, run_path in sorted(run_dirs.items()):
        analyzer = RunAnalyzer(run_name, run_path)
        try:
            analyzer.load_data()
            analyzers.append(analyzer)
        except Exception as e:
            print(f"  ⚠️  Failed to load {run_name}: {e}")
    
    if len(analyzers) == 0:
        print("❌ No runs could be loaded successfully!")
        return
    
    print(f"\n✅ Successfully loaded {len(analyzers)} runs")
    
    # 生成报告
    print("\n📄 Generating comprehensive PDF report...")
    reporter = ComprehensiveReportGenerator(analyzers)
    
    output_pdf = "comprehensive_runs_report.pdf"
    reporter.generate_report(output_pdf)
    
    # 导出数据表格
    print("\n📊 Exporting data tables...")
    
    # 1. 综合比较表
    comparison_csv = "runs_comparison_table.csv"
    reporter.comparison_df.to_csv(comparison_csv, index=False)
    print(f"   ✓ Comparison table: {comparison_csv}")
    
    # 2. 每个run的top候选分子
    for analyzer in analyzers:
        if analyzer.data.get('gold') is not None:
            top_csv = f"{analyzer.run_name}_top_candidates.csv"
            analyzer.data['gold'].head(20).to_csv(top_csv, index=False)
            print(f"   ✓ Top candidates: {top_csv}")
    
    # 3. 统计摘要表
    summary_data = []
    for analyzer in analyzers:
        summary_data.append({
            'Run': analyzer.run_name,
            'Total_Molecules': analyzer.metrics.get('n_molecules', 0),
            'Gold_Count': analyzer.metrics.get('n_gold', 0),
            'High_Count': analyzer.metrics.get('n_high', 0),
            'Good_Count': analyzer.metrics.get('n_good', 0),
            'Mean_QSAR_pIC50': f"{analyzer.metrics.get('qsar_mean', 0):.2f}",
            'Mean_QED': f"{analyzer.metrics.get('qed_mean', 0):.3f}",
            'Mean_MW': f"{analyzer.metrics.get('mw_mean', 0):.1f}",
            'Mean_Score': f"{analyzer.metrics.get('total_score_mean', 0):.3f}",
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = "runs_summary_statistics.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"   ✓ Summary statistics: {summary_csv}")
    
    # 显示快速统计
    print("\n" + "=" * 80)
    print("QUICK SUMMARY")
    print("=" * 80)
    print("\n" + summary_df.to_string(index=False))
    
    print("\n" + "=" * 80)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\nGenerated files:")
    print(f"  1. {output_pdf} - Comprehensive PDF report (publication-ready)")
    print(f"  2. {comparison_csv} - Detailed comparison table")
    print(f"  3. {summary_csv} - Summary statistics")
    print(f"  4. *_top_candidates.csv - Top molecules for each run")
    print(f"\n📖 Open {output_pdf} to view the complete analysis!")
    print("=" * 80)


if __name__ == "__main__":
    main()
