#!/usr/bin/env python3
"""
REINVENT4 终极HTML报告生成器 (最终修复版)
- 自动发现并比较所有runs
- 修复数据筛选与隔离问题
- 增强可视化与统计，兼容不同RDKit版本
- 移除特殊字符以保证最大兼容性
"""
import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime
import re
import warnings
import base64
from io import BytesIO
import matplotlib

# 使用 'Agg' 后端，这样在没有图形界面的服务器上也能运行
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# --- 全局配置 ---
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 120

# --- 数据加载与解析模块 ---

class ConfigParser:
    """解析TOML配置文件 (简化版)"""
    @staticmethod
    def parse_toml(filepath):
        config = {'components': []}
        try:
            with open(filepath, 'r', encoding='utf-8') as f: content = f.read()
            components = re.findall(r'\[stage\.scoring\.component\.(\w+)\]', content)
            config['components'] = components
        except Exception: pass
        return config

class MolecularFeatureAnalyzer:
    """分子特征分析器 (智能兼容版)"""
    @staticmethod
    def calculate_features(df):
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen
            from rdkit.Chem.Scaffolds import MurckoScaffold
        except ImportError: 
            print("    - Warning: RDKit not found, skipping feature calculation.")
            return None
        
        # --- 智能版本检测FractionCsp3 ---
        if hasattr(Lipinski, 'FractionCsp3'):
            frac_csp3_func = Lipinski.FractionCsp3
        elif hasattr(Descriptors, 'FractionCsp3'):
            frac_csp3_func = Descriptors.FractionCsp3
        else:
            frac_csp3_func = lambda x: np.nan
            print("    - Warning: FractionCsp3 not found in this RDKit version. Skipping.")

        features = []
        for _, row in df.iterrows():
            smiles = row.get('SMILES', '')
            if not smiles: continue
            mol = Chem.MolFromSmiles(smiles)
            if not mol: continue
            
            feat = {
                'SMILES': smiles,
                'MW': Descriptors.MolWt(mol), 'LogP': Crippen.MolLogP(mol), 'TPSA': Descriptors.TPSA(mol),
                'HBA': Lipinski.NumHAcceptors(mol), 'HBD': Lipinski.NumHDonors(mol), 'RotBonds': Descriptors.NumRotatableBonds(mol),
                'AromaticRings': Descriptors.NumAromaticRings(mol), 'FractionCSP3': frac_csp3_func(mol),
                'NumRings': Descriptors.RingCount(mol),
            }
            
            for col in ['DENV_Activity (raw)', 'DENV_Activity', 'pIC50']:
                if col in row: feat['pIC50'] = row[col]; break
            for col in ['QED (raw)', 'QED']:
                if col in row: feat['QED'] = row[col]; break
            for col in ['SA (raw)', 'SA']:
                if col in row: feat['SA'] = row[col]; break
            
            try: feat['Scaffold'] = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
            except: feat['Scaffold'] = None
            
            features.append(feat)
        return pd.DataFrame(features) if features else None

class RunAnalyzer:
    """分析单个run, 修复了数据隔离问题"""
    def __init__(self, run_name, run_dir):
        # 确保每个实例都有自己独立的数据结构
        self.run_name = run_name
        self.run_dir = run_dir
        self.data = {}
        self.config = {}
        self.metrics = {}
        self.molecular_features = None

    def analyze(self):
        print(f"  Analyzing {self.run_name}...")
        
        config_path_toml = os.path.join(self.run_dir, 'config.toml')
        config_path_json = glob.glob(os.path.join(self.run_dir, '*.json'))
        config_path = config_path_toml if os.path.exists(config_path_toml) else (config_path_json[0] if config_path_json else None)

        if config_path:
            self.config = ConfigParser.parse_toml(config_path)
            print(f"    + Config: {len(self.config['components'])} components found")

        csv_mapping = {
            'results': ['results_*.csv', 'results.csv'],
            'gold': ['candidates_gold.csv', 'candidates_金标准*.csv', '*gold*.csv', 'promising*.csv'],
        }
        
        for key, patterns in csv_mapping.items():
            for pattern in patterns:
                files = glob.glob(os.path.join(self.run_dir, pattern))
                if files:
                    try:
                        self.data[key] = pd.read_csv(files[0])
                        self.data[key].columns = self.data[key].columns.str.strip() # 清理列名空格
                        print(f"    + {key.capitalize()}: {len(self.data[key])} molecules loaded")
                        break
                    except Exception as e:
                        print(f"    - Error reading {files[0]}: {e}")
        
        self._analyze_features()

    def _analyze_features(self):
        df_gold = self.data.get('gold')
        if df_gold is None or df_gold.empty:
            print("    - 'gold' file not found or empty, skipping feature analysis.")
            self.metrics['Total'] = len(self.data.get('results', []))
            self.metrics['Gold'] = 0
            self.metrics['Success_Rate_%'] = 0
            return

        self.molecular_features = MolecularFeatureAnalyzer.calculate_features(df_gold)
        
        if self.molecular_features is not None and not self.molecular_features.empty:
            print(f"    + Feature Analysis: {len(self.molecular_features)} 'gold' molecules analyzed")
            for feat in ['pIC50', 'QED', 'SA', 'MW', 'LogP', 'TPSA', 'RotBonds', 'AromaticRings', 'FractionCSP3']:
                if feat in self.molecular_features.columns and self.molecular_features[feat].notna().any():
                    self.metrics[f'{feat}_mean'] = self.molecular_features[feat].mean()
            if 'Scaffold' in self.molecular_features.columns:
                self.metrics['Scaffolds'] = self.molecular_features['Scaffold'].nunique()

        self.metrics['Total'] = len(self.data.get('results', []))
        self.metrics['Gold'] = len(df_gold)
        self.metrics['Success_Rate_%'] = (self.metrics['Gold'] / self.metrics['Total'] * 100) if self.metrics['Total'] > 0 else 0

class HTMLReportGenerator:
    """生成高度美化的多Run对比HTML报告"""
    def __init__(self, analyzers):
        self.analyzers = sorted(analyzers, key=lambda a: a.metrics.get('Gold', 0), reverse=True)
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    def generate(self, output_file='ultimate_runs_report.html'):
        print("\nGenerating HTML report...")
        html = self._html_header()
        html += self._section_summary_table()
        html += self._section_visualizations()
        html += self._section_top_structures()
        html += self._html_footer()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"+ HTML report generated: {output_file}")
    
    def _html_header(self):
        return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8"><title>REINVENT4 Multi-Run Analysis</title>
    <style>
        body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; background-color: #f8f9fa; color: #212529; margin: 0; padding: 2rem; }}
        .container {{ max-width: 1800px; margin: auto; background: #fff; border-radius: 12px; box-shadow: 0 8px 30px rgba(0,0,0,0.1); padding: 2rem; }}
        h1, h2 {{ color: #0056b3; border-bottom: 2px solid #dee2e6; padding-bottom: 0.5rem; margin-bottom: 1.5rem; }}
        h1 {{ font-size: 2rem; }} h2 {{ font-size: 1.5rem; margin-top: 2.5rem; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 1rem; font-size: 0.9rem; }}
        th, td {{ padding: 0.75rem; text-align: left; border-bottom: 1px solid #dee2e6; }}
        th {{ background-color: #e9ecef; font-weight: 600; }}
        tbody tr:hover {{ background-color: #f1f3f5; }}
        .best-run td {{ background-color: #d4edda; font-weight: bold; }}
        .viz-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(500px, 1fr)); gap: 2rem; margin-top: 1.5rem; }}
        .plot-container {{ padding: 1rem; border: 1px solid #dee2e6; border-radius: 8px; background: #fff; }}
        img {{ max-width: 100%; height: auto; border: 1px solid #eee; }}
        .struct-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1.5rem; margin-top: 1rem; }}
        .struct-card {{ border: 1px solid #dee2e6; border-radius: 8px; padding: 1rem; text-align: center; background: #fdfdfd; }}
        .struct-card h4 {{ margin: 0 0 1rem 0; font-size: 1.1rem; color: #0056b3; }}
        .struct-card img {{ margin-bottom: 0.5rem; }}
        .struct-card p {{ font-size: 0.8rem; color: #495057; word-break: break-all; }}
        .struct-card-item {{ margin-top: 1rem; padding-top: 1rem; border-top: 1px dashed #ccc; }}
    </style>
</head>
<body><div class="container">
<h1>REINVENT4 Multi-Run Performance Analysis</h1>
<p>Generated on: {self.timestamp}</p>
"""
    
    def _section_summary_table(self):
        html = """
    <h2>Overall Performance Ranking (by 'Gold' Molecules Count)</h2>
    <table><thead>
        <tr><th>Rank</th><th>Run Name</th><th>Gold Molecules</th><th>Success Rate (%)</th><th>Unique Scaffolds</th><th>Avg pIC50</th><th>Avg QED</th><th>Avg SA Score</th><th>Avg MW</th></tr>
    </thead><tbody>
"""
        for i, a in enumerate(self.analyzers):
            row_class = 'best-run' if i == 0 and a.metrics.get('Gold', 0) > 0 else ''
            html += f"""
        <tr class="{row_class}">
            <td>{i+1}</td>
            <td><strong>{a.run_name}</strong></td>
            <td>{a.metrics.get('Gold', 0)}</td>
            <td>{a.metrics.get('Success_Rate_%', 0):.4f}</td>
            <td>{a.metrics.get('Scaffolds', 'N/A')}</td>
            <td>{a.metrics.get('pIC50_mean', float('nan')):.2f}</td>
            <td>{a.metrics.get('QED_mean', float('nan')):.2f}</td>
            <td>{a.metrics.get('SA_mean', float('nan')):.2f}</td>
            <td>{a.metrics.get('MW_mean', float('nan')):.1f}</td>
        </tr>"""
        html += "</tbody></table>"
        return html

    def _plot_to_base64(self, plot_func, **kwargs):
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_func(ax=ax, **kwargs)
        buf = BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight')
        plt.close(fig)
        return base64.b64encode(buf.getvalue()).decode('utf-8')

    def _section_visualizations(self):
        html = "<h2>Comparative Visualizations</h2><div class='viz-grid'>"
        
        def plot_gold(ax):
            names = [a.run_name for a in self.analyzers]
            golds = [a.metrics.get('Gold', 0) for a in self.analyzers]
            colors = sns.color_palette("viridis", len(names))
            bars = ax.bar(names, golds, color=colors)
            ax.set_title("Gold Standard Molecules per Run", fontsize=14, fontweight='bold')
            ax.set_ylabel("Count")
            ax.tick_params(axis='x', rotation=45, ha='right')
            ax.bar_label(bars, padding=3)
            ax.grid(axis='y', linestyle='--', alpha=0.7)

        html += f"<div class='plot-container'><img src='data:image/png;base64,{self._plot_to_base64(plot_gold)}'></div>"
        
        def plot_scatter(ax):
            colors = sns.color_palette("plasma", len(self.analyzers))
            for i, a in enumerate(self.analyzers):
                if a.molecular_features is not None and not a.molecular_features.empty:
                    if 'pIC50' in a.molecular_features.columns and 'QED' in a.molecular_features.columns:
                        ax.scatter(a.molecular_features['pIC50'], a.molecular_features['QED'], label=a.run_name, color=colors[i], alpha=0.6, edgecolors='k', s=50)
            ax.set_title("pIC50 vs QED of Gold Molecules", fontsize=14, fontweight='bold')
            ax.set_xlabel("Predicted pIC50")
            ax.set_ylabel("QED Score")
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(linestyle='--', alpha=0.7)

        html += f"<div class='plot-container'><img src='data:image/png;base64,{self._plot_to_base64(plot_scatter)}'></div>"
        
        html += "</div>"
        return html

    def _section_top_structures(self):
        html = "<h2>Top 3 'Gold' Molecules from Each Run</h2>"
        html += "<div class='struct-grid'>"
        for a in self.analyzers:
            html += f"<div class='struct-card'><h4>{a.run_name}</h4>"
            df_gold = a.data.get('gold')
            if df_gold is not None and not df_gold.empty:
                pic50_col = next((col for col in ['DENV_Activity (raw)', 'pIC50'] if col in df_gold.columns), None)
                top3 = df_gold.sort_values(by=pic50_col, ascending=False).head(3) if pic50_col else df_gold.head(3)

                for i, row in top3.iterrows():
                    img_html = self._generate_mol_image(row['SMILES'])
                    pic50 = row.get(pic50_col, np.nan)
                    qed = row.get('QED (raw)', np.nan)
                    sa = row.get('SA (raw)', np.nan)
                    html += f"<div class='struct-card-item'><strong>#{i+1}</strong><br>{img_html}"
                    html += f"<p>pIC50: {pic50:.2f} | QED: {qed:.2f} | SA: {sa:.2f}</p></div>"
            else:
                html += "<p style='color: #999;'>No 'gold' molecules found.</p>"
            html += "</div>"
        html += "</div>"
        return html
        
    def _generate_mol_image(self, smiles):
        try:
            from rdkit import Chem, Draw
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(250, 250))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                return f'<img src="data:image/png;base64,{img_str}">'
        except: pass
        return "<p>Image render error</p>"
        
    def _html_footer(self):
        return "</div></body></html>"


def main():
    print("================================================================================")
    print("REINVENT4 Ultimate HTML Report Generator")
    print("================================================================================")
    
    base_dir = "experiments/runs"
    if not os.path.isdir(base_dir):
        print(f"Error: Directory '{base_dir}' not found. Run this script from the project root.")
        return

    run_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)) and d.lower() != 'archive']
    
    print(f"\nFound {len(run_dirs)} runs to analyze: {', '.join(sorted(run_dirs))}")
    
    analyzers = []
    for run_name in sorted(run_dirs):
        analyzer = RunAnalyzer(run_name, os.path.join(base_dir, run_name))
        analyzer.analyze()
        analyzers.append(analyzer)
    
    print(f"\nSuccessfully analyzed {len(analyzers)} runs.")
    
    if analyzers:
        reporter = HTMLReportGenerator(analyzers)
        reporter.generate()
        
        # --- CSV Summary Generation ---
        summary_data = []
        for a in analyzers:
            row_data = {'Run': a.run_name}
            row_data.update(a.metrics)
            summary_data.append(row_data)
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv('ultimate_analysis_summary.csv', index=False)
        print(f"+ CSV summary saved to: ultimate_analysis_summary.csv")

        print("\n" + "="*50 + " QUICK SUMMARY " + "="*50)
        print(summary_df[['Run', 'Gold', 'Success_Rate_%', 'pIC50_mean', 'QED_mean', 'SA_mean']].to_string(index=False))
        print("="*115)
    
    print("\nAnalysis complete. Open 'ultimate_runs_report.html' in your browser.")

if __name__ == "__main__":
    try:
        from rdkit import Chem
    except ImportError:
        print("Error: RDKit is not installed in the current environment.")
        print("Please activate the correct conda environment (e.g., 'conda activate drug_discovery')")
    else:
        main()