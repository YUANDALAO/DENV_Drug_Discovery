#!/usr/bin/env python3
"""
REINVENT4 最终版分析工具
- 修复所有错误
- 添加分子结构图
- 添加统计分布
- 响应式设计（支持手机查看）
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
warnings.filterwarnings('ignore')

class ConfigParser:
    """解析TOML配置文件"""
    
    @staticmethod
    def parse_toml(filepath):
        config = {
            'run_type': None,
            'device': None,
            'prior_file': None,
            'agent_file': None,
            'batch_size': None,
            'learning_rate': None,
            'sigma': None,
            'scoring_type': None,
            'components': [],
            'component_weights': {},
        }
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            config['run_type'] = ConfigParser._extract(content, r'run_type\s*=\s*"([^"]+)"')
            config['device'] = ConfigParser._extract(content, r'device\s*=\s*"([^"]+)"')
            config['prior_file'] = ConfigParser._extract(content, r'prior_file\s*=\s*"([^"]+)"')
            config['agent_file'] = ConfigParser._extract(content, r'agent_file\s*=\s*"([^"]+)"')
            config['batch_size'] = ConfigParser._extract(content, r'batch_size\s*=\s*(\d+)', int)
            config['learning_rate'] = ConfigParser._extract(content, r'rate\s*=\s*([0-9.]+)', float)
            config['sigma'] = ConfigParser._extract(content, r'sigma\s*=\s*([0-9.]+)', float)
            
            scoring_match = re.search(r'\[(?:stage\.)?scoring\]\s*type\s*=\s*"([^"]+)"', content)
            if scoring_match:
                config['scoring_type'] = scoring_match.group(1)
            
            components = ConfigParser._extract_components(content)
            config['components'] = [c['name'] for c in components]
            config['component_weights'] = {c['name']: c['weight'] for c in components}
            
        except Exception as e:
            print(f"    ⚠️  配置解析错误: {e}")
        
        return config
    
    @staticmethod
    def _extract(content, pattern, dtype=str):
        match = re.search(pattern, content)
        if match:
            try:
                return dtype(match.group(1))
            except:
                return match.group(1)
        return None
    
    @staticmethod
    def _extract_components(content):
        components = []
        component_pattern = r'\[\[(?:stage\.)?scoring\.component\]\]\s*\[(?:stage\.)?scoring\.component\.(\w+)\]'
        
        for match in re.finditer(component_pattern, content):
            comp_name = match.group(1)
            comp_start = match.start()
            next_comp = re.search(r'\[\[(?:stage\.)?scoring\.component\]\]', content[comp_start+10:])
            comp_end = comp_start + next_comp.start() + 10 if next_comp else len(content)
            comp_section = content[comp_start:comp_end]
            
            weight_match = re.search(r'weight\s*=\s*([0-9.]+)', comp_section)
            weight = float(weight_match.group(1)) if weight_match else 1.0
            
            name_match = re.search(r'name\s*=\s*"([^"]+)"', comp_section)
            display_name = name_match.group(1) if name_match else comp_name
            
            components.append({
                'name': comp_name,
                'display_name': display_name,
                'weight': weight
            })
        
        return components


class RunAnalyzer:
    """分析单个run"""
    
    def __init__(self, run_name, run_dir):
        self.run_name = run_name
        self.run_dir = run_dir
        self.data = {}
        self.config = {}
        self.metrics = {}
    
    def load_data(self):
        """加载数据"""
        print(f"  分析 {self.run_name}...")
        
        # 加载配置
        config_path = os.path.join(self.run_dir, 'config.toml')
        if os.path.exists(config_path):
            self.config = ConfigParser.parse_toml(config_path)
            print(f"    ✓ 配置: {len(self.config['components'])} 个组件")
        
        # 智能查找CSV文件
        csv_mapping = {
            'results': ['results_*.csv', 'results.csv'],
            'gold': ['candidates_gold.csv', 'candidates_金标准*.csv', '*gold*.csv'],
            'high': ['candidates_high.csv', 'candidates_高标准*.csv', '*high*.csv'],
            'good': ['candidates_good.csv', 'candidates_中标准*.csv', '*good*.csv'],
        }
        
        for key, patterns in csv_mapping.items():
            for pattern in patterns:
                files = glob.glob(os.path.join(self.run_dir, pattern))
                if files:
                    try:
                        df = pd.read_csv(files[0])
                        self.data[key] = df
                        print(f"    ✓ {key}: {len(df)} 条")
                        break
                    except:
                        pass
        
        self._calculate_metrics()
    
    def _calculate_metrics(self):
        """计算统计指标"""
        m = {}
        
        m['n_total'] = len(self.data.get('results', []))
        m['n_gold'] = len(self.data.get('gold', []))
        m['n_high'] = len(self.data.get('high', []))
        m['n_good'] = len(self.data.get('good', []))
        
        df = self.data.get('gold')
        if df is None or len(df) == 0:
            df = self.data.get('results')
        
        if df is not None and len(df) > 0:
            df.columns = df.columns.str.strip()
            
            col_map = {
                'total_score': ['Score', 'total_score', 'score'],
                'qsar': ['DENV_Activity (raw)', 'DENV_Activity', 'QSAR_Score', 'qsar'],
                'qed': ['QED (raw)', 'QED', 'qed'],
                'mw': ['MW', 'Molecular weight', 'molecular_weight'],
                'sa': ['SA (raw)', 'SA', 'SAScore', 'sa_score'],
                'logp': ['LogP', 'SlogP (RDKit)', 'SlogP'],
                'tpsa': ['TPSA', 'tpsa'],
            }
            
            for metric, possible_cols in col_map.items():
                for col in possible_cols:
                    if col in df.columns:
                        vals = pd.to_numeric(df[col], errors='coerce').dropna()
                        if len(vals) > 0:
                            m[f'{metric}_mean'] = vals.mean()
                            m[f'{metric}_std'] = vals.std()
                            m[f'{metric}_min'] = vals.min()
                            m[f'{metric}_max'] = vals.max()
                        break
        
        self.metrics = m


class HTMLReportGenerator:
    """生成HTML报告（移动端友好）"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    def generate(self, output_file='runs_analysis_report.html'):
        """生成HTML报告"""
        print("\n📄 生成HTML报告...")
        
        html = self._html_header()
        html += self._section_summary()
        html += self._section_gold_comparison()
        html += self._section_structures()
        html += self._section_config_comparison()
        html += self._section_metrics_detailed()
        html += self._section_components()
        html += self._section_detailed()
        html += self._html_footer()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"✅ HTML报告: {output_file}")
        return output_file
    
    def _html_header(self):
        return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>REINVENT4 Runs Analysis</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            line-height: 1.6; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 10px;
        }}
        .container {{ 
            max-width: 1400px; 
            margin: 0 auto; 
            background: white; 
            padding: 20px; 
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
        }}
        h1 {{ 
            color: #2c3e50; 
            border-bottom: 3px solid #3498db; 
            padding-bottom: 10px;
            font-size: 1.8em;
        }}
        h2 {{ 
            color: #34495e; 
            margin-top: 30px; 
            border-left: 4px solid #3498db; 
            padding-left: 10px;
            font-size: 1.4em;
        }}
        h3 {{ color: #555; margin-top: 20px; font-size: 1.2em; }}
        
        table {{ 
            width: 100%; 
            border-collapse: collapse; 
            margin: 20px 0;
            overflow-x: auto;
            display: block;
        }}
        thead {{ display: table-header-group; }}
        tbody {{ display: table-row-group; }}
        th {{ 
            background: #3498db; 
            color: white; 
            padding: 10px; 
            text-align: left;
            position: sticky;
            top: 0;
            font-size: 0.9em;
        }}
        td {{ 
            padding: 8px; 
            border-bottom: 1px solid #ddd;
            font-size: 0.85em;
        }}
        tr:hover {{ background: #f5f5f5; }}
        
        .metric-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .metric-card {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            padding: 20px; 
            border-radius: 8px;
            text-align: center;
        }}
        .metric-value {{ 
            font-size: 2.5em; 
            font-weight: bold;
            margin: 10px 0;
        }}
        .metric-label {{ 
            font-size: 0.9em;
            opacity: 0.9;
        }}
        
        .stats-box {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
        }}
        .stats-box strong {{ color: #3498db; }}
        
        .highlight {{ background: #fffacd; font-weight: bold; }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 10px;
            font-size: 0.85em;
            font-weight: 600;
        }}
        .badge-gold {{ background: #ffd700; color: #000; }}
        .badge-high {{ background: #ff8c00; color: #fff; }}
        
        .structure-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }}
        .structure-card {{
            border: 2px solid #ddd;
            border-radius: 8px;
            padding: 10px;
            background: #f9f9f9;
        }}
        .structure-card h4 {{
            margin: 0 0 10px 0;
            color: #3498db;
            font-size: 1em;
        }}
        
        @media (max-width: 768px) {{
            .container {{ padding: 15px; }}
            h1 {{ font-size: 1.5em; }}
            h2 {{ font-size: 1.2em; }}
            table {{ font-size: 0.8em; }}
            .metric-grid {{ grid-template-columns: 1fr; }}
            th, td {{ padding: 6px; }}
        }}
        
        .info-box {{
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 15px;
            margin: 15px 0;
            border-radius: 4px;
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 REINVENT4 Multi-Run Analysis Report</h1>
    <p style="text-align: center; color: #7f8c8d; margin-bottom: 20px;">
        Generated: {self.timestamp} | Optimized for Mobile & Desktop
    </p>
"""
    
    def _section_summary(self):
        total_gold = sum(a.metrics.get('n_gold', 0) for a in self.analyzers)
        total_high = sum(a.metrics.get('n_high', 0) for a in self.analyzers)
        total_good = sum(a.metrics.get('n_good', 0) for a in self.analyzers)
        
        html = f"""
    <h2>1. 📊 执行摘要</h2>
    <div class="metric-grid">
        <div class="metric-card">
            <div class="metric-label">总运行数</div>
            <div class="metric-value">{len(self.analyzers)}</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">金标准分子</div>
            <div class="metric-value">{total_gold}</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">高质量分子</div>
            <div class="metric-value">{total_high}</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">良好分子</div>
            <div class="metric-value">{total_good}</div>
        </div>
    </div>
    
    <h3>Runs概览</h3>
    <table>
        <thead>
            <tr><th>Run</th><th>总数</th><th>Gold</th><th>High</th><th>Good</th><th>成功率</th></tr>
        </thead>
        <tbody>
"""
        for a in self.analyzers:
            n_total = a.metrics.get('n_total', 1)
            n_gold = a.metrics.get('n_gold', 0)
            rate = (n_gold/n_total*100) if n_total > 0 else 0
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{n_total:,}</td>
                <td class="highlight">{n_gold}</td>
                <td>{a.metrics.get('n_high', 0)}</td>
                <td>{a.metrics.get('n_good', 0)}</td>
                <td>{rate:.2f}%</td>
            </tr>
"""
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_gold_comparison(self):
        return """
    <h2>2. 🥇 金标准分子定义</h2>
    <div class="info-box">
        <strong>金标准分子必须同时满足:</strong><br>
        ✓ QSAR预测活性 (pIC50) ≥ 8.0 (IC50 ≤ 1 μM)<br>
        ✓ 类药性 (QED) ≥ 0.7<br>
        ✓ 分子量: 300-500 Da<br>
        ✓ 合成可及性 (SA Score) ≤ 4.0<br>
        ✓ LogP: 1-4
    </div>
"""
    
    def _section_structures(self):
        """分子结构对比"""
        html = """
    <h2>3. 🔬 Top金标准分子结构对比</h2>
    <div class="structure-grid">
"""
        
        for a in self.analyzers:
            gold_df = a.data.get('gold')
            if gold_df is not None and len(gold_df) > 0:
                top3 = gold_df.head(3)
                
                html += f"""
        <div class="structure-card">
            <h4>{a.run_name} (Top 3)</h4>
"""
                
                for idx, row in top3.iterrows():
                    smiles = row.get('SMILES', 'N/A')
                    score = row.get('Score', row.get('total_score', 'N/A'))
                    qsar = row.get('DENV_Activity (raw)', row.get('DENV_Activity', 'N/A'))
                    
                    # 尝试生成分子图
                    img_html = self._generate_molecule_image(smiles, a.run_name, idx)
                    
                    score_str = f"{score:.3f}" if isinstance(score, (int, float)) else str(score)
                    qsar_str = f"{qsar:.2f}" if isinstance(qsar, (int, float)) else str(qsar)
                    
                    html += f"""
            <div style="margin: 10px 0; padding: 10px; background: white; border-radius: 4px;">
                <strong>#{idx+1}</strong><br>
                {img_html}
                <small style="display: block; margin-top: 5px; word-break: break-all;">
                    SMILES: {smiles[:50]}{'...' if len(smiles) > 50 else ''}<br>
                    Score: {score_str} | pIC50: {qsar_str}
                </small>
            </div>
"""
                
                html += """
        </div>
"""
        
        html += """
    </div>
"""
        return html
    
    def _generate_molecule_image(self, smiles, run_name, idx):
        """生成分子图片（如果RDKit可用）"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(250, 250))
                
                # 转换为base64
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                
                return f'<img src="data:image/png;base64,{img_str}" style="width: 100%; max-width: 250px; height: auto;">'
        except:
            pass
        
        return '<div style="background: #eee; padding: 20px; text-align: center; color: #999;">结构图需要RDKit</div>'
    
    def _section_config_comparison(self):
        html = """
    <h2>4. ⚙️ 配置对比</h2>
    <table>
        <thead>
            <tr><th>Run</th><th>Prior模型</th><th>批次</th><th>学习率</th><th>Sigma</th><th>评分</th><th>组件数</th></tr>
        </thead>
        <tbody>
"""
        for a in self.analyzers:
            prior = os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'
            lr = a.config.get('learning_rate')
            lr_str = f"{lr:.4f}" if lr is not None else 'N/A'
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{prior}</td>
                <td>{a.config.get('batch_size', 'N/A')}</td>
                <td>{lr_str}</td>
                <td>{a.config.get('sigma', 'N/A')}</td>
                <td>{a.config.get('scoring_type', 'N/A')}</td>
                <td>{len(a.config.get('components', []))}</td>
            </tr>
"""
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_metrics_detailed(self):
        """详细性能指标（带统计分布）"""
        html = """
    <h2>5. 📈 详细性能指标 (Mean ± SD)</h2>
    
    <h3>活性预测 (QSAR pIC50)</h3>
    <table>
        <thead>
            <tr><th>Run</th><th>均值 ± 标准差</th><th>最小值</th><th>最大值</th></tr>
        </thead>
        <tbody>
"""
        for a in self.analyzers:
            mean = a.metrics.get('qsar_mean', 0)
            std = a.metrics.get('qsar_std', 0)
            min_val = a.metrics.get('qsar_min', 0)
            max_val = a.metrics.get('qsar_max', 0)
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{mean:.2f} ± {std:.2f}</td>
                <td>{min_val:.2f}</td>
                <td class="highlight">{max_val:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    
    <h3>类药性指标</h3>
    <table>
        <thead>
            <tr><th>Run</th><th>QED</th><th>分子量 (Da)</th><th>LogP</th><th>SA Score</th></tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            qed_m = a.metrics.get('qed_mean', 0)
            qed_s = a.metrics.get('qed_std', 0)
            mw_m = a.metrics.get('mw_mean', 0)
            mw_s = a.metrics.get('mw_std', 0)
            logp_m = a.metrics.get('logp_mean', 0)
            logp_s = a.metrics.get('logp_std', 0)
            sa_m = a.metrics.get('sa_mean', 0)
            sa_s = a.metrics.get('sa_std', 0)
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{qed_m:.2f} ± {qed_s:.2f}</td>
                <td>{mw_m:.1f} ± {mw_s:.1f}</td>
                <td>{logp_m:.2f} ± {logp_s:.2f}</td>
                <td>{sa_m:.2f} ± {sa_s:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_components(self):
        """评分组件权重对比"""
        all_comps = set()
        for a in self.analyzers:
            all_comps.update(a.config.get('components', []))
        
        html = """
    <h2>6. 🎯 评分组件权重对比</h2>
    <div style="overflow-x: auto;">
    <table>
        <thead>
            <tr><th>组件</th>
"""
        for a in self.analyzers:
            html += f"<th>{a.run_name}</th>"
        html += """
            </tr>
        </thead>
        <tbody>
"""
        
        for comp in sorted(all_comps):
            html += f"<tr><td><strong>{comp}</strong></td>"
            for a in self.analyzers:
                weight = a.config.get('component_weights', {}).get(comp)
                if weight is not None:
                    html += f"<td>{weight:.2f}</td>"
                else:
                    html += "<td style='color: #ccc;'>—</td>"
            html += "</tr>"
        
        html += """
        </tbody>
    </table>
    </div>
"""
        return html
    
    def _section_detailed(self):
        """详细信息"""
        html = """
    <h2>7. 📋 详细Run信息</h2>
"""
        for a in self.analyzers:
            html += f"""
    <div class="stats-box">
        <h3>{a.run_name}</h3>
        <p>
            <strong>总分子数:</strong> {a.metrics.get('n_total', 0):,} | 
            <strong>Gold:</strong> {a.metrics.get('n_gold', 0)} | 
            <strong>High:</strong> {a.metrics.get('n_high', 0)} | 
            <strong>Good:</strong> {a.metrics.get('n_good', 0)}<br>
            <strong>组件数:</strong> {len(a.config.get('components', []))} | 
            <strong>Prior:</strong> {os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'}
        </p>
    </div>
"""
        return html
    
    def _html_footer(self):
        return """
    <div style="margin-top: 50px; padding-top: 20px; border-top: 2px solid #ddd; text-align: center; color: #7f8c8d;">
        <p><strong>REINVENT4 Comprehensive Analysis Report</strong></p>
        <p style="font-size: 0.9em;">Generated by CHEN, Yuan | Optimized for Mobile & Desktop</p>
    </div>
</div>
</body>
</html>
"""


def main():
    print("="*80)
    print("REINVENT4 最终版分析")
    print("="*80)
    
    runs_dir = "experiments/runs"
    run_dirs = {}
    
    for item in os.listdir(runs_dir):
        if item.lower() == 'archive':
            continue
        item_path = os.path.join(runs_dir, item)
        if os.path.isdir(item_path):
            run_dirs[item] = item_path
    
    print(f"\n找到 {len(run_dirs)} 个runs")
    
    analyzers = []
    for run_name, run_path in sorted(run_dirs.items()):
        analyzer = RunAnalyzer(run_name, run_path)
        try:
            analyzer.load_data()
            analyzers.append(analyzer)
        except Exception as e:
            print(f"  ❌ {run_name} 失败: {e}")
    
    print(f"\n✅ 成功分析 {len(analyzers)} 个runs")
    
    if len(analyzers) > 0:
        reporter = HTMLReportGenerator(analyzers)
        reporter.generate()
        
        # CSV导出
        summary_data = []
        for a in analyzers:
            summary_data.append({
                'Run': a.run_name,
                'Total': a.metrics.get('n_total', 0),
                'Gold': a.metrics.get('n_gold', 0),
                'High': a.metrics.get('n_high', 0),
                'QSAR_mean±std': f"{a.metrics.get('qsar_mean', 0):.2f}±{a.metrics.get('qsar_std', 0):.2f}",
                'QSAR_max': f"{a.metrics.get('qsar_max', 0):.2f}",
                'QED_mean±std': f"{a.metrics.get('qed_mean', 0):.2f}±{a.metrics.get('qed_std', 0):.2f}",
                'MW_mean±std': f"{a.metrics.get('mw_mean', 0):.1f}±{a.metrics.get('mw_std', 0):.1f}",
                'Components': len(a.config.get('components', []))
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv('runs_summary.csv', index=False)
        print("✅ CSV: runs_summary.csv")
        print("\n" + df.to_string(index=False))
    
    print("\n" + "="*80)
    print("📱 HTML报告已优化，支持手机和电脑查看")
    print("🔬 包含分子结构图（需要RDKit）和完整统计分布")
    print("="*80)

if __name__ == "__main__":
    main()
