#!/usr/bin/env python3
"""
REINVENT4 修复版分析工具
"""
import os
import glob
import pandas as pd
import numpy as np
from datetime import datetime
import re
import warnings
warnings.filterwarnings('ignore')

class ConfigParser:
    """解析TOML配置文件"""
    
    @staticmethod
    def parse_toml(filepath):
        """解析TOML配置文件"""
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
            config['component_details'] = components
            
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
        
        # 计算指标
        self._calculate_metrics()
    
    def _calculate_metrics(self):
        """计算统计指标 - 修复DataFrame判断"""
        m = {}
        
        # 分子数量
        m['n_total'] = len(self.data.get('results', []))
        m['n_gold'] = len(self.data.get('gold', []))
        m['n_high'] = len(self.data.get('high', []))
        m['n_good'] = len(self.data.get('good', []))
        
        # 从数据中提取指标 - 修复判断
        df = self.data.get('gold')
        if df is None or len(df) == 0:  # 修复：使用 is None 和 len()
            df = self.data.get('results')
        
        if df is not None and len(df) > 0:  # 修复：正确的判断方式
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
    """生成HTML报告"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    def generate(self, output_file='runs_analysis_report.html'):
        """生成HTML报告"""
        print("\n📄 生成HTML报告...")
        
        html = self._html_header()
        html += self._section_summary()
        html += self._section_gold_comparison()
        html += self._section_config_comparison()
        html += self._section_metrics()
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
        body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; border-left: 4px solid #3498db; padding-left: 10px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th {{ background: #3498db; color: white; padding: 12px; text-align: left; }}
        td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background: #f5f5f5; }}
        .badge {{ display: inline-block; padding: 3px 8px; border-radius: 10px; font-size: 0.9em; }}
        .badge-gold {{ background: #ffd700; color: #000; }}
        .badge-high {{ background: #ff8c00; color: #fff; }}
        .highlight {{ background: #fffacd; font-weight: bold; }}
        .metric-card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                       color: white; padding: 20px; border-radius: 8px; margin: 10px; 
                       display: inline-block; min-width: 200px; }}
        .metric-value {{ font-size: 2em; font-weight: bold; }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 REINVENT4 Multi-Run Analysis Report</h1>
    <p style="text-align: center; color: #7f8c8d;">Generated: {self.timestamp}</p>
"""
    
    def _section_summary(self):
        total_gold = sum(a.metrics.get('n_gold', 0) for a in self.analyzers)
        total_high = sum(a.metrics.get('n_high', 0) for a in self.analyzers)
        
        html = f"""
    <h2>1. 执行摘要</h2>
    <div style="text-align: center;">
        <div class="metric-card">
            <div>总运行数</div>
            <div class="metric-value">{len(self.analyzers)}</div>
        </div>
        <div class="metric-card">
            <div>金标准分子</div>
            <div class="metric-value">{total_gold}</div>
        </div>
        <div class="metric-card">
            <div>高质量分子</div>
            <div class="metric-value">{total_high}</div>
        </div>
    </div>
    
    <table>
        <tr><th>Run</th><th>总分子数</th><th>Gold</th><th>High</th><th>Good</th><th>成功率</th></tr>
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
        html += "</table>"
        return html
    
    def _section_gold_comparison(self):
        return """
    <h2>2. 金标准分子对比</h2>
    <p><strong>金标准定义:</strong> Score ≥ 0.80, pIC50 ≥ 6.0, QED ≥ 0.5, MW: 250-600 Da</p>
"""
    
    def _section_config_comparison(self):
        html = """
    <h2>3. 配置对比</h2>
    <table>
        <tr><th>Run</th><th>Prior模型</th><th>批次大小</th><th>学习率</th><th>Sigma</th><th>评分类型</th></tr>
"""
        for a in self.analyzers:
            prior = os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'
            html += f"""
        <tr>
            <td><strong>{a.run_name}</strong></td>
            <td>{prior}</td>
            <td>{a.config.get('batch_size', 'N/A')}</td>
            <td>{a.config.get('learning_rate', 'N/A')}</td>
            <td>{a.config.get('sigma', 'N/A')}</td>
            <td>{a.config.get('scoring_type', 'N/A')}</td>
        </tr>
"""
        html += "</table>"
        return html
    
    def _section_metrics(self):
        html = """
    <h2>4. 性能指标</h2>
    <table>
        <tr><th>Run</th><th>平均pIC50</th><th>最大pIC50</th><th>平均QED</th><th>平均MW</th><th>平均分数</th></tr>
"""
        for a in self.analyzers:
            html += f"""
        <tr>
            <td><strong>{a.run_name}</strong></td>
            <td>{a.metrics.get('qsar_mean', 0):.2f}</td>
            <td class="highlight">{a.metrics.get('qsar_max', 0):.2f}</td>
            <td>{a.metrics.get('qed_mean', 0):.3f}</td>
            <td>{a.metrics.get('mw_mean', 0):.1f}</td>
            <td>{a.metrics.get('total_score_mean', 0):.3f}</td>
        </tr>
"""
        html += "</table>"
        return html
    
    def _section_components(self):
        all_comps = set()
        for a in self.analyzers:
            all_comps.update(a.config.get('components', []))
        
        html = """
    <h2>5. 评分组件权重</h2>
    <table>
        <tr><th>组件</th>
"""
        for a in self.analyzers:
            html += f"<th>{a.run_name}</th>"
        html += "</tr>"
        
        for comp in sorted(all_comps):
            html += f"<tr><td><strong>{comp}</strong></td>"
            for a in self.analyzers:
                weight = a.config.get('component_weights', {}).get(comp)
                html += f"<td>{weight:.2f if weight else '—'}</td>"
            html += "</tr>"
        
        html += "</table>"
        return html
    
    def _section_detailed(self):
        html = """
    <h2>6. 详细信息</h2>
"""
        for a in self.analyzers:
            html += f"""
    <h3>{a.run_name}</h3>
    <p>总分子: {a.metrics.get('n_total', 0):,} | 
       Gold: {a.metrics.get('n_gold', 0)} | 
       组件数: {len(a.config.get('components', []))}</p>
"""
        return html
    
    def _html_footer(self):
        return """
</div>
</body>
</html>
"""


def main():
    print("="*80)
    print("REINVENT4 修复版分析")
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
            import traceback
            traceback.print_exc()
    
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
                'QSAR_mean': f"{a.metrics.get('qsar_mean', 0):.2f}",
                'QSAR_max': f"{a.metrics.get('qsar_max', 0):.2f}",
                'Components': len(a.config.get('components', []))
            })
        
        pd.DataFrame(summary_data).to_csv('runs_summary.csv', index=False)
        print("✅ CSV: runs_summary.csv")
        print("\n" + pd.DataFrame(summary_data).to_string(index=False))
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()
