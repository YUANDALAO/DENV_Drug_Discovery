#!/usr/bin/env python3
"""
REINVENT4 改进的综合分析工具
- 正确解析config.toml
- 生成HTML报告
- 排除archive文件夹
- 修复字体和编码问题
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from pathlib import Path
import json
import re
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# 设置matplotlib支持中文和避免乱码
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 9

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
            
            # 基本配置
            config['run_type'] = ConfigParser._extract(content, r'run_type\s*=\s*"([^"]+)"')
            config['device'] = ConfigParser._extract(content, r'device\s*=\s*"([^"]+)"')
            config['prior_file'] = ConfigParser._extract(content, r'prior_file\s*=\s*"([^"]+)"')
            config['agent_file'] = ConfigParser._extract(content, r'agent_file\s*=\s*"([^"]+)"')
            
            # 数值参数
            config['batch_size'] = ConfigParser._extract(content, r'batch_size\s*=\s*(\d+)', int)
            config['learning_rate'] = ConfigParser._extract(content, r'rate\s*=\s*([0-9.]+)', float)
            config['sigma'] = ConfigParser._extract(content, r'sigma\s*=\s*([0-9.]+)', float)
            
            # Scoring类型
            scoring_match = re.search(r'\[(?:stage\.)?scoring\]\s*type\s*=\s*"([^"]+)"', content)
            if scoring_match:
                config['scoring_type'] = scoring_match.group(1)
            
            # 提取所有scoring components
            components = ConfigParser._extract_components(content)
            config['components'] = [c['name'] for c in components]
            config['component_weights'] = {c['name']: c['weight'] for c in components}
            config['component_details'] = components
            
        except Exception as e:
            print(f"    ⚠️  配置解析错误: {e}")
        
        return config
    
    @staticmethod
    def _extract(content, pattern, dtype=str):
        """提取单个值"""
        match = re.search(pattern, content)
        if match:
            try:
                return dtype(match.group(1))
            except:
                return match.group(1)
        return None
    
    @staticmethod
    def _extract_components(content):
        """提取所有scoring组件及其权重"""
        components = []
        
        # 匹配所有component块
        # 模式：[[stage.scoring.component]] 或 [[scoring.component]]
        component_pattern = r'\[\[(?:stage\.)?scoring\.component\]\]\s*\[(?:stage\.)?scoring\.component\.(\w+)\]'
        
        for match in re.finditer(component_pattern, content):
            comp_name = match.group(1)
            comp_start = match.start()
            
            # 找到这个component的结束位置（下一个component或文件结尾）
            next_comp = re.search(r'\[\[(?:stage\.)?scoring\.component\]\]', content[comp_start+10:])
            comp_end = comp_start + next_comp.start() + 10 if next_comp else len(content)
            
            comp_section = content[comp_start:comp_end]
            
            # 提取weight
            weight_match = re.search(r'weight\s*=\s*([0-9.]+)', comp_section)
            weight = float(weight_match.group(1)) if weight_match else 1.0
            
            # 提取name（用户定义的名称）
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
            'results': ['results_*.csv', 'results.csv', 'output*.csv'],
            'gold': ['candidates_gold.csv', '*gold*.csv'],
            'high': ['candidates_high.csv', '*high*.csv'],
            'good': ['candidates_good.csv', '*good*.csv'],
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
        """计算统计指标"""
        m = {}
        
        # 分子数量
        m['n_total'] = len(self.data.get('results', []))
        m['n_gold'] = len(self.data.get('gold', []))
        m['n_high'] = len(self.data.get('high', []))
        m['n_good'] = len(self.data.get('good', []))
        
        # 从数据中提取指标
        df = self.data.get('gold') or self.data.get('results')
        
        if df is not None and len(df) > 0:
            # 标准化列名
            df.columns = df.columns.str.strip()
            
            # 定义列名映射
            col_map = {
                'total_score': ['total_score', 'Score', 'score'],
                'qsar': ['DENV_Activity', 'QSAR_Score', 'qsar', 'pIC50'],
                'qed': ['QED', 'qed'],
                'mw': ['MW', 'Molecular weight', 'molecular_weight'],
                'sa': ['SA', 'SAScore', 'sa_score'],
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
        """生成完整HTML报告"""
        print("\n📄 生成HTML报告...")
        
        html = self._create_html_template()
        
        # 添加各个部分
        html += self._section_summary()
        html += self._section_gold_comparison()
        html += self._section_config_comparison()
        html += self._section_metrics_comparison()
        html += self._section_component_weights()
        html += self._section_detailed_runs()
        
        html += self._html_footer()
        
        # 写入文件
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"✅ HTML报告已生成: {output_file}")
        return output_file
    
    def _create_html_template(self):
        """创建HTML模板"""
        return f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>REINVENT4 Runs Analysis Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            padding: 40px;
            border-radius: 15px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
        }}
        
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 10px;
            font-size: 2.5em;
            border-bottom: 4px solid #3498db;
            padding-bottom: 15px;
        }}
        
        h2 {{
            color: #34495e;
            margin-top: 40px;
            margin-bottom: 20px;
            padding-left: 15px;
            border-left: 5px solid #3498db;
            font-size: 1.8em;
        }}
        
        h3 {{
            color: #555;
            margin-top: 25px;
            margin-bottom: 15px;
            font-size: 1.3em;
        }}
        
        .subtitle {{
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 10px;
            font-size: 1.1em;
        }}
        
        .timestamp {{
            text-align: center;
            color: #95a5a6;
            margin-bottom: 40px;
            font-style: italic;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            font-size: 0.95em;
        }}
        
        th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 14px;
            text-align: left;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            font-size: 0.85em;
        }}
        
        td {{
            padding: 12px;
            border-bottom: 1px solid #ecf0f1;
        }}
        
        tr:hover {{
            background-color: #f8f9fa;
        }}
        
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        
        .metric-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        .metric-card {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
            transition: transform 0.3s;
        }}
        
        .metric-card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 5px 20px rgba(0,0,0,0.15);
        }}
        
        .metric-card h4 {{
            color: #2c3e50;
            margin-bottom: 10px;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #3498db;
        }}
        
        .metric-label {{
            font-size: 0.85em;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        
        .highlight {{
            background-color: #fff3cd;
            font-weight: bold;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
            margin: 2px;
        }}
        
        .badge-gold {{
            background-color: #ffd700;
            color: #000;
        }}
        
        .badge-high {{
            background-color: #ff8c00;
            color: #fff;
        }}
        
        .badge-good {{
            background-color: #4CAF50;
            color: #fff;
        }}
        
        .component-list {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin: 10px 0;
        }}
        
        .component-tag {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 6px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 500;
        }}
        
        .config-diff {{
            background-color: #ffe6e6;
            padding: 3px 6px;
            border-radius: 4px;
        }}
        
        .toc {{
            background-color: #f8f9fa;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 40px;
            border-left: 4px solid #3498db;
        }}
        
        .toc h3 {{
            margin-top: 0;
            color: #2c3e50;
        }}
        
        .toc ul {{
            list-style-type: none;
            padding-left: 20px;
        }}
        
        .toc li {{
            margin: 8px 0;
        }}
        
        .toc a {{
            color: #3498db;
            text-decoration: none;
            transition: color 0.3s;
        }}
        
        .toc a:hover {{
            color: #2980b9;
            text-decoration: underline;
        }}
        
        .warning {{
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
        }}
        
        .success {{
            background-color: #d4edda;
            border-left: 4px solid #28a745;
            padding: 15px;
            margin: 20px 0;
            border-radius: 5px;
        }}
        
        img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 3px 10px rgba(0,0,0,0.1);
            margin: 20px 0;
        }}
        
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            .container {{
                box-shadow: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 REINVENT4 Multi-Run Analysis Report</h1>
        <div class="subtitle">De Novo Drug Design for DENV NS2B-NS3 Protease Inhibitors</div>
        <div class="timestamp">Generated: {self.timestamp}</div>
        
        <div class="toc">
            <h3>📑 Table of Contents</h3>
            <ul>
                <li><a href="#summary">1. Executive Summary</a></li>
                <li><a href="#gold">2. Gold Standard Molecules Comparison</a></li>
                <li><a href="#config">3. Configuration Analysis</a></li>
                <li><a href="#metrics">4. Performance Metrics</a></li>
                <li><a href="#components">5. Scoring Component Weights</a></li>
                <li><a href="#details">6. Detailed Run Information</a></li>
            </ul>
        </div>
"""
    
    def _section_summary(self):
        """执行摘要"""
        total_gold = sum(a.metrics.get('n_gold', 0) for a in self.analyzers)
        total_high = sum(a.metrics.get('n_high', 0) for a in self.analyzers)
        total_good = sum(a.metrics.get('n_good', 0) for a in self.analyzers)
        
        html = f"""
        <h2 id="summary">1. Executive Summary</h2>
        
        <div class="metric-grid">
            <div class="metric-card">
                <h4>Total Runs Analyzed</h4>
                <div class="metric-value">{len(self.analyzers)}</div>
                <div class="metric-label">Comparison Groups</div>
            </div>
            <div class="metric-card">
                <h4>Gold Candidates</h4>
                <div class="metric-value" style="color: #ffd700;">{total_gold}</div>
                <div class="metric-label">High Quality Molecules</div>
            </div>
            <div class="metric-card">
                <h4>High Quality</h4>
                <div class="metric-value" style="color: #ff8c00;">{total_high}</div>
                <div class="metric-label">Promising Molecules</div>
            </div>
            <div class="metric-card">
                <h4>Good Quality</h4>
                <div class="metric-value" style="color: #4CAF50;">{total_good}</div>
                <div class="metric-label">Acceptable Molecules</div>
            </div>
        </div>
        
        <h3>Runs Overview</h3>
        <table>
            <tr>
                <th>Run Name</th>
                <th>Total Molecules</th>
                <th>Gold</th>
                <th>High</th>
                <th>Good</th>
                <th>Success Rate</th>
            </tr>
"""
        
        for a in self.analyzers:
            n_total = a.metrics.get('n_total', 1)
            n_gold = a.metrics.get('n_gold', 0)
            success_rate = (n_gold / n_total * 100) if n_total > 0 else 0
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{n_total:,}</td>
                <td><span class="badge badge-gold">{n_gold}</span></td>
                <td><span class="badge badge-high">{a.metrics.get('n_high', 0)}</span></td>
                <td><span class="badge badge-good">{a.metrics.get('n_good', 0)}</span></td>
                <td>{success_rate:.2f}%</td>
            </tr>
"""
        
        html += """
        </table>
"""
        return html
    
    def _section_gold_comparison(self):
        """金标准分子对比"""
        html = f"""
        <h2 id="gold">2. Gold Standard Molecules Comparison</h2>
        
        <div class="success">
            <strong>Gold Standard Criteria:</strong><br>
            ✓ Total Score ≥ 0.80<br>
            ✓ QSAR predicted pIC50 ≥ 6.0<br>
            ✓ QED ≥ 0.5 (drug-likeness)<br>
            ✓ Molecular Weight: 250-600 Da<br>
            ✓ Pass structural alerts<br>
            ✓ Synthetic Accessibility (SA) ≤ 5.5
        </div>
        
        <h3>Candidate Distribution by Quality</h3>
        <table>
            <tr>
                <th>Run Name</th>
                <th>Gold Count</th>
                <th>High Count</th>
                <th>Good Count</th>
                <th>Total</th>
                <th>Gold %</th>
            </tr>
"""
        
        for a in self.analyzers:
            n_gold = a.metrics.get('n_gold', 0)
            n_high = a.metrics.get('n_high', 0)
            n_good = a.metrics.get('n_good', 0)
            n_total = a.metrics.get('n_total', 1)
            gold_pct = (n_gold / n_total * 100) if n_total > 0 else 0
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td class="highlight">{n_gold}</td>
                <td>{n_high}</td>
                <td>{n_good}</td>
                <td>{n_total}</td>
                <td>{gold_pct:.2f}%</td>
            </tr>
"""
        
        html += """
        </table>
"""
        return html
    
    def _section_config_comparison(self):
        """配置对比"""
        html = f"""
        <h2 id="config">3. Configuration Analysis</h2>
        
        <h3>Model and Training Configuration</h3>
        <table>
            <tr>
                <th>Run</th>
                <th>Prior Model</th>
                <th>Agent Model</th>
                <th>Batch Size</th>
                <th>Learning Rate</th>
                <th>Sigma</th>
                <th>Scoring Type</th>
            </tr>
"""
        
        # 收集所有值用于高亮差异
        all_priors = [a.config.get('prior_file') for a in self.analyzers]
        all_batch = [a.config.get('batch_size') for a in self.analyzers]
        all_lr = [a.config.get('learning_rate') for a in self.analyzers]
        all_sigma = [a.config.get('sigma') for a in self.analyzers]
        
        for a in self.analyzers:
            prior = a.config.get('prior_file', 'N/A')
            agent = a.config.get('agent_file', 'N/A')
            batch = a.config.get('batch_size', 'N/A')
            lr = a.config.get('learning_rate', 'N/A')
            sigma = a.config.get('sigma', 'N/A')
            scoring = a.config.get('scoring_type', 'N/A')
            
            # 简化模型路径显示
            if prior != 'N/A':
                prior = os.path.basename(prior)
            if agent != 'N/A':
                agent = os.path.basename(agent)
            
            # 高亮不同的值
            batch_class = 'config-diff' if len(set(all_batch)) > 1 else ''
            lr_class = 'config-diff' if len(set(all_lr)) > 1 else ''
            sigma_class = 'config-diff' if len(set(all_sigma)) > 1 else ''
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{prior}</td>
                <td>{agent}</td>
                <td class="{batch_class}">{batch}</td>
                <td class="{lr_class}">{lr if lr == 'N/A' else f'{lr:.4f}'}</td>
                <td class="{sigma_class}">{sigma}</td>
                <td>{scoring}</td>
            </tr>
"""
        
        html += """
        </table>
        
        <div class="warning">
            <strong>⚠️ Note:</strong> Highlighted cells (pink background) indicate parameters that differ across runs.
        </div>
"""
        return html
    
    def _section_metrics_comparison(self):
        """性能指标对比"""
        html = f"""
        <h2 id="metrics">4. Performance Metrics Comparison</h2>
        
        <h3>QSAR Activity Prediction (pIC50)</h3>
        <table>
            <tr>
                <th>Run</th>
                <th>Mean ± SD</th>
                <th>Min</th>
                <th>Max</th>
                <th>Median</th>
            </tr>
"""
        
        for a in self.analyzers:
            qsar_mean = a.metrics.get('qsar_mean', 0)
            qsar_std = a.metrics.get('qsar_std', 0)
            qsar_min = a.metrics.get('qsar_min', 0)
            qsar_max = a.metrics.get('qsar_max', 0)
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{qsar_mean:.2f} ± {qsar_std:.2f}</td>
                <td>{qsar_min:.2f}</td>
                <td class="highlight">{qsar_max:.2f}</td>
                <td>{qsar_mean:.2f}</td>
            </tr>
"""
        
        html += """
        </table>
        
        <h3>Drug-likeness Properties</h3>
        <table>
            <tr>
                <th>Run</th>
                <th>QED</th>
                <th>MW (Da)</th>
                <th>LogP</th>
                <th>SA Score</th>
                <th>Total Score</th>
            </tr>
"""
        
        for a in self.analyzers:
            qed = a.metrics.get('qed_mean', 0)
            mw = a.metrics.get('mw_mean', 0)
            logp = a.metrics.get('logp_mean', 0)
            sa = a.metrics.get('sa_mean', 0)
            score = a.metrics.get('total_score_mean', 0)
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{qed:.3f}</td>
                <td>{mw:.1f}</td>
                <td>{logp:.2f}</td>
                <td>{sa:.2f}</td>
                <td class="highlight">{score:.3f}</td>
            </tr>
"""
        
        html += """
        </table>
"""
        return html
    
    def _section_component_weights(self):
        """评分组件权重对比"""
        html = f"""
        <h2 id="components">5. Scoring Component Weights Comparison</h2>
        
        <p>This section compares the scoring components and their weights across different runs.</p>
"""
        
        # 收集所有组件
        all_components = set()
        for a in self.analyzers:
            all_components.update(a.config.get('components', []))
        
        all_components = sorted(list(all_components))
        
        if all_components:
            html += """
        <h3>Component Usage Matrix</h3>
        <table>
            <tr>
                <th>Component</th>
"""
            
            for a in self.analyzers:
                html += f"<th>{a.run_name}</th>"
            
            html += "</tr>"
            
            for comp in all_components:
                html += f"<tr><td><strong>{comp}</strong></td>"
                
                for a in self.analyzers:
                    comp_weights = a.config.get('component_weights', {})
                    if comp in comp_weights:
                        weight = comp_weights[comp]
                        html += f'<td class="highlight">{weight:.2f}</td>'
                    else:
                        html += '<td style="color: #ccc;">—</td>'
                
                html += "</tr>"
            
            html += """
        </table>
        
        <h3>Total Components Count</h3>
        <table>
            <tr>
                <th>Run</th>
                <th>Number of Components</th>
                <th>Total Weight Sum</th>
            </tr>
"""
            
            for a in self.analyzers:
                n_comp = len(a.config.get('components', []))
                total_weight = sum(a.config.get('component_weights', {}).values())
                
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{n_comp}</td>
                <td>{total_weight:.2f}</td>
            </tr>
"""
            
            html += """
        </table>
"""
        else:
            html += '<p class="warning">No component information available.</p>'
        
        return html
    
    def _section_detailed_runs(self):
        """详细run信息"""
        html = f"""
        <h2 id="details">6. Detailed Run Information</h2>
"""
        
        for a in self.analyzers:
            # 获取top分子
            top_molecules = []
            if a.data.get('gold') is not None and len(a.data['gold']) > 0:
                df = a.data['gold'].head(10)
                for idx, row in df.iterrows():
                    smiles = row.get('SMILES', 'N/A')
                    score = row.get('total_score', row.get('Score', 'N/A'))
                    qsar = row.get('DENV_Activity', row.get('QSAR_Score', 'N/A'))
                    
                    top_molecules.append({
                        'rank': idx + 1,
                        'smiles': smiles,
                        'score': score,
                        'qsar': qsar
                    })
            
            html += f"""
        <h3>{a.run_name}</h3>
        
        <div class="metric-grid">
            <div class="metric-card">
                <h4>Total Molecules</h4>
                <div class="metric-value">{a.metrics.get('n_total', 0)}</div>
            </div>
            <div class="metric-card">
                <h4>Gold Candidates</h4>
                <div class="metric-value" style="color: #ffd700;">{a.metrics.get('n_gold', 0)}</div>
            </div>
            <div class="metric-card">
                <h4>Mean QSAR pIC50</h4>
                <div class="metric-value" style="color: #9b59b6;">{a.metrics.get('qsar_mean', 0):.2f}</div>
            </div>
            <div class="metric-card">
                <h4>Mean Total Score</h4>
                <div class="metric-value" style="color: #e74c3c;">{a.metrics.get('total_score_mean', 0):.3f}</div>
            </div>
        </div>
        
        <h4>Configuration Details</h4>
        <table>
            <tr>
                <th>Parameter</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Run Type</td>
                <td>{a.config.get('run_type', 'N/A')}</td>
            </tr>
            <tr>
                <td>Device</td>
                <td>{a.config.get('device', 'N/A')}</td>
            </tr>
            <tr>
                <td>Prior Model</td>
                <td>{os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'}</td>
            </tr>
            <tr>
                <td>Batch Size</td>
                <td>{a.config.get('batch_size', 'N/A')}</td>
            </tr>
            <tr>
                <td>Learning Rate</td>
                <td>{a.config.get('learning_rate', 'N/A')}</td>
            </tr>
            <tr>
                <td>Sigma</td>
                <td>{a.config.get('sigma', 'N/A')}</td>
            </tr>
            <tr>
                <td>Scoring Type</td>
                <td>{a.config.get('scoring_type', 'N/A')}</td>
            </tr>
        </table>
        
        <h4>Scoring Components</h4>
        <div class="component-list">
"""
            
            for comp in a.config.get('components', []):
                weight = a.config.get('component_weights', {}).get(comp, 1.0)
                html += f'<span class="component-tag">{comp} (w={weight:.2f})</span>'
            
            html += """
        </div>
"""
            
            if top_molecules:
                html += """
        <h4>Top 10 Gold Candidates</h4>
        <table>
            <tr>
                <th>Rank</th>
                <th>SMILES</th>
                <th>Total Score</th>
                <th>QSAR pIC50</th>
            </tr>
"""
                
                for mol in top_molecules:
                    smiles_short = mol['smiles'][:60] + '...' if len(mol['smiles']) > 60 else mol['smiles']
                    score_str = f"{mol['score']:.4f}" if isinstance(mol['score'], (int, float)) else mol['score']
                    qsar_str = f"{mol['qsar']:.2f}" if isinstance(mol['qsar'], (int, float)) else mol['qsar']
                    
                    html += f"""
            <tr>
                <td><strong>#{mol['rank']}</strong></td>
                <td><code>{smiles_short}</code></td>
                <td>{score_str}</td>
                <td>{qsar_str}</td>
            </tr>
"""
                
                html += """
        </table>
"""
            else:
                html += '<p class="warning">No gold candidates found for this run.</p>'
            
            html += '<hr style="margin: 40px 0; border: none; border-top: 2px solid #ecf0f1;">'
        
        return html
    
    def _html_footer(self):
        """HTML页脚"""
        return """
        <div style="margin-top: 60px; padding-top: 30px; border-top: 3px solid #3498db; text-align: center; color: #7f8c8d;">
            <p><strong>REINVENT4 Comprehensive Analysis Report</strong></p>
            <p>Generated by Automated Analysis Pipeline</p>
            <p style="font-size: 0.9em; margin-top: 10px;">
                For questions or issues, please contact the computational chemistry team.
            </p>
        </div>
    </div>
</body>
</html>
"""


def main():
    """主函数"""
    print("=" * 80)
    print("REINVENT4 IMPROVED RUNS ANALYSIS")
    print("=" * 80)
    
    # 查找所有runs（排除archive）
    print("\n📁 Scanning for run directories...")
    runs_dir = "experiments/runs"
    
    if not os.path.exists(runs_dir):
        print(f"❌ Error: {runs_dir} not found!")
        return
    
    run_dirs = {}
    for item in os.listdir(runs_dir):
        # 排除archive文件夹
        if item.lower() == 'archive':
            print(f"   ⊗ Skipping: {item} (excluded)")
            continue
        
        item_path = os.path.join(runs_dir, item)
        if os.path.isdir(item_path):
            run_dirs[item] = item_path
    
    print(f"   Found {len(run_dirs)} run directories")
    for run_name in sorted(run_dirs.keys()):
        print(f"     • {run_name}")
    
    if len(run_dirs) == 0:
        print("❌ No run directories found!")
        return
    
    # 分析所有runs
    print("\n📊 Analyzing runs...")
    analyzers = []
    
    for run_name, run_path in sorted(run_dirs.items()):
        analyzer = RunAnalyzer(run_name, run_path)
        try:
            analyzer.load_data()
            analyzers.append(analyzer)
        except Exception as e:
            print(f"  ⚠️  Failed to load {run_name}: {e}")
    
    if len(analyzers) == 0:
        print("❌ No runs could be loaded!")
        return
    
    print(f"\n✅ Successfully analyzed {len(analyzers)} runs")
    
    # 生成HTML报告
    print("\n📄 Generating reports...")
    reporter = HTMLReportGenerator(analyzers)
    html_file = reporter.generate(output_file='runs_analysis_report.html')
    
    # 导出CSV数据
    print("\n📊 Exporting data tables...")
    
    # 1. 摘要表
    summary_data = []
    for a in analyzers:
        summary_data.append({
            'Run_Name': a.run_name,
            'Total_Molecules': a.metrics.get('n_total', 0),
            'Gold_Count': a.metrics.get('n_gold', 0),
            'High_Count': a.metrics.get('n_high', 0),
            'Good_Count': a.metrics.get('n_good', 0),
            'Mean_QSAR_pIC50': f"{a.metrics.get('qsar_mean', 0):.2f}",
            'Max_QSAR_pIC50': f"{a.metrics.get('qsar_max', 0):.2f}",
            'Mean_QED': f"{a.metrics.get('qed_mean', 0):.3f}",
            'Mean_MW': f"{a.metrics.get('mw_mean', 0):.1f}",
            'Mean_Total_Score': f"{a.metrics.get('total_score_mean', 0):.3f}",
            'Prior_Model': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
            'Batch_Size': a.config.get('batch_size', 'N/A'),
            'Learning_Rate': a.config.get('learning_rate', 'N/A'),
            'Sigma': a.config.get('sigma', 'N/A'),
            'Num_Components': len(a.config.get('components', [])),
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_csv = 'runs_summary.csv'
    summary_df.to_csv(summary_csv, index=False)
    print(f"   ✓ {summary_csv}")
    
    # 2. 组件权重对比表
    all_components = set()
    for a in analyzers:
        all_components.update(a.config.get('components', []))
    
    if all_components:
        comp_data = []
        for comp in sorted(all_components):
            row = {'Component': comp}
            for a in analyzers:
                weight = a.config.get('component_weights', {}).get(comp, None)
                row[a.run_name] = weight if weight is not None else ''
            comp_data.append(row)
        
        comp_df = pd.DataFrame(comp_data)
        comp_csv = 'component_weights_comparison.csv'
        comp_df.to_csv(comp_csv, index=False)
        print(f"   ✓ {comp_csv}")
    
    # 3. 配置对比表
    config_data = []
    for a in analyzers:
        config_data.append({
            'Run': a.run_name,
            'Run_Type': a.config.get('run_type', 'N/A'),
            'Device': a.config.get('device', 'N/A'),
            'Prior_Model': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
            'Agent_Model': os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A',
            'Batch_Size': a.config.get('batch_size', 'N/A'),
            'Learning_Rate': a.config.get('learning_rate', 'N/A'),
            'Sigma': a.config.get('sigma', 'N/A'),
            'Scoring_Type': a.config.get('scoring_type', 'N/A'),
            'Num_Components': len(a.config.get('components', [])),
        })
    
    config_df = pd.DataFrame(config_data)
    config_csv = 'configuration_comparison.csv'
    config_df.to_csv(config_csv, index=False)
    print(f"   ✓ {config_csv}")
    
    # 显示快速摘要
    print("\n" + "=" * 80)
    print("QUICK SUMMARY")
    print("=" * 80)
    print("\n" + summary_df.to_string(index=False))
    
    print("\n" + "=" * 80)
    print("✅ ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"\n📁 Generated files:")
    print(f"  1. {html_file} - Interactive HTML report (open in browser)")
    print(f"  2. {summary_csv} - Summary statistics table")
    print(f"  3. {comp_csv} - Component weights comparison")
    print(f"  4. {config_csv} - Configuration comparison")
    print(f"\n🌐 Open {html_file} in your browser to view the full report!")
    print("=" * 80)


if __name__ == "__main__":
    main()
