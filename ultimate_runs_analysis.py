#!/usr/bin/env python3
"""
REINVENT4 终极分析工具
- 正确的金标准定义 (pIC50>8, QED>0.7, SA<4, MW 300-500, LogP 1-4)
- 毒性组件分析
- Agent模型对比
- 核心骨架分析
- 金标准分子深度比较
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
from collections import Counter
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
            'smiles_file': None,
            'batch_size': None,
            'learning_rate': None,
            'sigma': None,
            'scoring_type': None,
            'components': [],
            'component_weights': {},
            'has_toxicity': False,
        }
        
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            
            config['run_type'] = ConfigParser._extract(content, r'run_type\s*=\s*"([^"]+)"')
            config['device'] = ConfigParser._extract(content, r'device\s*=\s*"([^"]+)"')
            config['prior_file'] = ConfigParser._extract(content, r'prior_file\s*=\s*"([^"]+)"')
            config['agent_file'] = ConfigParser._extract(content, r'agent_file\s*=\s*"([^"]+)"')
            config['smiles_file'] = ConfigParser._extract(content, r'smiles_file\s*=\s*"([^"]+)"')
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
            
            # 检测毒性相关组件
            toxicity_keywords = ['toxicity', 'toxic', 'alert', 'admet', 'safety']
            for comp in config['components']:
                if any(kw in comp.lower() for kw in toxicity_keywords):
                    config['has_toxicity'] = True
                    break
            
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
        self.gold_molecules = []
        self.scaffolds = set()
    
    def load_data(self):
        """加载数据"""
        print(f"  分析 {self.run_name}...")
        
        # 加载配置
        config_path = os.path.join(self.run_dir, 'config.toml')
        if os.path.exists(config_path):
            self.config = ConfigParser.parse_toml(config_path)
            print(f"    ✓ 配置: {len(self.config['components'])} 个组件" + 
                  (" [含毒性]" if self.config.get('has_toxicity') else ""))
        
        # 智能查找CSV文件
        csv_mapping = {
            'results': ['results_*.csv', 'results.csv'],
            'gold': ['candidates_gold.csv', 'candidates_金标准*.csv', 'promising*.csv'],
            'high': ['candidates_high.csv', 'candidates_高标准*.csv'],
            'good': ['candidates_good.csv', 'candidates_中标准*.csv'],
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
        self._apply_gold_standard()
        self._extract_scaffolds()
    
    def _calculate_metrics(self):
        """计算统计指标"""
        m = {}
        
        m['n_total'] = len(self.data.get('results', []))
        m['n_gold_original'] = len(self.data.get('gold', []))
        m['n_high'] = len(self.data.get('high', []))
        m['n_good'] = len(self.data.get('good', []))
        
        df = self.data.get('gold')
        if df is None or len(df) == 0:
            df = self.data.get('results')
        
        if df is not None and len(df) > 0:
            df.columns = df.columns.str.strip()
            
            col_map = {
                'total_score': ['Score', 'total_score', 'score'],
                'qsar': ['DENV_Activity (raw)', 'DENV_Activity', 'QSAR_Score', 'qsar', 'pIC50'],
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
    
    def _apply_gold_standard(self):
        """应用新的金标准筛选"""
        df = self.data.get('gold')
        if df is None or len(df) == 0:
            df = self.data.get('results')
        
        if df is None or len(df) == 0:
            self.metrics['n_gold_strict'] = 0
            return
        
        df = df.copy()
        df.columns = df.columns.str.strip()
        
        # 提取各个指标
        def get_col_value(df, possible_names):
            for name in possible_names:
                if name in df.columns:
                    return pd.to_numeric(df[name], errors='coerce')
            return pd.Series([np.nan] * len(df))
        
        pic50 = get_col_value(df, ['DENV_Activity (raw)', 'DENV_Activity', 'pIC50', 'QSAR_Score'])
        qed = get_col_value(df, ['QED (raw)', 'QED', 'qed'])
        sa = get_col_value(df, ['SA (raw)', 'SA', 'SAScore', 'sa_score'])
        mw = get_col_value(df, ['MW', 'Molecular weight', 'molecular_weight'])
        logp = get_col_value(df, ['LogP', 'SlogP (RDKit)', 'SlogP'])
        smiles = df.get('SMILES', pd.Series(['']*len(df)))
        
        # 应用金标准: pIC50>8.0, QED>0.7, SA<4.0, MW 300-500, LogP 1-4
        mask = (
            (pic50 > 8.0) &
            (qed > 0.7) &
            (sa < 4.0) &
            (mw >= 300) & (mw <= 500) &
            (logp >= 1) & (logp <= 4)
        )
        
        gold_strict = df[mask]
        self.metrics['n_gold_strict'] = len(gold_strict)
        
        # 保存金标准分子用于后续分析
        if len(gold_strict) > 0:
            for idx, row in gold_strict.iterrows():
                self.gold_molecules.append({
                    'run': self.run_name,
                    'smiles': row.get('SMILES', ''),
                    'pic50': pic50.loc[idx] if idx in pic50.index else np.nan,
                    'qed': qed.loc[idx] if idx in qed.index else np.nan,
                    'sa': sa.loc[idx] if idx in sa.index else np.nan,
                    'mw': mw.loc[idx] if idx in mw.index else np.nan,
                    'logp': logp.loc[idx] if idx in logp.index else np.nan,
                })
    
    def _extract_scaffolds(self):
        """提取核心骨架"""
        try:
            from rdkit import Chem
            from rdkit.Chem.Scaffolds import MurckoScaffold
            
            for mol_data in self.gold_molecules:
                smiles = mol_data['smiles']
                if smiles:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        try:
                            scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                            scaffold_smiles = Chem.MolToSmiles(scaffold)
                            self.scaffolds.add(scaffold_smiles)
                        except:
                            pass
        except ImportError:
            print("    ⚠️  RDKit未安装，跳过骨架分析")


class HTMLReportGenerator:
    """生成HTML报告"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    def generate(self, output_file='ultimate_runs_report.html'):
        """生成完整报告"""
        print("\n📄 生成终极HTML报告...")
        
        html = self._html_header()
        html += self._section_gold_standard_definition()
        html += self._section_summary()
        html += self._section_gold_strict_comparison()
        html += self._section_structures()
        html += self._section_config_comparison()
        html += self._section_toxicity_analysis()
        html += self._section_scaffold_analysis()
        html += self._section_metrics_detailed()
        html += self._section_components()
        html += self._section_gold_molecules_comparison()
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
    <title>REINVENT4 Ultimate Analysis</title>
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
            border-bottom: 4px solid #3498db; 
            padding-bottom: 10px;
            font-size: 1.8em;
        }}
        h2 {{ 
            color: #34495e; 
            margin-top: 30px; 
            border-left: 5px solid #3498db; 
            padding-left: 15px;
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
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; 
            padding: 12px; 
            text-align: left;
            position: sticky;
            top: 0;
            font-size: 0.9em;
        }}
        td {{ 
            padding: 10px; 
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
            border-radius: 10px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
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
        
        .gold-standard-box {{
            background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 4px 10px rgba(0,0,0,0.15);
        }}
        .gold-standard-box h3 {{
            color: #000;
            margin-top: 0;
        }}
        .criterion {{
            background: white;
            padding: 10px;
            margin: 5px 0;
            border-radius: 5px;
            font-weight: 600;
            color: #2c3e50;
        }}
        
        .toxicity-badge {{
            background: #e74c3c;
            color: white;
            padding: 4px 10px;
            border-radius: 15px;
            font-size: 0.85em;
            font-weight: 600;
        }}
        
        .scaffold-box {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            border-left: 4px solid #3498db;
        }}
        
        .comparison-box {{
            background: #e8f4f8;
            padding: 20px;
            border-radius: 10px;
            margin: 20px 0;
        }}
        
        .highlight {{ background: #fffacd; font-weight: bold; }}
        .best {{ background: #90EE90; font-weight: bold; }}
        
        @media (max-width: 768px) {{
            .container {{ padding: 15px; }}
            h1 {{ font-size: 1.5em; }}
            h2 {{ font-size: 1.2em; }}
            table {{ font-size: 0.75em; }}
            .metric-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 REINVENT4 Ultimate Multi-Run Analysis</h1>
    <p style="text-align: center; color: #7f8c8d; margin-bottom: 20px;">
        Generated: {self.timestamp} | With Gold Standard Filtering & Deep Comparison
    </p>
"""
    
    def _section_gold_standard_definition(self):
        return """
    <div class="gold-standard-box">
        <h3>⭐ 金标准分子定义 (Strict Gold Standard)</h3>
        <p style="color: #000; margin-bottom: 15px;">
            金标准分子必须<strong>同时满足</strong>以下所有条件：
        </p>
        <div class="criterion">✓ 预测活性: pIC50 > 8.0 (IC50 < 10 nM)</div>
        <div class="criterion">✓ 类药性: QED > 0.7</div>
        <div class="criterion">✓ 合成可及性: SA Score < 4.0 (易于合成)</div>
        <div class="criterion">✓ 分子量: 300-500 Da</div>
        <div class="criterion">✓ 亲脂性: LogP 1-4 (适中)</div>
        <p style="color: #000; margin-top: 15px; font-size: 0.9em;">
            <strong>注意:</strong> 此标准比原始筛选更严格，旨在获得最优质的候选分子
        </p>
    </div>
"""
    
    def _section_summary(self):
        total_gold = sum(a.metrics.get('n_gold_strict', 0) for a in self.analyzers)
        total_original = sum(a.metrics.get('n_gold_original', 0) for a in self.analyzers)
        
        html = f"""
    <h2>1. 📊 执行摘要</h2>
    <div class="metric-grid">
        <div class="metric-card">
            <div class="metric-label">总运行数</div>
            <div class="metric-value">{len(self.analyzers)}</div>
        </div>
        <div class="metric-card" style="background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);">
            <div class="metric-label">严格金标准分子</div>
            <div class="metric-value">{total_gold}</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">原始Gold分子</div>
            <div class="metric-value">{total_original}</div>
        </div>
        <div class="metric-card">
            <div class="metric-label">筛选率</div>
            <div class="metric-value">{(total_gold/total_original*100) if total_original > 0 else 0:.1f}%</div>
        </div>
    </div>
"""
        return html
    
    def _section_gold_strict_comparison(self):
        html = """
    <h2>2. 🥇 严格金标准筛选对比</h2>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>原始Gold</th>
                <th>严格Gold</th>
                <th>通过率</th>
                <th>平均pIC50</th>
                <th>平均QED</th>
                <th>平均SA</th>
            </tr>
        </thead>
        <tbody>
"""
        
        best_pass_rate = 0
        best_run = None
        
        for a in self.analyzers:
            n_orig = a.metrics.get('n_gold_original', 0)
            n_strict = a.metrics.get('n_gold_strict', 0)
            pass_rate = (n_strict/n_orig*100) if n_orig > 0 else 0
            
            if pass_rate > best_pass_rate:
                best_pass_rate = pass_rate
                best_run = a.run_name
            
            # 计算gold分子的平均指标
            gold_mols = [m for m in a.gold_molecules]
            avg_pic50 = np.mean([m['pic50'] for m in gold_mols]) if gold_mols else 0
            avg_qed = np.mean([m['qed'] for m in gold_mols]) if gold_mols else 0
            avg_sa = np.mean([m['sa'] for m in gold_mols]) if gold_mols else 0
            
            row_class = 'best' if a.run_name == best_run and n_strict > 0 else ''
            
            html += f"""
            <tr class="{row_class}">
                <td><strong>{a.run_name}</strong></td>
                <td>{n_orig}</td>
                <td class="highlight">{n_strict}</td>
                <td>{pass_rate:.1f}%</td>
                <td>{avg_pic50:.2f}</td>
                <td>{avg_qed:.2f}</td>
                <td>{avg_sa:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    <p style="color: #27ae60; font-weight: 600;">
        ✓ 绿色高亮表示通过率最高的run
    </p>
"""
        return html
    
    def _section_structures(self):
        """分子结构对比 - Top严格金标准分子"""
        html = """
    <h2>3. 🔬 Top严格金标准分子结构</h2>
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px;">
"""
        
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                html += f"""
        <div style="border: 3px solid #FFD700; border-radius: 10px; padding: 15px; background: #fffef0;">
            <h4 style="color: #000; margin: 0 0 15px 0;">{a.run_name} (Top 3)</h4>
"""
                
                for i, mol in enumerate(a.gold_molecules[:3]):
                    img_html = self._generate_molecule_image(mol['smiles'], a.run_name, i)
                    
                    html += f"""
            <div style="margin: 15px 0; padding: 15px; background: white; border-radius: 8px; border: 1px solid #ddd;">
                <strong>#{i+1}</strong><br>
                {img_html}
                <div style="margin-top: 10px; font-size: 0.85em;">
                    <strong>pIC50:</strong> {mol['pic50']:.2f}<br>
                    <strong>QED:</strong> {mol['qed']:.2f}<br>
                    <strong>SA:</strong> {mol['sa']:.2f}<br>
                    <strong>MW:</strong> {mol['mw']:.1f}<br>
                    <strong>LogP:</strong> {mol['logp']:.2f}
                </div>
                <div style="margin-top: 5px; font-size: 0.75em; word-break: break-all; color: #666;">
                    {mol['smiles'][:60]}{'...' if len(mol['smiles']) > 60 else ''}
                </div>
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
        """生成分子图片"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(280, 280))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                return f'<img src="data:image/png;base64,{img_str}" style="width: 100%; max-width: 280px; height: auto; border-radius: 4px;">'
        except:
            pass
        
        return '<div style="background: #f0f0f0; padding: 40px; text-align: center; color: #999; border-radius: 4px;">需要RDKit显示结构</div>'
    
    def _section_config_comparison(self):
        html = """
    <h2>4. ⚙️ 详细配置对比</h2>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>Prior模型</th>
                <th>Agent模型</th>
                <th>骨架文件</th>
                <th>批次</th>
                <th>学习率</th>
                <th>Sigma</th>
                <th>组件数</th>
                <th>毒性组件</th>
            </tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            prior = os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'
            agent = os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A'
            smiles_file = os.path.basename(a.config.get('smiles_file', 'N/A')) if a.config.get('smiles_file') else 'N/A'
            lr = a.config.get('learning_rate')
            lr_str = f"{lr:.4f}" if lr is not None else 'N/A'
            has_tox = a.config.get('has_toxicity', False)
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{prior}</td>
                <td><strong>{agent}</strong></td>
                <td>{smiles_file}</td>
                <td>{a.config.get('batch_size', 'N/A')}</td>
                <td>{lr_str}</td>
                <td>{a.config.get('sigma', 'N/A')}</td>
                <td>{len(a.config.get('components', []))}</td>
                <td>{'<span class="toxicity-badge">含毒性</span>' if has_tox else '—'}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_toxicity_analysis(self):
        """毒性组件分析"""
        runs_with_tox = [a for a in self.analyzers if a.config.get('has_toxicity')]
        runs_without_tox = [a for a in self.analyzers if not a.config.get('has_toxicity')]
        
        html = """
    <h2>5. ☠️ 毒性组件分析</h2>
"""
        
        if runs_with_tox:
            avg_with = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_with_tox])
            html += f"""
    <div class="comparison-box">
        <h3>含毒性组件的Runs ({len(runs_with_tox)}个)</h3>
        <p><strong>平均严格金标准数:</strong> {avg_with:.1f}</p>
        <ul>
"""
            for a in runs_with_tox:
                html += f"<li><strong>{a.run_name}</strong>: {a.metrics.get('n_gold_strict', 0)} 个严格金标准</li>"
            
            html += """
        </ul>
    </div>
"""
        
        if runs_without_tox:
            avg_without = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_without_tox])
            html += f"""
    <div class="comparison-box">
        <h3>不含毒性组件的Runs ({len(runs_without_tox)}个)</h3>
        <p><strong>平均严格金标准数:</strong> {avg_without:.1f}</p>
        <ul>
"""
            for a in runs_without_tox:
                html += f"<li><strong>{a.run_name}</strong>: {a.metrics.get('n_gold_strict', 0)} 个严格金标准</li>"
            
            html += """
        </ul>
    </div>
"""
        
        if runs_with_tox and runs_without_tox:
            improvement = ((avg_with - avg_without) / avg_without * 100) if avg_without > 0 else 0
            comparison_text = f"含毒性组件的runs平均产生 <strong>{improvement:+.1f}%</strong> 的严格金标准分子"
            
            html += f"""
    <div style="background: #fff3cd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffc107;">
        <strong>💡 洞察:</strong> {comparison_text}
    </div>
"""
        
        return html
    
    def _section_scaffold_analysis(self):
        """核心骨架分析"""
        html = """
    <h2>6. 🧩 核心骨架分析</h2>
"""
        
        # 统计每个run的骨架数量
        scaffold_data = []
        all_scaffolds = set()
        
        for a in self.analyzers:
            n_scaffolds = len(a.scaffolds)
            scaffold_data.append({
                'run': a.run_name,
                'n_scaffolds': n_scaffolds,
                'scaffolds': a.scaffolds
            })
            all_scaffolds.update(a.scaffolds)
        
        html += f"""
    <div class="scaffold-box">
        <h3>骨架多样性统计</h3>
        <p><strong>所有runs共产生 {len(all_scaffolds)} 个不同的Murcko骨架</strong></p>
    </div>
    
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>独特骨架数</th>
                <th>严格金标准数</th>
                <th>骨架多样性比</th>
            </tr>
        </thead>
        <tbody>
"""
        
        for data in scaffold_data:
            n_gold = next((a.metrics.get('n_gold_strict', 0) for a in self.analyzers if a.run_name == data['run']), 0)
            diversity_ratio = (data['n_scaffolds'] / n_gold) if n_gold > 0 else 0
            
            html += f"""
            <tr>
                <td><strong>{data['run']}</strong></td>
                <td class="highlight">{data['n_scaffolds']}</td>
                <td>{n_gold}</td>
                <td>{diversity_ratio:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    
    <h3>骨架共享分析</h3>
"""
        
        # 分析runs之间的骨架重叠
        if len(self.analyzers) >= 2:
            html += '<div style="margin: 20px 0;">'
            
            for i, a1 in enumerate(self.analyzers):
                for a2 in self.analyzers[i+1:]:
                    shared = a1.scaffolds & a2.scaffolds
                    unique_a1 = a1.scaffolds - a2.scaffolds
                    unique_a2 = a2.scaffolds - a1.scaffolds
                    
                    html += f"""
    <div class="scaffold-box">
        <strong>{a1.run_name} vs {a2.run_name}:</strong><br>
        • 共享骨架: {len(shared)}<br>
        • {a1.run_name}独有: {len(unique_a1)}<br>
        • {a2.run_name}独有: {len(unique_a2)}
    </div>
"""
            
            html += '</div>'
        
        return html
    
    def _section_metrics_detailed(self):
        """详细性能指标"""
        html = """
    <h2>7. 📈 详细性能指标 (仅严格金标准分子)</h2>
    
    <h3>活性预测 (pIC50)</h3>
    <table>
        <thead>
            <tr><th>Run</th><th>均值 ± 标准差</th><th>最小值</th><th>最大值</th><th>范围</th></tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                pic50_vals = [m['pic50'] for m in a.gold_molecules if not np.isnan(m['pic50'])]
                if pic50_vals:
                    mean_val = np.mean(pic50_vals)
                    std_val = np.std(pic50_vals)
                    min_val = np.min(pic50_vals)
                    max_val = np.max(pic50_vals)
                    range_val = max_val - min_val
                    
                    html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{mean_val:.2f} ± {std_val:.2f}</td>
                <td>{min_val:.2f}</td>
                <td class="highlight">{max_val:.2f}</td>
                <td>{range_val:.2f}</td>
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
            if len(a.gold_molecules) > 0:
                qed_vals = [m['qed'] for m in a.gold_molecules if not np.isnan(m['qed'])]
                mw_vals = [m['mw'] for m in a.gold_molecules if not np.isnan(m['mw'])]
                logp_vals = [m['logp'] for m in a.gold_molecules if not np.isnan(m['logp'])]
                sa_vals = [m['sa'] for m in a.gold_molecules if not np.isnan(m['sa'])]
                
                qed_str = f"{np.mean(qed_vals):.2f} ± {np.std(qed_vals):.2f}" if qed_vals else 'N/A'
                mw_str = f"{np.mean(mw_vals):.1f} ± {np.std(mw_vals):.1f}" if mw_vals else 'N/A'
                logp_str = f"{np.mean(logp_vals):.2f} ± {np.std(logp_vals):.2f}" if logp_vals else 'N/A'
                sa_str = f"{np.mean(sa_vals):.2f} ± {np.std(sa_vals):.2f}" if sa_vals else 'N/A'
                
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{qed_str}</td>
                <td>{mw_str}</td>
                <td>{logp_str}</td>
                <td>{sa_str}</td>
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
    <h2>8. 🎯 评分组件权重详细对比</h2>
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
    
    <h3>组件总权重对比</h3>
    <table>
        <thead>
            <tr><th>Run</th><th>组件总数</th><th>权重总和</th><th>平均权重</th></tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            n_comp = len(a.config.get('components', []))
            total_weight = sum(a.config.get('component_weights', {}).values())
            avg_weight = total_weight / n_comp if n_comp > 0 else 0
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{n_comp}</td>
                <td>{total_weight:.2f}</td>
                <td>{avg_weight:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_gold_molecules_comparison(self):
        """金标准分子深度比较"""
        html = """
    <h2>9. 💎 严格金标准分子深度比较</h2>
    
    <h3>最佳分子对比</h3>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>SMILES</th>
                <th>pIC50</th>
                <th>QED</th>
                <th>SA</th>
                <th>MW</th>
                <th>LogP</th>
            </tr>
        </thead>
        <tbody>
"""
        
        # 找出每个run的最佳分子（按pIC50）
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                best_mol = max(a.gold_molecules, key=lambda x: x['pic50'] if not np.isnan(x['pic50']) else 0)
                smiles_short = best_mol['smiles'][:50] + '...' if len(best_mol['smiles']) > 50 else best_mol['smiles']
                
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td style="font-family: monospace; font-size: 0.8em;">{smiles_short}</td>
                <td class="highlight">{best_mol['pic50']:.2f}</td>
                <td>{best_mol['qed']:.2f}</td>
                <td>{best_mol['sa']:.2f}</td>
                <td>{best_mol['mw']:.1f}</td>
                <td>{best_mol['logp']:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    
    <h3>优缺点分析</h3>
"""
        
        # 分析每个run的优缺点
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                pic50_vals = [m['pic50'] for m in a.gold_molecules if not np.isnan(m['pic50'])]
                qed_vals = [m['qed'] for m in a.gold_molecules if not np.isnan(m['qed'])]
                sa_vals = [m['sa'] for m in a.gold_molecules if not np.isnan(m['sa'])]
                
                avg_pic50 = np.mean(pic50_vals) if pic50_vals else 0
                avg_qed = np.mean(qed_vals) if qed_vals else 0
                avg_sa = np.mean(sa_vals) if sa_vals else 0
                
                strengths = []
                weaknesses = []
                
                # 判断优点
                if avg_pic50 > 8.5:
                    strengths.append("非常高的活性预测")
                if avg_qed > 0.75:
                    strengths.append("优秀的类药性")
                if avg_sa < 3.5:
                    strengths.append("极佳的合成可及性")
                if len(a.scaffolds) > 5:
                    strengths.append("骨架多样性好")
                if a.config.get('has_toxicity'):
                    strengths.append("包含毒性筛选")
                
                # 判断缺点
                if avg_pic50 < 8.3:
                    weaknesses.append("活性预测相对较低")
                if avg_qed < 0.73:
                    weaknesses.append("类药性有待提升")
                if avg_sa > 3.7:
                    weaknesses.append("合成难度较高")
                if len(a.scaffolds) < 3:
                    weaknesses.append("骨架多样性不足")
                if len(a.gold_molecules) < 10:
                    weaknesses.append("严格金标准产量较低")
                
                html += f"""
    <div class="comparison-box">
        <h4>{a.run_name}</h4>
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 15px;">
            <div>
                <strong style="color: #27ae60;">✓ 优点:</strong>
                <ul style="margin: 5px 0 0 20px;">
"""
                
                for s in strengths:
                    html += f"<li>{s}</li>"
                
                if not strengths:
                    html += "<li>暂无显著优势</li>"
                
                html += """
                </ul>
            </div>
            <div>
                <strong style="color: #e74c3c;">✗ 需改进:</strong>
                <ul style="margin: 5px 0 0 20px;">
"""
                
                for w in weaknesses:
                    html += f"<li>{w}</li>"
                
                if not weaknesses:
                    html += "<li>表现全面优秀</li>"
                
                html += """
                </ul>
            </div>
        </div>
    </div>
"""
        
        # 总体推荐
        best_run = max(self.analyzers, key=lambda x: x.metrics.get('n_gold_strict', 0))
        
        html += f"""
    <div style="background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%); 
                padding: 20px; border-radius: 10px; margin-top: 30px; color: #000;">
        <h3 style="margin-top: 0;">🏆 推荐策略</h3>
        <p><strong>最佳Run:</strong> {best_run.run_name}</p>
        <p><strong>产生严格金标准数:</strong> {best_run.metrics.get('n_gold_strict', 0)}</p>
        <p><strong>推荐理由:</strong> 该run在严格金标准筛选下产生了最多的高质量候选分子，
        综合考虑活性、类药性、合成可及性等多个维度表现优异。</p>
    </div>
"""
        
        return html
    
    def _html_footer(self):
        return """
    <div style="margin-top: 60px; padding-top: 30px; border-top: 3px solid #3498db; 
                text-align: center; color: #7f8c8d;">
        <p><strong>REINVENT4 Ultimate Comprehensive Analysis Report</strong></p>
        <p style="font-size: 0.9em; margin-top: 10px;">
            With Strict Gold Standard Filtering, Toxicity Analysis, 
            Scaffold Diversity & Deep Molecular Comparison
        </p>
        <p style="font-size: 0.85em; margin-top: 15px; color: #95a5a6;">
            📱 Optimized for Mobile & Desktop | 🔬 Research Grade Analysis
        </p>
    </div>
</div>
</body>
</html>
"""


def main():
    print("="*80)
    print("REINVENT4 终极分析工具")
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
            gold_mols = a.gold_molecules
            avg_pic50 = np.mean([m['pic50'] for m in gold_mols if not np.isnan(m['pic50'])]) if gold_mols else 0
            avg_qed = np.mean([m['qed'] for m in gold_mols if not np.isnan(m['qed'])]) if gold_mols else 0
            
            summary_data.append({
                'Run': a.run_name,
                'Original_Gold': a.metrics.get('n_gold_original', 0),
                'Strict_Gold': a.metrics.get('n_gold_strict', 0),
                'Pass_Rate_%': f"{(a.metrics.get('n_gold_strict', 0)/a.metrics.get('n_gold_original', 1)*100):.1f}",
                'Avg_pIC50': f"{avg_pic50:.2f}",
                'Avg_QED': f"{avg_qed:.2f}",
                'N_Scaffolds': len(a.scaffolds),
                'Has_Toxicity': 'Yes' if a.config.get('has_toxicity') else 'No',
                'Prior': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
                'Agent': os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A',
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv('ultimate_runs_summary.csv', index=False)
        print("\n✅ CSV: ultimate_runs_summary.csv")
        print("\n" + df.to_string(index=False))
    
    print("\n" + "="*80)
    print("🏆 终极分析完成！")
    print("📊 包含:")
    print("   • 严格金标准筛选 (pIC50>8, QED>0.7, SA<4, MW 300-500, LogP 1-4)")
    print("   • 毒性组件分析")
    print("   • Agent模型对比")
    print("   • 核心骨架多样性分析")
    print("   • 深度分子比较与优缺点分析")
    print("="*80)

if __name__ == "__main__":
    main()
