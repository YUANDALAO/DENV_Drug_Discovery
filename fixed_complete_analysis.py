#!/usr/bin/env python3
"""
REINVENT4 完整终极分析工具 - 修复版
修复了格式字符串错误
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

# ... (前面的ConfigParser和RunAnalyzer类保持不变，直接从上一个脚本复制)

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
            'toxicity_components': [],
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
            
            toxicity_names = {'Toxicity_Alerts', 'PAINS_Filter', 'Metabolic_Stability'}
            for comp_detail in components:
                display_name = comp_detail.get('display_name', '')
                if display_name in toxicity_names:
                    config['has_toxicity'] = True
                    if display_name not in config['toxicity_components']:
                        config['toxicity_components'].append(display_name)
            
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
        print(f"  分析 {self.run_name}...")
        
        config_path = os.path.join(self.run_dir, 'config.toml')
        if os.path.exists(config_path):
            self.config = ConfigParser.parse_toml(config_path)
            tox_info = ""
            if self.config.get('has_toxicity'):
                tox_comps = self.config.get('toxicity_components', [])
                tox_info = f" [🧪 {', '.join(tox_comps)}]"
            print(f"    ✓ 配置: {len(self.config['components'])} 个组件{tox_info}")
        
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
                'total_score': ['Score', 'total_score'],
                'qsar': ['DENV_Activity (raw)', 'DENV_Activity'],
                'qed': ['QED (raw)', 'QED'],
                'mw': ['MW (raw)', 'MW'],
                'sa': ['SA (raw)', 'SA'],
                'logp': ['LogP (raw)', 'LogP'],
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
        df = self.data.get('gold')
        if df is None or len(df) == 0:
            self.metrics['n_gold_strict'] = 0
            return
        
        df = df.copy()
        df.columns = df.columns.str.strip()
        
        def get_col(possible_names):
            for name in possible_names:
                if name in df.columns:
                    return pd.to_numeric(df[name], errors='coerce')
            return pd.Series([np.nan] * len(df))
        
        pic50 = get_col(['DENV_Activity (raw)', 'DENV_Activity'])
        qed = get_col(['QED (raw)', 'QED'])
        sa = get_col(['SA (raw)', 'SA'])
        mw = get_col(['MW (raw)', 'MW'])
        logp = get_col(['LogP (raw)', 'LogP'])
        
        mask = (
            (pic50 >= 8.0) &
            (qed >= 0.7) &
            (sa <= 4.0) &
            (mw >= 300) & (mw <= 500) &
            (logp >= 1) & (logp <= 4)
        )
        
        gold_strict = df[mask]
        self.metrics['n_gold_strict'] = len(gold_strict)
        
        if len(gold_strict) > 0:
            for idx, row in gold_strict.iterrows():
                self.gold_molecules.append({
                    'run': self.run_name,
                    'smiles': row.get('SMILES', ''),
                    'score': row.get('Score', np.nan),
                    'pic50': pic50.loc[idx] if idx in pic50.index else np.nan,
                    'qed': qed.loc[idx] if idx in qed.index else np.nan,
                    'sa': sa.loc[idx] if idx in sa.index else np.nan,
                    'mw': mw.loc[idx] if idx in mw.index else np.nan,
                    'logp': logp.loc[idx] if idx in logp.index else np.nan,
                })
    
    def _extract_scaffolds(self):
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
            pass


class HTMLReportGenerator:
    """生成HTML报告"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    def generate(self, output_file='complete_runs_report.html'):
        print("\n📄 生成完整HTML报告...")
        
        html = self._html_header()
        html += self._section_gold_standard()
        html += self._section_summary()
        html += self._section_comparison()
        html += self._section_structures()
        html += self._section_toxicity()
        html += self._section_scaffolds()
        html += self._section_config()
        html += self._section_components()
        html += self._html_footer()
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        
        print(f"✅ HTML报告: {output_file}")
        return output_file
    
    def _html_header(self):
        return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>REINVENT4 Analysis</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: Arial, sans-serif; line-height: 1.6; 
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 10px; }}
        .container {{ max-width: 1400px; margin: 0 auto; background: white; 
                      padding: 20px; border-radius: 12px; box-shadow: 0 10px 40px rgba(0,0,0,0.2); }}
        h1 {{ color: #2c3e50; border-bottom: 4px solid #3498db; padding-bottom: 10px; font-size: 1.8em; }}
        h2 {{ color: #34495e; margin-top: 30px; border-left: 5px solid #3498db; 
              padding-left: 15px; font-size: 1.4em; }}
        h3 {{ color: #555; margin-top: 20px; font-size: 1.2em; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.9em; }}
        th {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
              color: white; padding: 12px; text-align: left; font-size: 0.85em; }}
        td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background: #f5f5f5; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); 
                        gap: 15px; margin: 20px 0; }}
        .metric-card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        color: white; padding: 20px; border-radius: 10px; text-align: center; }}
        .metric-value {{ font-size: 2.5em; font-weight: bold; margin: 10px 0; }}
        .gold-box {{ background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%); 
                     padding: 20px; border-radius: 10px; margin: 20px 0; color: #000; }}
        .criterion {{ background: white; padding: 10px; margin: 5px 0; border-radius: 5px; font-weight: 600; }}
        .toxicity-badge {{ background: #e74c3c; color: white; padding: 4px 10px; 
                           border-radius: 15px; font-size: 0.85em; }}
        .highlight {{ background: #fffacd; font-weight: bold; }}
        .best {{ background: #90EE90; font-weight: bold; }}
        @media (max-width: 768px) {{
            .container {{ padding: 15px; }}
            h1 {{ font-size: 1.5em; }}
            table {{ font-size: 0.75em; }}
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 REINVENT4 Complete Analysis</h1>
    <p style="text-align: center; color: #7f8c8d; margin-bottom: 20px;">
        Generated: {self.timestamp}
    </p>
"""
    
    def _section_gold_standard(self):
        return """
    <div class="gold-box">
        <h3>⭐ 金标准分子定义</h3>
        <p style="margin-bottom: 15px;">必须<strong>同时满足</strong>以下条件：</p>
        <div class="criterion">✓ pIC50 ≥ 8.0 (IC50 ≤ 10 nM)</div>
        <div class="criterion">✓ QED ≥ 0.7</div>
        <div class="criterion">✓ SA ≤ 4.0</div>
        <div class="criterion">✓ MW: 300-500 Da</div>
        <div class="criterion">✓ LogP: 1-4</div>
    </div>
"""
    
    def _section_summary(self):
        total_gold = sum(a.metrics.get('n_gold_strict', 0) for a in self.analyzers)
        total_orig = sum(a.metrics.get('n_gold_original', 0) for a in self.analyzers)
        
        pass_rate = (total_gold/total_orig*100) if total_orig > 0 else 0
        
        html = f"""
    <h2>1. 📊 执行摘要</h2>
    <div class="metric-grid">
        <div class="metric-card">
            <div style="opacity: 0.9;">总运行数</div>
            <div class="metric-value">{len(self.analyzers)}</div>
        </div>
        <div class="metric-card">
            <div style="opacity: 0.9;">原始Gold</div>
            <div class="metric-value">{total_orig}</div>
        </div>
        <div class="metric-card" style="background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);">
            <div style="opacity: 0.9; color: #000;">严格金标准</div>
            <div class="metric-value" style="color: #000;">{total_gold}</div>
        </div>
        <div class="metric-card">
            <div style="opacity: 0.9;">通过率</div>
            <div class="metric-value">{pass_rate:.1f}%</div>
        </div>
    </div>
"""
        return html
    
    def _section_comparison(self):
        html = """
    <h2>2. 🥇 严格金标准对比</h2>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>原始Gold</th>
                <th>严格Gold</th>
                <th>通过率</th>
                <th>骨架数</th>
                <th>Avg pIC50</th>
                <th>Avg QED</th>
                <th>Avg SA</th>
                <th>毒性组件</th>
            </tr>
        </thead>
        <tbody>
"""
        
        best_gold = max((a.metrics.get('n_gold_strict', 0) for a in self.analyzers), default=0)
        
        for a in self.analyzers:
            n_orig = a.metrics.get('n_gold_original', 0)
            n_strict = a.metrics.get('n_gold_strict', 0)
            rate = (n_strict/n_orig*100) if n_orig > 0 else 0
            
            gold_mols = a.gold_molecules
            if gold_mols:
                pic50_vals = [m['pic50'] for m in gold_mols if not np.isnan(m['pic50'])]
                qed_vals = [m['qed'] for m in gold_mols if not np.isnan(m['qed'])]
                sa_vals = [m['sa'] for m in gold_mols if not np.isnan(m['sa'])]
                
                avg_pic50 = np.mean(pic50_vals) if pic50_vals else 0
                avg_qed = np.mean(qed_vals) if qed_vals else 0
                avg_sa = np.mean(sa_vals) if sa_vals else 0
                
                pic50_str = f"{avg_pic50:.2f}"
                qed_str = f"{avg_qed:.2f}"
                sa_str = f"{avg_sa:.2f}"
            else:
                pic50_str = '—'
                qed_str = '—'
                sa_str = '—'
            
            tox_badge = ''
            if a.config.get('has_toxicity'):
                tox_comps = ', '.join(a.config.get('toxicity_components', []))
                if len(tox_comps) > 20:
                    tox_comps = tox_comps[:17] + '...'
                tox_badge = f'<span class="toxicity-badge">{tox_comps}</span>'
            
            row_class = 'best' if n_strict == best_gold and n_strict > 0 else ''
            
            html += f"""
            <tr class="{row_class}">
                <td><strong>{a.run_name}</strong></td>
                <td>{n_orig}</td>
                <td class="highlight">{n_strict}</td>
                <td>{rate:.1f}%</td>
                <td>{len(a.scaffolds)}</td>
                <td>{pic50_str}</td>
                <td>{qed_str}</td>
                <td>{sa_str}</td>
                <td>{tox_badge if tox_badge else '—'}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    <p style="color: #27ae60; font-weight: 600;">✓ 绿色高亮表示金标准数量最多的run</p>
"""
        return html
    
    def _section_structures(self):
        html = """
    <h2>3. 🔬 Top金标准分子结构</h2>
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px;">
"""
        
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                top3 = sorted(a.gold_molecules, key=lambda x: x['pic50'] if not np.isnan(x['pic50']) else 0, reverse=True)[:3]
                
                html += f"""
        <div style="border: 3px solid #FFD700; border-radius: 10px; padding: 15px; background: #fffef0;">
            <h4 style="margin: 0 0 15px 0;">{a.run_name} (Top 3)</h4>
"""
                
                for i, mol in enumerate(top3):
                    img_html = self._generate_mol_image(mol['smiles'])
                    
                    html += f"""
            <div style="margin: 15px 0; padding: 15px; background: white; border-radius: 8px;">
                <strong>#{i+1}</strong><br>
                {img_html}
                <div style="margin-top: 10px; font-size: 0.85em;">
                    <strong>pIC50:</strong> {mol['pic50']:.2f} | 
                    <strong>QED:</strong> {mol['qed']:.2f} | 
                    <strong>SA:</strong> {mol['sa']:.2f}<br>
                    <strong>MW:</strong> {mol['mw']:.1f} | 
                    <strong>LogP:</strong> {mol['logp']:.2f}
                </div>
                <div style="margin-top: 5px; font-size: 0.7em; word-break: break-all; color: #666;">
                    {mol['smiles'][:60] + '...' if len(mol['smiles']) > 60 else mol['smiles']}
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
    
    def _generate_mol_image(self, smiles):
        try:
            from rdkit import Chem
            from rdkit.Chem import Draw
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(250, 250))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                return f'<img src="data:image/png;base64,{img_str}" style="width: 100%; max-width: 250px;">'
        except:
            pass
        return '<div style="background: #f0f0f0; padding: 40px; text-align: center; color: #999;">需要RDKit</div>'
    
    def _section_toxicity(self):
        runs_with_tox = [a for a in self.analyzers if a.config.get('has_toxicity')]
        runs_without_tox = [a for a in self.analyzers if not a.config.get('has_toxicity')]
        
        html = """
    <h2>4. 🧪 毒性组件影响分析</h2>
"""
        
        if runs_with_tox and runs_without_tox:
            avg_with = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_with_tox])
            avg_without = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_without_tox])
            
            impact_pct = (avg_with-avg_without)/avg_without*100 if avg_without > 0 else 0
            
            html += f"""
    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0;">
        <div style="background: #fff3cd; padding: 20px; border-radius: 10px;">
            <h3>含毒性组件 ({len(runs_with_tox)}个)</h3>
            <p style="font-size: 2em; font-weight: bold; color: #856404;">{avg_with:.1f}</p>
            <p>平均严格金标准数</p>
        </div>
        <div style="background: #d4edda; padding: 20px; border-radius: 10px;">
            <h3>不含毒性组件 ({len(runs_without_tox)}个)</h3>
            <p style="font-size: 2em; font-weight: bold; color: #155724;">{avg_without:.1f}</p>
            <p>平均严格金标准数</p>
        </div>
    </div>
    <p style="background: #e3f2fd; padding: 15px; border-radius: 8px;">
        <strong>💡 结论:</strong> 毒性组件使金标准数量变化 <strong>{impact_pct:+.1f}%</strong>
    </p>
"""
        
        return html
    
    def _section_scaffolds(self):
        all_scaffolds = set()
        for a in self.analyzers:
            all_scaffolds.update(a.scaffolds)
        
        html = f"""
    <h2>5. 🧩 核心骨架多样性</h2>
    <p style="background: #f8f9fa; padding: 15px; border-radius: 8px;">
        <strong>总计:</strong> {len(all_scaffolds)} 个不同的Murcko骨架
    </p>
    
    <table>
        <thead>
            <tr><th>Run</th><th>骨架数</th><th>金标准数</th><th>多样性比</th></tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            n_gold = a.metrics.get('n_gold_strict', 0)
            n_scaffolds = len(a.scaffolds)
            ratio = (n_scaffolds / n_gold) if n_gold > 0 else 0
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td class="highlight">{n_scaffolds}</td>
                <td>{n_gold}</td>
                <td>{ratio:.2f}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_config(self):
        html = """
    <h2>6. ⚙️ 配置详细对比</h2>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>Prior</th>
                <th>Agent</th>
                <th>Batch</th>
                <th>LR</th>
                <th>Sigma</th>
                <th>组件数</th>
            </tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            prior = os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A'
            agent = os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A'
            lr = a.config.get('learning_rate')
            lr_str = f'{lr:.5f}' if lr else 'N/A'
            
            html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{prior}</td>
                <td>{agent}</td>
                <td>{a.config.get('batch_size', 'N/A')}</td>
                <td>{lr_str}</td>
                <td>{a.config.get('sigma', 'N/A')}</td>
                <td>{len(a.config.get('components', []))}</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
"""
        return html
    
    def _section_components(self):
        all_comps = set()
        for a in self.analyzers:
            all_comps.update(a.config.get('components', []))
        
        html = """
    <h2>7. 🎯 评分组件权重对比</h2>
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
    
    def _html_footer(self):
        return """
    <div style="margin-top: 60px; padding-top: 30px; border-top: 3px solid #3498db; text-align: center; color: #7f8c8d;">
        <p><strong>REINVENT4 Complete Analysis Report</strong></p>
        <p style="font-size: 0.9em; margin-top: 10px;">
            Gold Standard Filtering, Toxicity Analysis & Scaffold Diversity
        </p>
    </div>
</div>
</body>
</html>
"""


def main():
    print("="*80)
    print("REINVENT4 完整终极分析")
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
        # 生成HTML报告
        reporter = HTMLReportGenerator(analyzers)
        html_file = reporter.generate()
        
        # 生成CSV摘要
        print("\n📊 生成CSV摘要...")
        summary_data = []
        for a in analyzers:
            gold_mols = a.gold_molecules
            if gold_mols:
                pic50_vals = [m['pic50'] for m in gold_mols if not np.isnan(m['pic50'])]
                qed_vals = [m['qed'] for m in gold_mols if not np.isnan(m['qed'])]
                sa_vals = [m['sa'] for m in gold_mols if not np.isnan(m['sa'])]
                
                avg_pic50 = np.mean(pic50_vals) if pic50_vals else 0
                avg_qed = np.mean(qed_vals) if qed_vals else 0
                avg_sa = np.mean(sa_vals) if sa_vals else 0
                
                pic50_str = f"{avg_pic50:.2f}"
                qed_str = f"{avg_qed:.2f}"
                sa_str = f"{avg_sa:.2f}"
            else:
                pic50_str = '—'
                qed_str = '—'
                sa_str = '—'
            
            summary_data.append({
                'Run': a.run_name,
                'Original_Gold': a.metrics.get('n_gold_original', 0),
                'Strict_Gold': a.metrics.get('n_gold_strict', 0),
                'Scaffolds': len(a.scaffolds),
                'Pass_Rate_%': f"{(a.metrics.get('n_gold_strict', 0)/a.metrics.get('n_gold_original', 1)*100):.1f}",
                'Toxicity': 'Yes' if a.config.get('has_toxicity') else 'No',
                'Tox_Components': ', '.join(a.config.get('toxicity_components', [])) if a.config.get('toxicity_components') else '—',
                'Avg_pIC50': pic50_str,
                'Avg_QED': qed_str,
                'Avg_SA': sa_str,
                'Prior': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
                'Agent': os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A',
                'Batch_Size': a.config.get('batch_size', 'N/A'),
                'Learning_Rate': a.config.get('learning_rate', 'N/A'),
                'Sigma': a.config.get('sigma', 'N/A'),
                'N_Components': len(a.config.get('components', [])),
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv('complete_analysis_summary.csv', index=False)
        print("✅ CSV摘要: complete_analysis_summary.csv")
        
        # 打印快速摘要
        print("\n" + "="*80)
        print("📋 快速摘要")
        print("="*80)
        print("\n" + df[['Run', 'Original_Gold', 'Strict_Gold', 'Scaffolds', 'Toxicity']].to_string(index=False))
        
        # 找出最佳run
        if df['Strict_Gold'].astype(int).max() > 0:
            best_idx = df['Strict_Gold'].astype(int).idxmax()
            best_run = df.loc[best_idx]
            print(f"\n🏆 最佳Run: {best_run['Run']}")
            print(f"   严格金标准: {best_run['Strict_Gold']}")
            print(f"   骨架多样性: {best_run['Scaffolds']}")
            print(f"   毒性组件: {best_run['Toxicity']}")
        
        # 毒性影响分析
        with_tox = df[df['Toxicity'] == 'Yes']
        without_tox = df[df['Toxicity'] == 'No']
        
        if len(with_tox) > 0 and len(without_tox) > 0:
            avg_with = with_tox['Strict_Gold'].astype(int).mean()
            avg_without = without_tox['Strict_Gold'].astype(int).mean()
            
            print(f"\n🧪 毒性组件影响:")
            print(f"   含毒性组件: 平均 {avg_with:.1f} 个严格金标准")
            print(f"   不含毒性组件: 平均 {avg_without:.1f} 个严格金标准")
            if avg_without > 0:
                print(f"   影响: {(avg_with-avg_without)/avg_without*100:+.1f}%")
    
    print("\n" + "="*80)
    print("✅ 分析完成！")
    print("="*80)
    print(f"\n📁 生成的文件:")
    print(f"  1. complete_runs_report.html - 交互式HTML报告")
    print(f"  2. complete_analysis_summary.csv - 详细数据表")
    print(f"\n🌐 在浏览器中打开 complete_runs_report.html 查看完整报告！")
    print("   (Windows: explorer.exe complete_runs_report.html)")
    print("="*80)

if __name__ == "__main__":
    main()
