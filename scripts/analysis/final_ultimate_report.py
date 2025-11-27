#!/usr/bin/env python3
"""
REINVENT4 终极HTML报告 - 修复版
兼容不同RDKit版本
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
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 100

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


class MolecularFeatureAnalyzer:
    """分子特征分析器 - 兼容多版本RDKit"""
    
    @staticmethod
    def calculate_features(df):
        """计算分子特征"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen
            from rdkit.Chem.Scaffolds import MurckoScaffold
        except ImportError:
            print("    ⚠️  RDKit未安装，跳过分子特征计算")
            return None
        
        features = []
        
        for idx, row in df.iterrows():
            smiles = row.get('SMILES', '')
            if not smiles:
                continue
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            try:
                feat = {
                    'SMILES': smiles,
                    'MW': Descriptors.MolWt(mol),
                    'LogP': Crippen.MolLogP(mol),
                    'TPSA': Descriptors.TPSA(mol),
                    'HBA': Lipinski.NumHAcceptors(mol),
                    'HBD': Lipinski.NumHDonors(mol),
                    'RotBonds': Descriptors.NumRotatableBonds(mol),
                    'AromaticRings': Descriptors.NumAromaticRings(mol),
                    'HeavyAtoms': Lipinski.HeavyAtomCount(mol),
                    'NumRings': Descriptors.RingCount(mol),
                }
                
                # FractionCsp3 - 兼容不同版本
                try:
                    # 新版本在 Lipinski 模块
                    feat['FractionCSP3'] = Lipinski.FractionCsp3(mol)
                except AttributeError:
                    try:
                        # 旧版本在 Descriptors 模块
                        feat['FractionCSP3'] = Descriptors.FractionCsp3(mol)
                    except:
                        # 手动计算
                        num_csp3 = sum(1 for atom in mol.GetAtoms() if atom.GetHybridization() == Chem.HybridizationType.SP3 and atom.GetSymbol() == 'C')
                        num_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
                        feat['FractionCSP3'] = num_csp3 / num_c if num_c > 0 else 0
                
                # 从CSV提取数据
                for col in ['DENV_Activity (raw)', 'DENV_Activity', 'pIC50']:
                    if col in row:
                        feat['pIC50'] = row[col]
                        break
                
                for col in ['QED (raw)', 'QED']:
                    if col in row:
                        feat['QED'] = row[col]
                        break
                
                for col in ['SA (raw)', 'SA']:
                    if col in row:
                        feat['SA'] = row[col]
                        break
                
                # 毒性相关基团
                feat['HasNitro'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[N+](=O)[O-]'))
                feat['HasAzo'] = mol.HasSubstructMatch(Chem.MolFromSmarts('N=N'))
                feat['HasHalogen2'] = mol.HasSubstructMatch(Chem.MolFromSmarts('[F,Cl,Br,I][C,c][F,Cl,Br,I]'))
                feat['HasQuinone'] = mol.HasSubstructMatch(Chem.MolFromSmarts('C1=CC(=O)C=CC1=O'))
                
                toxicity_flags = [feat['HasNitro'], feat['HasAzo'], feat['HasHalogen2'], feat['HasQuinone']]
                feat['ToxicityRisk'] = sum(toxicity_flags)
                
                features.append(feat)
            except Exception as e:
                # 跳过有问题的分子
                continue
        
        return pd.DataFrame(features) if features else None


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
        self.molecular_features = None
    
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
        self._analyze_molecular_features()
    
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
    
    def _analyze_molecular_features(self):
        """分析分子特征"""
        df = self.data.get('gold')
        if df is not None and len(df) > 0:
            self.molecular_features = MolecularFeatureAnalyzer.calculate_features(df)
            if self.molecular_features is not None:
                print(f"    ✓ 分子特征: {len(self.molecular_features)} 个分子")


class HTMLReportGenerator:
    """生成HTML报告"""
    
    def __init__(self, analyzers):
        self.analyzers = analyzers
        self.timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    def generate(self, output_file='ultimate_runs_report.html'):
        print("\n📄 生成终极HTML报告...")
        
        html = self._html_header()
        html += self._section_gold_standard()
        html += self._section_summary()
        html += self._section_comparison()
        html += self._section_molecular_features()
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
    <title>REINVENT4 Ultimate Analysis</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', 'DejaVu Sans', Arial, sans-serif; line-height: 1.6; 
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 10px; }}
        .container {{ max-width: 1600px; margin: 0 auto; background: white; 
                      padding: 20px; border-radius: 12px; box-shadow: 0 10px 40px rgba(0,0,0,0.2); }}
        h1 {{ color: #2c3e50; border-bottom: 4px solid #3498db; padding-bottom: 10px; font-size: 1.8em; }}
        h2 {{ color: #34495e; margin-top: 30px; border-left: 5px solid #3498db; 
              padding-left: 15px; font-size: 1.4em; }}
        h3 {{ color: #555; margin-top: 20px; font-size: 1.2em; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.85em; }}
        th {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
              color: white; padding: 10px; text-align: left; font-size: 0.9em; }}
        td {{ padding: 8px; border-bottom: 1px solid #ddd; }}
        tr:hover {{ background: #f5f5f5; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); 
                        gap: 15px; margin: 20px 0; }}
        .metric-card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                        color: white; padding: 18px; border-radius: 10px; text-align: center; }}
        .metric-value {{ font-size: 2.2em; font-weight: bold; margin: 8px 0; }}
        .gold-box {{ background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%); 
                     padding: 20px; border-radius: 10px; margin: 20px 0; color: #000; }}
        .criterion {{ background: white; padding: 10px; margin: 5px 0; border-radius: 5px; font-weight: 600; }}
        .toxicity-badge {{ background: #e74c3c; color: white; padding: 3px 8px; 
                           border-radius: 12px; font-size: 0.8em; font-weight: 600; }}
        .highlight {{ background: #fffacd; font-weight: bold; }}
        .best {{ background: #90EE90; font-weight: bold; }}
        .feature-table {{ font-size: 0.8em; }}
        .feature-table th {{ font-size: 0.85em; padding: 8px; }}
        .feature-table td {{ padding: 6px; }}
        @media (max-width: 768px) {{
            .container {{ padding: 12px; }}
            h1 {{ font-size: 1.4em; }}
            table {{ font-size: 0.7em; }}
        }}
    </style>
</head>
<body>
<div class="container">
    <h1>🧬 REINVENT4 Ultimate Multi-Run Analysis</h1>
    <p style="text-align: center; color: #7f8c8d; margin-bottom: 20px; font-size: 0.95em;">
        Generated: {self.timestamp} | Complete Molecular Feature Analysis
    </p>
"""
    
    def _section_gold_standard(self):
        return """
    <div class="gold-box">
        <h3 style="margin-top: 0;">⭐ Gold Standard Criteria</h3>
        <p style="margin-bottom: 12px;">Must satisfy <strong>ALL</strong> conditions:</p>
        <div class="criterion">✓ pIC50 >= 8.0 (IC50 <= 10 nM)</div>
        <div class="criterion">✓ QED >= 0.7 (Drug-likeness)</div>
        <div class="criterion">✓ SA <= 4.0 (Synthetic Accessibility)</div>
        <div class="criterion">✓ MW: 300-500 Da</div>
        <div class="criterion">✓ LogP: 1-4</div>
    </div>
"""
    
    def _section_summary(self):
        total_mols = sum(a.metrics.get('n_total', 0) for a in self.analyzers)
        total_gold = sum(a.metrics.get('n_gold_strict', 0) for a in self.analyzers)
        
        success_rate = (total_gold / total_mols * 100) if total_mols > 0 else 0
        
        html = f"""
    <h2>1. 📊 Executive Summary</h2>
    <div class="metric-grid">
        <div class="metric-card">
            <div style="opacity: 0.9; font-size: 0.9em;">Total Runs</div>
            <div class="metric-value">{len(self.analyzers)}</div>
        </div>
        <div class="metric-card">
            <div style="opacity: 0.9; font-size: 0.9em;">Total Molecules</div>
            <div class="metric-value">{total_mols:,}</div>
        </div>
        <div class="metric-card" style="background: linear-gradient(135deg, #FFD700 0%, #FFA500 100%);">
            <div style="opacity: 0.9; font-size: 0.9em; color: #000;">Gold Standard</div>
            <div class="metric-value" style="color: #000;">{total_gold}</div>
        </div>
        <div class="metric-card">
            <div style="opacity: 0.9; font-size: 0.9em;">Success Rate</div>
            <div class="metric-value">{success_rate:.3f}%</div>
        </div>
    </div>
"""
        return html
    
    def _section_comparison(self):
        html = """
    <h2>2. 🥇 Gold Standard Comparison</h2>
    <table>
        <thead>
            <tr>
                <th>Run</th>
                <th>Total Molecules</th>
                <th>Gold Standard</th>
                <th>Success Rate</th>
                <th>Scaffolds</th>
                <th>Avg pIC50</th>
                <th>Avg QED</th>
                <th>Avg SA</th>
                <th>Toxicity</th>
            </tr>
        </thead>
        <tbody>
"""
        
        best_gold = max((a.metrics.get('n_gold_strict', 0) for a in self.analyzers), default=0)
        
        for a in self.analyzers:
            n_total = a.metrics.get('n_total', 0)
            n_gold = a.metrics.get('n_gold_strict', 0)
            success_rate = (n_gold / n_total * 100) if n_total > 0 else 0
            
            gold_mols = a.gold_molecules
            if gold_mols:
                pic50_vals = [m['pic50'] for m in gold_mols if not np.isnan(m['pic50'])]
                qed_vals = [m['qed'] for m in gold_mols if not np.isnan(m['qed'])]
                sa_vals = [m['sa'] for m in gold_mols if not np.isnan(m['sa'])]
                
                pic50_str = f"{np.mean(pic50_vals):.2f}" if pic50_vals else '—'
                qed_str = f"{np.mean(qed_vals):.2f}" if qed_vals else '—'
                sa_str = f"{np.mean(sa_vals):.2f}" if sa_vals else '—'
            else:
                pic50_str = qed_str = sa_str = '—'
            
            tox_badge = ''
            if a.config.get('has_toxicity'):
                tox_comps = ', '.join(a.config.get('toxicity_components', []))
                if len(tox_comps) > 15:
                    tox_comps = tox_comps[:12] + '...'
                tox_badge = f'<span class="toxicity-badge">{tox_comps}</span>'
            
            row_class = 'best' if n_gold == best_gold and n_gold > 0 else ''
            
            html += f"""
            <tr class="{row_class}">
                <td><strong>{a.run_name}</strong></td>
                <td>{n_total:,}</td>
                <td class="highlight">{n_gold}</td>
                <td>{success_rate:.4f}%</td>
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
    <p style="color: #27ae60; font-weight: 600; font-size: 0.9em;">✓ Green highlight = Highest gold standard count</p>
"""
        return html
    
    def _section_molecular_features(self):
        """深度分子特征对比"""
        html = """
    <h2>3. 🔬 Molecular Feature Analysis (Mean ± SD)</h2>
    <p style="background: #e3f2fd; padding: 12px; border-radius: 6px; font-size: 0.9em;">
        Detailed molecular properties of gold standard candidates from each run
    </p>
    
    <table class="feature-table">
        <thead>
            <tr>
                <th>Run</th>
                <th>pIC50</th>
                <th>QED</th>
                <th>SA</th>
                <th>MW</th>
                <th>LogP</th>
                <th>TPSA</th>
                <th>HBA</th>
                <th>HBD</th>
                <th>RotBonds</th>
                <th>CSP3</th>
                <th>ToxRisk</th>
            </tr>
        </thead>
        <tbody>
"""
        
        for a in self.analyzers:
            if a.molecular_features is not None and len(a.molecular_features) > 0:
                mf = a.molecular_features
                
                def format_stat(col):
                    if col in mf.columns:
                        vals = mf[col].dropna()
                        if len(vals) > 0:
                            return f"{vals.mean():.2f}±{vals.std():.2f}"
                    return '—'
                
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{format_stat('pIC50')}</td>
                <td>{format_stat('QED')}</td>
                <td>{format_stat('SA')}</td>
                <td>{format_stat('MW')}</td>
                <td>{format_stat('LogP')}</td>
                <td>{format_stat('TPSA')}</td>
                <td>{format_stat('HBA')}</td>
                <td>{format_stat('HBD')}</td>
                <td>{format_stat('RotBonds')}</td>
                <td>{format_stat('FractionCSP3')}</td>
                <td>{format_stat('ToxicityRisk')}</td>
            </tr>
"""
            else:
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td colspan="11" style="text-align: center; color: #999;">No feature data</td>
            </tr>
"""
        
        html += """
        </tbody>
    </table>
    <p style="font-size: 0.85em; color: #666; margin-top: 10px;">
        <strong>Note:</strong> CSP3 = Fraction of sp³ carbons; ToxRisk = Toxicity risk score (0-4)
    </p>
"""
        return html
    
    def _section_structures(self):
        html = """
    <h2>4. 💎 Top Gold Standard Structures</h2>
    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 15px;">
"""
        
        for a in self.analyzers:
            if len(a.gold_molecules) > 0:
                top3 = sorted(a.gold_molecules, key=lambda x: x['pic50'] if not np.isnan(x['pic50']) else 0, reverse=True)[:3]
                
                html += f"""
        <div style="border: 3px solid #FFD700; border-radius: 8px; padding: 12px; background: #fffef0;">
            <h4 style="margin: 0 0 12px 0; font-size: 1em;">{a.run_name} (Top 3)</h4>
"""
                
                for i, mol in enumerate(top3):
                    img_html = self._generate_mol_image(mol['smiles'])
                    
                    html += f"""
            <div style="margin: 12px 0; padding: 12px; background: white; border-radius: 6px;">
                <strong>#{i+1}</strong><br>
                {img_html}
                <div style="margin-top: 8px; font-size: 0.8em;">
                    <strong>pIC50:</strong> {mol['pic50']:.2f} | 
                    <strong>QED:</strong> {mol['qed']:.2f} | 
                    <strong>SA:</strong> {mol['sa']:.2f}<br>
                    <strong>MW:</strong> {mol['mw']:.1f} |
                    <strong>LogP:</strong> {mol['logp']:.2f}
                </div>
                <div style="margin-top: 5px; font-size: 0.7em; word-break: break-all; color: #666;">
                    {mol['smiles'][:50] + '...' if len(mol['smiles']) > 50 else mol['smiles']}
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
                img = Draw.MolToImage(mol, size=(230, 230))
                buffered = BytesIO()
                img.save(buffered, format="PNG")
                img_str = base64.b64encode(buffered.getvalue()).decode()
                return f'<img src="data:image/png;base64,{img_str}" style="width: 100%; max-width: 230px;">'
        except:
            pass
        return '<div style="background: #f0f0f0; padding: 35px; text-align: center; color: #999; font-size: 0.85em;">RDKit required</div>'
    
    def _section_toxicity(self):
        runs_with_tox = [a for a in self.analyzers if a.config.get('has_toxicity')]
        runs_without_tox = [a for a in self.analyzers if not a.config.get('has_toxicity')]
        
        html = """
    <h2>5. 🧪 Toxicity Component Impact Analysis</h2>
"""
        
        if runs_with_tox and runs_without_tox:
            avg_with = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_with_tox])
            avg_without = np.mean([a.metrics.get('n_gold_strict', 0) for a in runs_without_tox])
            
            impact_pct = (avg_with-avg_without)/avg_without*100 if avg_without > 0 else 0
            
            html += f"""
    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 18px; margin: 20px 0;">
        <div style="background: #fff3cd; padding: 18px; border-radius: 8px;">
            <h3 style="margin-top: 0; font-size: 1.1em;">With Toxicity ({len(runs_with_tox)} runs)</h3>
            <p style="font-size: 1.8em; font-weight: bold; color: #856404; margin: 8px 0;">{avg_with:.1f}</p>
            <p style="font-size: 0.9em;">Avg gold standard molecules</p>
        </div>
        <div style="background: #d4edda; padding: 18px; border-radius: 8px;">
            <h3 style="margin-top: 0; font-size: 1.1em;">Without Toxicity ({len(runs_without_tox)} runs)</h3>
            <p style="font-size: 1.8em; font-weight: bold; color: #155724; margin: 8px 0;">{avg_without:.1f}</p>
            <p style="font-size: 0.9em;">Avg gold standard molecules</p>
        </div>
    </div>
    <p style="background: #e3f2fd; padding: 12px; border-radius: 6px; font-size: 0.9em;">
        <strong>💡 Conclusion:</strong> Toxicity components changed gold standard count by <strong>{impact_pct:+.1f}%</strong>
    </p>
"""
        
        # 显示哪些runs有毒性组件
        if runs_with_tox:
            html += """
    <h3 style="font-size: 1.1em;">Runs with Toxicity Components:</h3>
    <table class="feature-table">
        <thead>
            <tr><th>Run</th><th>Toxicity Components</th><th>Gold Count</th></tr>
        </thead>
        <tbody>
"""
            for a in runs_with_tox:
                tox_comps = ', '.join(a.config.get('toxicity_components', []))
                html += f"""
            <tr>
                <td><strong>{a.run_name}</strong></td>
                <td>{tox_comps}</td>
                <td>{a.metrics.get('n_gold_strict', 0)}</td>
            </tr>
"""
            html += """
        </tbody>
    </table>
"""
        
        return html
    
    def _section_scaffolds(self):
        all_scaffolds = set()
        for a in self.analyzers:
            all_scaffolds.update(a.scaffolds)
        
        html = f"""
    <h2>6. 🧩 Scaffold Diversity Analysis</h2>
    <p style="background: #f8f9fa; padding: 12px; border-radius: 6px; font-size: 0.9em;">
        <strong>Total:</strong> {len(all_scaffolds)} unique Murcko scaffolds across all runs
    </p>
    
    <table class="feature-table">
        <thead>
            <tr><th>Run</th><th>Unique Scaffolds</th><th>Gold Count</th><th>Diversity Ratio</th></tr>
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
    <h2>7. ⚙️ Configuration Comparison</h2>
    <table class="feature-table">
        <thead>
            <tr>
                <th>Run</th>
                <th>Prior Model</th>
                <th>Agent Model</th>
                <th>Batch</th>
                <th>Learning Rate</th>
                <th>Sigma</th>
                <th>Components</th>
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
    <h2>8. 🎯 Scoring Component Weights</h2>
    <div style="overflow-x: auto;">
    <table class="feature-table">
        <thead>
            <tr><th>Component</th>
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
    <div style="margin-top: 50px; padding-top: 25px; border-top: 3px solid #3498db; text-align: center; color: #7f8c8d;">
        <p style="font-size: 1.1em;"><strong>REINVENT4 Ultimate Analysis Report</strong></p>
        <p style="font-size: 0.9em; margin-top: 8px;">
            Complete Molecular Feature Analysis | Gold Standard Filtering | Toxicity Impact
        </p>
        <p style="font-size: 0.85em; margin-top: 12px; color: #95a5a6;">
            📱 Mobile & Desktop Optimized | 🔬 Research Grade Analysis
        </p>
    </div>
</div>
</body>
</html>
"""


def main():
    print("="*80)
    print("REINVENT4 终极HTML报告生成器")
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
                
                pic50_str = f"{np.mean(pic50_vals):.2f}±{np.std(pic50_vals):.2f}" if pic50_vals else '—'
                qed_str = f"{np.mean(qed_vals):.2f}±{np.std(qed_vals):.2f}" if qed_vals else '—'
                sa_str = f"{np.mean(sa_vals):.2f}±{np.std(sa_vals):.2f}" if sa_vals else '—'
            else:
                pic50_str = qed_str = sa_str = '—'
            
            n_total = a.metrics.get('n_total', 0)
            n_gold = a.metrics.get('n_gold_strict', 0)
            success_rate = (n_gold / n_total * 100) if n_total > 0 else 0
            
            summary_data.append({
                'Run': a.run_name,
                'Total_Molecules': n_total,
                'Gold_Standard': n_gold,
                'Success_Rate_%': f"{success_rate:.4f}",
                'Scaffolds': len(a.scaffolds),
                'Toxicity': 'Yes' if a.config.get('has_toxicity') else 'No',
                'Tox_Components': '; '.join(a.config.get('toxicity_components', [])) if a.config.get('toxicity_components') else '—',
                'pIC50': pic50_str,
                'QED': qed_str,
                'SA': sa_str,
                'Prior': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
                'Agent': os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A',
                'Batch_Size': a.config.get('batch_size', 'N/A'),
                'Learning_Rate': a.config.get('learning_rate', 'N/A'),
                'Sigma': a.config.get('sigma', 'N/A'),
                'N_Components': len(a.config.get('components', [])),
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv('ultimate_analysis_summary.csv', index=False)
        print("✅ CSV摘要: ultimate_analysis_summary.csv")
        
        # 打印快速摘要
        print("\n" + "="*80)
        print("📋 快速摘要")
        print("="*80)
        print("\n" + df[['Run', 'Total_Molecules', 'Gold_Standard', 'Success_Rate_%', 'Toxicity']].to_string(index=False))
        
        # 最佳run
        if df['Gold_Standard'].astype(int).max() > 0:
            best_idx = df['Gold_Standard'].astype(int).idxmax()
            best_run = df.loc[best_idx]
            print(f"\n🏆 Best Run: {best_run['Run']}")
            print(f"   Gold Standard: {best_run['Gold_Standard']}")
            print(f"   Success Rate: {best_run['Success_Rate_%']}%")
            print(f"   Toxicity Components: {best_run['Tox_Components']}")
    
    print("\n" + "="*80)
    print("✅ 分析完成！")
    print("="*80)
    print(f"\n📁 Generated Files:")
    print(f"  1. ultimate_runs_report.html - Interactive HTML report")
    print(f"  2. ultimate_analysis_summary.csv - Data summary")
    print(f"\n🌐 Open ultimate_runs_report.html in browser!")
    print("   (Windows: explorer.exe ultimate_runs_report.html)")
    print("="*80)

if __name__ == "__main__":
    main()


# HTMLReportGenerator类保持不变，使用之前的代码
# (太长了，这里省略，直接使用前面定义的HTML生成器)

# ... 继续使用之前的HTMLReportGenerator和main函数 ...
