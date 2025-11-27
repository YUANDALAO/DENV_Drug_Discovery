#!/usr/bin/env python3
"""
REINVENT4 最终工作版本
使用原始的宽松金标准定义（匹配unified_candidate_analysis.py）
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
            
            # 检测毒性组件（更精确）
            toxicity_names = {'Toxicity_Alerts', 'PAINS_Filter', 'Metabolic_Stability'}
            for comp_detail in components:
                display_name = comp_detail.get('display_name', '')
                if display_name in toxicity_names:
                    config['has_toxicity'] = True
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
        """加载数据"""
        print(f"  分析 {self.run_name}...")
        
        # 加载配置
        config_path = os.path.join(self.run_dir, 'config.toml')
        if os.path.exists(config_path):
            self.config = ConfigParser.parse_toml(config_path)
            tox_info = ""
            if self.config.get('has_toxicity'):
                tox_comps = self.config.get('toxicity_components', [])
                tox_info = f" [🧪 毒性: {', '.join(tox_comps)}]"
            print(f"    ✓ 配置: {len(self.config['components'])} 个组件{tox_info}")
        
        # 加载CSV文件
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
                'qsar': ['DENV_Activity (raw)', 'DENV_Activity', 'QSAR_Score'],
                'qed': ['QED (raw)', 'QED', 'qed'],
                'mw': ['MW (raw)', 'MW', 'Molecular weight'],
                'sa': ['SA (raw)', 'SA', 'SAScore'],
                'logp': ['LogP (raw)', 'LogP', 'SlogP (RDKit)'],
                'tpsa': ['TPSA (raw)', 'TPSA', 'tpsa'],
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
        """使用unified_candidate_analysis.py的金标准"""
        df = self.data.get('gold')
        if df is None or len(df) == 0:
            self.metrics['n_gold_strict'] = 0
            return
        
        df = df.copy()
        df.columns = df.columns.str.strip()
        
        # 提取各个指标（优先使用raw值）
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
        
        # 金标准: pIC50≥8.0, QED≥0.7, SA≤4.0, MW 300-500, LogP 1-4
        mask = (
            (pic50 >= 8.0) &
            (qed >= 0.7) &
            (sa <= 4.0) &
            (mw >= 300) & (mw <= 500) &
            (logp >= 1) & (logp <= 4)
        )
        
        gold_strict = df[mask]
        self.metrics['n_gold_strict'] = len(gold_strict)
        
        # 保存金标准分子
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
            pass


def main():
    print("="*80)
    print("REINVENT4 最终工作版本分析")
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
    
    # 生成摘要
    if len(analyzers) > 0:
        print("\n" + "="*80)
        print("📊 金标准分子摘要 (pIC50≥8, QED≥0.7, SA≤4, MW 300-500, LogP 1-4)")
        print("="*80)
        
        summary_data = []
        for a in analyzers:
            gold_mols = a.gold_molecules
            avg_pic50 = np.mean([m['pic50'] for m in gold_mols if not np.isnan(m['pic50'])]) if gold_mols else 0
            avg_qed = np.mean([m['qed'] for m in gold_mols if not np.isnan(m['qed'])]) if gold_mols else 0
            avg_sa = np.mean([m['sa'] for m in gold_mols if not np.isnan(m['sa'])]) if gold_mols else 0
            
            has_tox = '✓' if a.config.get('has_toxicity') else '—'
            tox_list = ', '.join(a.config.get('toxicity_components', []))[:30] if a.config.get('toxicity_components') else '—'
            
            summary_data.append({
                'Run': a.run_name,
                'Original_Gold': a.metrics.get('n_gold_original', 0),
                'Strict_Gold': a.metrics.get('n_gold_strict', 0),
                'Scaffolds': len(a.scaffolds),
                'Toxicity': has_tox,
                'Avg_pIC50': f"{avg_pic50:.2f}" if gold_mols else '—',
                'Avg_QED': f"{avg_qed:.2f}" if gold_mols else '—',
                'Avg_SA': f"{avg_sa:.2f}" if gold_mols else '—',
                'Prior': os.path.basename(a.config.get('prior_file', 'N/A')) if a.config.get('prior_file') else 'N/A',
                'Agent': os.path.basename(a.config.get('agent_file', 'N/A')) if a.config.get('agent_file') else 'N/A',
            })
        
        df = pd.DataFrame(summary_data)
        df.to_csv('gold_standard_summary.csv', index=False)
        
        print("\n" + df[['Run', 'Original_Gold', 'Strict_Gold', 'Scaffolds', 'Toxicity']].to_string(index=False))
        
        print(f"\n🏆 最佳Run: {df.loc[df['Strict_Gold'].idxmax(), 'Run']} "
              f"({df['Strict_Gold'].max()} 个严格金标准)")
        
        # 毒性分析
        with_tox = df[df['Toxicity'] == '✓']
        without_tox = df[df['Toxicity'] == '—']
        
        if len(with_tox) > 0 and len(without_tox) > 0:
            avg_with = with_tox['Strict_Gold'].mean()
            avg_without = without_tox['Strict_Gold'].mean()
            print(f"\n🧪 毒性组件影响:")
            print(f"   含毒性组件: 平均 {avg_with:.1f} 个严格金标准")
            print(f"   不含毒性组件: 平均 {avg_without:.1f} 个严格金标准")
        
        print(f"\n📁 详细结果: gold_standard_summary.csv")
    
    print("\n" + "="*80)

if __name__ == "__main__":
    main()
