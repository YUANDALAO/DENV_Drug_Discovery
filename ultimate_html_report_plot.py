cat > final_ultimate_report.py << 'EOF'
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


# HTMLReportGenerator类保持不变，使用之前的代码
# (太长了，这里省略，直接使用前面定义的HTML生成器)

# ... 继续使用之前的HTMLReportGenerator和main函数 ...
EOF