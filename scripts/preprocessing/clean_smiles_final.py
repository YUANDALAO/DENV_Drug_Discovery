import pandas as pd
from rdkit import Chem
import re

def tokenize_smiles(smiles):
    """简单的SMILES tokenization"""
    # 匹配常见的SMILES tokens
    pattern = r'\[[^\]]+\]|Br|Cl|[BCNOSPFIHbcnops()]|\d|[=#\-+\\\/]|%\d{2}'
    tokens = re.findall(pattern, smiles)
    return set(tokens)

def is_strictly_compatible(smiles):
    """严格检查是否只包含REINVENT支持的token"""
    # REINVENT确切支持的token集合
    allowed_tokens = {
        'c', 'N', '-', 's', '3', '4', '=', '%10', 'O', 'F', '8', 'S', '6', ')', 
        '9', 'Cl', '^', 'n', '[O-]', 'o', '#', '[S+]', '[nH]', '7', '(', '[n+]', 
        'C', '[N+]', '[N-]', '1', 'Br', '5', '2', '$'
    }
    
    # 简单检查：避免复杂的token
    forbidden_substrings = [
        '[n-]', '[Na+]', '[K+]', '[Cl-]', '[Br-]', '[I-]', '[F-]',
        '[Ca+2]', '[Mg+2]', '[Fe+2]', '[Fe+3]', '[Zn+2]',
        '|', '.'  # 避免混合物和盐
    ]
    
    for forbidden in forbidden_substrings:
        if forbidden in smiles:
            return False
    
    # 检查长度和复杂度
    if len(smiles) > 200:  # 避免过长的SMILES
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # 原子数限制
        if mol.GetNumAtoms() > 60:
            return False
        
        # 只允许常见原子
        atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        common_atoms = {'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'H'}
        if not atoms.issubset(common_atoms):
            return False
        
        return True
    except:
        return False

def final_clean_smiles(input_file, output_file):
    """最终的SMILES清理"""
    df = pd.read_csv(input_file)
    
    print(f"开始最终清理，原始数据: {len(df)} 个化合物")
    
    valid_smiles = []
    rejected = []
    
    for i, smiles in enumerate(df['Smiles']):
        if pd.isna(smiles):
            continue
            
        smiles_str = str(smiles).strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                if is_strictly_compatible(canonical_smiles):
                    valid_smiles.append(canonical_smiles)
                else:
                    rejected.append(canonical_smiles[:60] + "...")
                    if len(rejected) <= 5:  # 只显示前5个被拒绝的
                        print(f"拒绝: {canonical_smiles[:60]}...")
        except Exception as e:
            rejected.append(f"解析错误: {smiles_str[:30]}...")
    
    # 去重
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    # 保存
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"\n最终清理结果:")
    print(f"  有效SMILES: {len(unique_smiles)}")
    print(f"  被拒绝: {len(rejected)}")
    print(f"  成功率: {len(unique_smiles)/(len(unique_smiles)+len(rejected))*100:.1f}%")
    print(f"  保存到: {output_file}")

if __name__ == "__main__":
    final_clean_smiles("data/NS3.csv", "data/denv_smiles_final.tsv")
