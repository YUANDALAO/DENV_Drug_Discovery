import pandas as pd
from rdkit import Chem
import re

def is_compatible_with_reinvent(smiles):
    """严格检查SMILES是否与REINVENT模型兼容"""
    
    # 检查是否包含任何离子
    forbidden_patterns = [
        '[Cl-]', '[Br-]', '[I-]', '[F-]',  # 阴离子
        '[Na+]', '[K+]', '[Ca+2]', '[Mg+2]',  # 阳离子
        '[Li+]', '[Zn+2]', '[Fe+2]', '[Fe+3]',  # 其他金属离子
    ]
    
    for pattern in forbidden_patterns:
        if pattern in smiles:
            return False
    
    # 检查是否包含管道符（表示混合物）
    if '|' in smiles:
        return False
    
    # 检查是否包含点号（表示盐或混合物）
    if '.' in smiles:
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # 检查分子大小
        if mol.GetNumAtoms() > 80:  # 降低原子数限制
            return False
        
        # 检查是否包含不常见原子
        atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        allowed_atoms = {'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'H'}
        if not atoms.issubset(allowed_atoms):
            return False
        
        return True
    except:
        return False

def clean_smiles_for_reinvent(input_file, output_file):
    """为REINVENT严格清理SMILES数据"""
    df = pd.read_csv(input_file)
    
    print(f"原始数据: {len(df)} 个化合物")
    
    valid_smiles = []
    incompatible_count = 0
    
    for smiles in df['Smiles']:
        if pd.isna(smiles):
            continue
            
        smiles_str = str(smiles).strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                # 严格检查兼容性
                if is_compatible_with_reinvent(canonical_smiles):
                    valid_smiles.append(canonical_smiles)
                else:
                    incompatible_count += 1
                    print(f"跳过不兼容的SMILES: {smiles_str[:50]}...")
        except:
            incompatible_count += 1
    
    # 去重
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    # 保存
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"\n最终结果:")
    print(f"  兼容的SMILES: {len(unique_smiles)}")
    print(f"  不兼容/跳过: {incompatible_count}")
    print(f"  保存到: {output_file}")

if __name__ == "__main__":
    clean_smiles_for_reinvent("data/NS3.csv", "data/denv_smiles_strict.tsv")
