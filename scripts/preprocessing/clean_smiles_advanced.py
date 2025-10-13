import pandas as pd
from rdkit import Chem

def is_compatible_with_reinvent(smiles):
    """检查SMILES是否与REINVENT模型兼容"""
    # REINVENT支持的token集合
    allowed_tokens = {
        '[N+]', '[S+]', '#', 'C', '[N-]', 'c', '-', '[nH]', '3', 'Br', '1', ')', 
        'n', 'S', '9', '$', 'Cl', '=', '6', '2', '%10', 'N', 's', '7', '4', 'O', 
        '5', '8', '^', '[O-]', 'o', 'F', '(', '[n+]'
    }
    
    # 检查是否包含不支持的token
    forbidden_patterns = ['[Cl-]', '[Br-]', '[I-]', '[F-]']
    for pattern in forbidden_patterns:
        if pattern in smiles:
            return False
    
    # 检查分子复杂度
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # 避免过于复杂的分子
        if mol.GetNumAtoms() > 100:  # 原子数限制
            return False
            
        return True
    except:
        return False

def clean_smiles_for_reinvent(input_file, output_file):
    """为REINVENT清理SMILES数据"""
    df = pd.read_csv(input_file)
    
    print(f"原始数据: {len(df)} 个化合物")
    
    valid_smiles = []
    incompatible_count = 0
    
    for smiles in df['Smiles']:
        if pd.isna(smiles):
            continue
            
        smiles_str = str(smiles)
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                # 检查REINVENT兼容性
                if is_compatible_with_reinvent(canonical_smiles):
                    valid_smiles.append(canonical_smiles)
                else:
                    incompatible_count += 1
        except:
            incompatible_count += 1
    
    # 去重
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    # 保存
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"兼容的SMILES: {len(unique_smiles)}")
    print(f"不兼容/跳过: {incompatible_count}")
    print(f"保存到: {output_file}")

if __name__ == "__main__":
    clean_smiles_for_reinvent("data/NS3.csv", "data/denv_smiles_compatible.tsv")
