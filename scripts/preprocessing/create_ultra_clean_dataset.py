import pandas as pd
from rdkit import Chem

def is_ultra_compatible(smiles):
    """超严格兼容性检查"""
    # 明确禁止的token
    forbidden_tokens = ['[n-]', '[s+]', '[N-]', '[S-]', '[O+]', '[C-]', '[c-]']
    
    for token in forbidden_tokens:
        if token in smiles:
            return False
    
    # 避免复杂的芳香环标记
    if any(x in smiles for x in ['[nH+]', '[s+]', '[o+]']):
        return False
    
    # 避免立体化学标记
    if any(x in smiles for x in ['@', '/', '\\']):
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
            
        # 严格的分子大小限制
        if mol.GetNumAtoms() > 50:
            return False
            
        # 避免不常见原子
        atoms = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        allowed_atoms = {'C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'H'}
        if not atoms.issubset(allowed_atoms):
            return False
            
        return True
    except:
        return False

def create_ultra_clean_dataset(input_file, output_file):
    df = pd.read_csv(input_file)
    
    valid_smiles = []
    
    for smiles in df['Smiles']:
        if pd.isna(smiles):
            continue
            
        smiles_str = str(smiles).strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                # 基本过滤
                if ('.' not in canonical_smiles and 
                    '|' not in canonical_smiles and
                    len(canonical_smiles) < 100 and
                    is_ultra_compatible(canonical_smiles)):
                    valid_smiles.append(canonical_smiles)
        except:
            continue
    
    # 去重并保存
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"超清洁数据集: {len(unique_smiles)} 个SMILES")
    return len(unique_smiles)

# 创建超清洁数据集
count = create_ultra_clean_dataset("data/NS3.csv", "data/denv_ultra_clean.tsv")
