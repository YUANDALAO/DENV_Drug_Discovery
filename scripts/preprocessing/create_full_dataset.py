import pandas as pd
from rdkit import Chem

def create_full_dataset(input_file, output_file, max_molecules=800):
    """创建适中大小的训练数据集"""
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
                
                # 基本过滤：保留相对简单的分子
                if ('.' not in canonical_smiles and 
                    '|' not in canonical_smiles and
                    len(canonical_smiles) < 120 and
                    mol.GetNumAtoms() < 60 and
                    'B(' not in canonical_smiles):  # 避免硼化合物
                    valid_smiles.append(canonical_smiles)
        except:
            continue
    
    # 去重并限制数量
    unique_smiles = list(dict.fromkeys(valid_smiles))[:max_molecules]
    
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"完整数据集: {len(unique_smiles)} 个SMILES已保存到 {output_file}")

create_full_dataset("data/NS3.csv", "data/denv_training_set.tsv", 800)
