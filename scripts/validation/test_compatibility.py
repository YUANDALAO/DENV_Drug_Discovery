import pandas as pd
from rdkit import Chem

def simple_clean(input_file, output_file):
    """简化但有效的清理方法"""
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
                
                # 基本过滤：避免复杂结构
                if ('.' not in canonical_smiles and 
                    '|' not in canonical_smiles and
                    len(canonical_smiles) < 150 and
                    mol.GetNumAtoms() < 70):
                    valid_smiles.append(canonical_smiles)
        except:
            continue
    
    # 去重并限制数量进行测试
    unique_smiles = list(dict.fromkeys(valid_smiles))[:100]  # 先用100个测试
    
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"测试数据集: {len(unique_smiles)} 个SMILES")

simple_clean("data/NS3.csv", "data/denv_test.tsv")
