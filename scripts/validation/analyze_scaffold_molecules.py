import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def analyze_generated():
    df = pd.read_csv('scaffold_decorated_molecules.smi')
    
    # 检查是否包含目标骨架
    target_scaffolds = [
        "O=C([*])C1C([*])NC([*])C1",
        "O=C([*])C1C([*])C([*])C1"
    ]
    
    scaffold_counts = {s: 0 for s in target_scaffolds}
    valid_mols = []
    
    for idx, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # 计算性质
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            
            valid_mols.append({
                'SMILES': smiles,
                'MW': mw,
                'LogP': logp
            })
    
    print(f"生成分析:")
    print(f"  总生成: {len(df)}")
    print(f"  有效分子: {len(valid_mols)}")
    print(f"  平均分子量: {sum(m['MW'] for m in valid_mols)/len(valid_mols):.1f}")
    print(f"  平均LogP: {sum(m['LogP'] for m in valid_mols)/len(valid_mols):.1f}")
    
    # 保存前10个示例
    print("\n前10个生成的分子:")
    for i in range(min(10, len(valid_mols))):
        print(f"  {i+1}. {valid_mols[i]['SMILES'][:60]}...")

analyze_generated()
