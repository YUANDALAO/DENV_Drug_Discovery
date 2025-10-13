import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def has_scaffold(mol, scaffold_smarts):
    """检查分子是否包含特定骨架"""
    scaffold = Chem.MolFromSmarts(scaffold_smarts)
    if scaffold is None:
        print(f"警告：骨架解析失败 {scaffold_smarts}")
        return False
    return mol.HasSubstructMatch(scaffold)

def filter_molecules_by_scaffolds(input_file, output_file):
    """筛选包含目标骨架的分子"""
    
    # 定义目标骨架的SMARTS模式（更宽松的匹配）
    target_scaffolds = [
        "O=CC1C([*])NC([*])C1",  # 吡咯烷骨架
        "O=CC1C([*])C([*])C1",   # 环丁烷骨架
    ]
    
    df = pd.read_csv(input_file)
    matched_mols = []
    
    for idx, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol:
            for scaffold_smarts in target_scaffolds:
                if has_scaffold(mol, scaffold_smarts):
                    matched_mols.append({
                        'SMILES': row['SMILES'],
                        'NLL': row.get('NLL', 'N/A'),
                        'Matched_Scaffold': scaffold_smarts
                    })
                    break
    
    result_df = pd.DataFrame(matched_mols)
    result_df.to_csv(output_file, index=False)
    
    print(f"筛选结果:")
    print(f"  输入分子: {len(df)}")
    print(f"  匹配目标骨架: {len(matched_mols)}")
    print(f"  匹配率: {len(matched_mols)/len(df)*100:.1f}%")
    print(f"  保存到: {output_file}")

# 筛选已生成的1000个分子
filter_molecules_by_scaffolds(
    'denv_candidates_1000.smi', 
    'scaffold_matched_molecules.smi'
)
