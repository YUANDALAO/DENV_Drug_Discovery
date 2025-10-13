import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def validate_molecules(file_path):
    df = pd.read_csv(file_path)
    
    valid_mols = []
    invalid_mols = []
    
    for i, row in df.iterrows():
        smiles = row['SMILES']
        
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # 检查Lipinski规则
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            # Lipinski规则检查
            lipinski_pass = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
            
            if lipinski_pass:
                valid_mols.append((smiles, mw, logp, hbd, hba))
            else:
                invalid_mols.append((smiles, mw, logp, hbd, hba, "Lipinski违反"))
        else:
            invalid_mols.append((smiles, 0, 0, 0, 0, "RDKit解析失败"))
    
    print(f"验证结果: {len(valid_mols)}个有效, {len(invalid_mols)}个无效")
    
    if invalid_mols:
        print("\n无效分子:")
        for mol_info in invalid_mols:
            print(f"  {mol_info[0][:50]}... - {mol_info[5]}")
    
    return valid_mols, invalid_mols

validate_molecules("generated_molecules.smi")
