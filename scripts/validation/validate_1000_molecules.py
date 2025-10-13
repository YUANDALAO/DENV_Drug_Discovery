import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import numpy as np

def validate_large_set(file_path):
    df = pd.read_csv(file_path)
    
    valid_count = 0
    lipinski_pass = 0
    mw_values = []
    logp_values = []
    
    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol:
            valid_count += 1
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            
            mw_values.append(mw)
            logp_values.append(logp)
            
            if mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10:
                lipinski_pass += 1
    
    print(f"生成统计:")
    print(f"  总生成: {len(df)} 个分子")
    print(f"  有效分子: {valid_count} ({valid_count/len(df)*100:.1f}%)")
    print(f"  Lipinski通过: {lipinski_pass} ({lipinski_pass/len(df)*100:.1f}%)")
    print(f"  分子量范围: {np.mean(mw_values):.1f}±{np.std(mw_values):.1f}")
    print(f"  LogP范围: {np.mean(logp_values):.1f}±{np.std(logp_values):.1f}")

validate_large_set("denv_candidates_1000.smi")
