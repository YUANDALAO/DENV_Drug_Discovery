from rdkit import Chem
from rdkit.Chem import Descriptors
import itertools

def generate_scaffold_library():
    """在指定骨架上系统性地生成R基团组合"""
    
    scaffolds = [
        ("O=C([*])C1C([*])NC([*])C1", "Pyrrolidine"),
        ("O=C([*])C1C([*])C([*])C1", "Cyclobutane"),
    ]
    
    # 精选的药物样R基团
    r_groups = {
        'small': ["C", "CC", "CCC"],
        'aromatic': ["c1ccccc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1", "c1ccc(C)cc1"],
        'hetero': ["c1ccncc1", "c1ccoc1", "c1csccc1"],
        'polar': ["O", "OC", "N", "NC", "C(=O)C"],
        'halogen': ["F", "Cl", "Br"],
    }
    
    all_r = (r_groups['small'] + r_groups['aromatic'][:3] + 
             r_groups['hetero'][:2] + r_groups['polar'][:3] + 
             r_groups['halogen'])
    
    molecules = []
    
    for scaffold_smarts, name in scaffolds:
        count = 0
        for r1, r2, r3 in itertools.product(all_r[:7], all_r[:6], all_r[:5]):
            try:
                smiles = scaffold_smarts.replace("[*]", f"({r1})", 1)
                smiles = smiles.replace("[*]", f"({r2})", 1)
                smiles = smiles.replace("[*]", f"({r3})", 1)
                
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    canonical = Chem.MolToSmiles(mol)
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    
                    molecules.append({
                        'SMILES': canonical,
                        'Scaffold': name,
                        'MW': mw,
                        'LogP': logp
                    })
                    count += 1
                    
                    if count >= 250:  # 每个骨架250个
                        break
            except:
                continue
        
        if count >= 250:
            break
    
    # 保存
    import pandas as pd
    df = pd.DataFrame(molecules)
    df.to_csv('scaffold_focused_library.csv', index=False)
    
    print(f"生成的骨架定向分子库:")
    print(f"  总分子数: {len(df)}")
    print(f"  平均MW: {df['MW'].mean():.1f}")
    print(f"  平均LogP: {df['LogP'].mean():.1f}")
    
    lipinski = df[(df['MW'] <= 500) & (df['LogP'] <= 5)]
    print(f"  Lipinski通过: {len(lipinski)} ({len(lipinski)/len(df)*100:.1f}%)")
    
    print(f"\n前5个分子示例:")
    for i in range(5):
        print(f"  {df.iloc[i]['SMILES'][:65]}...")

generate_scaffold_library()
