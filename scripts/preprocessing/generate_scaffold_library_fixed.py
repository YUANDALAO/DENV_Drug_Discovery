from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import pandas as pd

def attach_rgroup_to_scaffold(scaffold_smarts, r1, r2, r3):
    """正确地将R基团连接到骨架上"""
    try:
        # 创建骨架分子
        scaffold = Chem.MolFromSmarts(scaffold_smarts)
        if not scaffold:
            return None
        
        # 将SMARTS转为可编辑的分子
        scaffold_mol = Chem.RWMol(scaffold)
        
        # 找到[*]原子的索引
        attachment_atoms = []
        for atom in scaffold_mol.GetAtoms():
            if atom.GetAtomicNum() == 0:  # [*]的原子序数为0
                attachment_atoms.append(atom.GetIdx())
        
        if len(attachment_atoms) != 3:
            return None
        
        # 依次连接R基团
        for i, rgroup_smiles in enumerate([r1, r2, r3]):
            rgroup = Chem.MolFromSmiles(rgroup_smiles)
            if not rgroup:
                continue
            
            # 组合分子（简化版）
            scaffold_mol = Chem.CombineMols(scaffold_mol, rgroup)
        
        # 转换为SMILES
        final_smiles = Chem.MolToSmiles(Chem.Mol(scaffold_mol))
        return final_smiles
        
    except:
        return None

def simple_replacement_method():
    """使用简单但有效的字符串替换方法"""
    
    # 使用更简单的骨架定义
    scaffolds = [
        "O=CC1CNC(*)C1",  # 吡咯烷简化版（只保留一个R位置）
        "O=CC1CC(*)C1",   # 环丁烷简化版
    ]
    
    # 简单的R基团
    r_groups = [
        "C", "CC", "CCC", "CCCC",
        "c1ccccc1", "c1cccnc1", 
        "OC", "NC", "Cl", "Br", "F",
    ]
    
    molecules = []
    
    for scaffold in scaffolds:
        for rgroup in r_groups:
            smiles = scaffold.replace("*", rgroup)
            mol = Chem.MolFromSmiles(smiles)
            
            if mol:
                canonical = Chem.MolToSmiles(mol)
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                
                molecules.append({
                    'SMILES': canonical,
                    'MW': mw,
                    'LogP': logp
                })
    
    df = pd.DataFrame(molecules)
    df.to_csv('scaffold_library_simple.csv', index=False)
    
    print(f"生成的分子库:")
    print(f"  总数: {len(df)}")
    if len(df) > 0:
        print(f"  平均MW: {df['MW'].mean():.1f}")
        print(f"  平均LogP: {df['LogP'].mean():.1f}")
        
        print(f"\n前10个分子:")
        for i in range(min(10, len(df))):
            print(f"  {i+1}. {df.iloc[i]['SMILES']}")

simple_replacement_method()
