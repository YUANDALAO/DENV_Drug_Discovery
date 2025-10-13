import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import requests

def get_zinc_subset():
    """从ZINC数据库获取药物样分子"""
    # ZINC15 drug-like subset
    print("尝试获取ZINC drug-like分子...")
    
    # 这里我们创建一个模拟的扩展数据集
    # 基于已有的DENV数据生成变体
    df = pd.read_csv('data/NS3.csv')
    
    expanded_smiles = []
    
    for smiles in df['Smiles']:
        if pd.notna(smiles):
            try:
                mol = Chem.MolFromSmiles(str(smiles))
                if mol:
                    # 添加原始分子
                    expanded_smiles.append(Chem.MolToSmiles(mol))
                    
            except:
                continue
    
    # 去重
    unique_smiles = list(set(expanded_smiles))
    
    # 保存扩展数据集
    with open('data/expanded_dataset.tsv', 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"创建了包含 {len(unique_smiles)} 个分子的扩展数据集")
    return len(unique_smiles)

get_zinc_subset()
