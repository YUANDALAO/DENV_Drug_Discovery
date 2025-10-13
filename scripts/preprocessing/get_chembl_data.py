import pandas as pd
import requests
from rdkit import Chem
import time

def download_antiviral_compounds():
    """从ChEMBL下载抗病毒化合物数据"""
    
    # 使用ChEMBL FTP服务下载大数据集
    print("从ChEMBL数据库获取抗病毒化合物...")
    
    # 这里我们使用一个预处理的方法
    # 实际项目中可以从ChEMBL FTP下载完整数据
    
    # 或者从其他开源数据集获取
    urls = [
        "https://raw.githubusercontent.com/molecularsets/moses/master/data/dataset_v1.csv",
        # 添加其他可靠的数据源
    ]
    
    all_smiles = []
    
    for url in urls:
        try:
            print(f"下载: {url}")
            df = pd.read_csv(url, nrows=50000)  # 限制大小
            if 'SMILES' in df.columns:
                smiles_col = 'SMILES'
            elif 'smiles' in df.columns:
                smiles_col = 'smiles'
            else:
                smiles_col = df.columns[0]  # 假设第一列是SMILES
            
            valid_smiles = []
            for smiles in df[smiles_col]:
                if pd.notna(smiles):
                    mol = Chem.MolFromSmiles(str(smiles))
                    if mol and mol.GetNumAtoms() < 80:
                        valid_smiles.append(Chem.MolToSmiles(mol))
            
            all_smiles.extend(valid_smiles)
            print(f"获得 {len(valid_smiles)} 个有效分子")
            
        except Exception as e:
            print(f"下载失败: {e}")
            continue
    
    # 去重并保存
    unique_smiles = list(set(all_smiles))
    
    with open('data/large_dataset.tsv', 'w') as f:
        for smiles in unique_smiles[:20000]:  # 限制为20K分子
            f.write(f"{smiles}\n")
    
    print(f"保存了 {min(len(unique_smiles), 20000)} 个独特分子到 data/large_dataset.tsv")
    return len(unique_smiles)

# 运行下载
count = download_antiviral_compounds()
