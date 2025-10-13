from rdkit import Chem
import pandas as pd

def extract_rgroups_from_molecules(smiles_file):
    """从现有分子中提取常见R基团"""
    
    # 常见药物R基团
    common_rgroups = [
        # 简单烷基
        "C", "CC", "CCC", "C(C)C",
        # 芳香基团
        "c1ccccc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1", "c1ccc(Br)cc1",
        "c1ccc(OC)cc1", "c1ccc(C)cc1",
        # 杂环
        "c1ccncc1", "c1cnccc1", "c1ccoc1", "c1ccsc1",
        # 含氧基团
        "OC", "OCC", "C(=O)C", "C(=O)OC",
        # 含氮基团
        "NC", "NCC", "N(C)C", "NC(=O)C",
        # 卤素
        "F", "Cl", "Br",
        # 其他
        "CF3", "CN", "C#N", "S(=O)(=O)C",
    ]
    
    # 保存R基团库
    with open('data/decorations.smi', 'w') as f:
        for rgroup in common_rgroups:
            f.write(f"{rgroup}\n")
    
    print(f"创建了 {len(common_rgroups)} 个R基团")
    return len(common_rgroups)

extract_rgroups_from_molecules('data/denv_ultra_clean.tsv')
