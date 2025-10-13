from rdkit import Chem

# 最简单直接的方法：只匹配核心子结构，不管环闭合
print("最终方案：匹配关键子结构模式\n")

# 吡咯烷关键特征: C(N环)芳环 + 酰胺
# 环丁烷关键特征: C(环)芳环 + 酰胺

smarts_final = [
    ("Pyrr-A", "C(N)a"),  # C连N和芳香原子
    ("Pyrr-B", "C(=O)N"),  # 酰胺
    ("Cyclo-A", "[C;R]a"),  # 环C连芳香
    ("Cyclo-B", "C(=O)N"),  # 酰胺
]

seeds = {
    "Pyrr": "CC1CC(C(N1)c1ccccc1)C(N)=O",
    "Cyclo": "CC1CC(C1c1ccccc1)C(N)=O",
}

print("单个模式测试:")
for name, smarts in smarts_final:
    typ = "Pyrr" if "Pyrr" in name else "Cyclo"
    mol = Chem.MolFromSmiles(seeds[typ])
    pat = Chem.MolFromSmarts(smarts)
    match = mol.HasSubstructMatch(pat) if pat else False
    print(f"  {'✓' if match else '✗'} {name:10s} {smarts:15s}")

# 组合：用逗号分隔的SMARTS表示"包含这个或那个"
# 但我们要"同时包含"，需要用递归SMARTS
print("\n使用递归SMARTS $(...):")

PYRR_RECURSIVE = "[C;$(C(N)a);$(C~C~C(=O)N)],[C;$(C(=O)N);$(C~C(N)a)]"
CYCLO_RECURSIVE = "[C;R;$(Ca);$(C~C(=O)N)],[C;$(C(=O)N);$(C~Ca)]"

print(f"吡咯烷: {PYRR_RECURSIVE}")
print(f"环丁烷: {CYCLO_RECURSIVE}")

# 最简单的方法：手动检查两个子结构
print("\n最简方法：分别检查两个子结构是否都匹配\n")

def check_scaffold(smiles, scaffold_type):
    mol = Chem.MolFromSmiles(smiles)
    
    if scaffold_type == "Pyrr":
        pat1 = Chem.MolFromSmarts("C(N)a")  # C连N和芳香
        pat2 = Chem.MolFromSmarts("C(=O)N")  # 酰胺
    else:
        pat1 = Chem.MolFromSmarts("[C;R]a")  # 环C连芳香  
        pat2 = Chem.MolFromSmarts("C(=O)N")  # 酰胺
    
    match1 = mol.HasSubstructMatch(pat1)
    match2 = mol.HasSubstructMatch(pat2)
    
    return match1 and match2

all_seeds = [
    ("Pyrr", "CC1CC(C(N1)c1ccccc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(F)cc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(Cl)cc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(OC)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccccc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(F)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(Cl)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(OC)cc1)C(N)=O"),
]

print("完整验证（手动双匹配法）:")
for typ, smi in all_seeds:
    result = check_scaffold(smi, typ)
    print(f"  {'✓' if result else '✗'} {typ:5s} | {smi}")

print("\n" + "="*80)
print("结论：REINVENT配置中使用CustomAlerts或两个独立的MatchingSubstructure")
print("  吡咯烷: 同时匹配 C(N)a 和 C(=O)N")
print("  环丁烷: 同时匹配 [C;R]a 和 C(=O)N")
