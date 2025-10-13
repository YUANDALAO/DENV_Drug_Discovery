from rdkit import Chem

# 种子SMILES分析
pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
cyclo_seed = "CC1CC(C1c1ccccc1)C(N)=O"

mol_p = Chem.MolFromSmiles(pyrr_seed)
mol_c = Chem.MolFromSmiles(cyclo_seed)

print("吡咯烷种子结构分析:")
print(f"  SMILES: {pyrr_seed}")
print(f"  原子数: {mol_p.GetNumAtoms()}")
print("  结构: 5元环(C-C-C-C-N), N上连苯基, C3上连酰胺, C2上连甲基")

print("\n环丁烷种子结构分析:")
print(f"  SMILES: {cyclo_seed}")
print(f"  原子数: {mol_c.GetNumAtoms()}")
print("  结构: 4元环(C-C-C-C), C1上连苯基, C2上连酰胺, C4上连甲基")

# 正确的SMARTS (更宽松的匹配)
print("\n修正方案:")
print("\n吡咯烷型 - 核心特征:")
print("  1. 5元环，包含1个N")
print("  2. 环上某个C连接酰胺 C(=O)N")
print("  3. 环的N连接芳环")
print("  建议SMARTS: [c,n,o,s]N1[C,c][C,c][C,c]([#6](=[O])[#7])[C,c]1")

print("\n环丁烷型 - 核心特征:")
print("  1. 4元全C环")
print("  2. 环上某个C连接酰胺 C(=O)N")
print("  3. 环上某个C连接芳环")
print("  建议SMARTS: [c,n,o,s][C,c]1[C,c][C,c]([#6](=[O])[#7])[C,c]1")
