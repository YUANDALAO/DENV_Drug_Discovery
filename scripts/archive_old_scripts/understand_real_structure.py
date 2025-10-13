from rdkit import Chem
from rdkit.Chem import Draw

pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
cyclo_seed = "CC1CC(C1c1ccccc1)C(N)=O"

mol_p = Chem.MolFromSmiles(pyrr_seed)
mol_c = Chem.MolFromSmiles(cyclo_seed)

print("真实结构分析:\n")
print("吡咯烷: CC1CC(C(N1)c1ccccc1)C(N)=O")
print("  解读: C-C1-C-C(-C(-N1-苯基)-酰胺)-C-1")
print("  关键: 芳环连在环碳C4上，C4连N，形成 C(N环)芳环 结构")

print("\n环丁烷: CC1CC(C1c1ccccc1)C(N)=O")  
print("  解读: C-C1-C-C(-C1-苯基)-酰胺")
print("  关键: 芳环直接连环碳C1")

# 正确的SMARTS应该是
print("\n正确的SMARTS策略:")
print("\n吡咯烷型 - 必须有:")
print("  1. 5元环含N")
print("  2. 环上某C连酰胺")
print("  3. 环上某C同时连N和芳环 [C;R]([N;R])[a]")

print("\n环丁烷型 - 必须有:")
print("  1. 4元全C环")
print("  2. 环上某C连酰胺")
print("  3. 环上某C连芳环")

# 测试新SMARTS
new_smarts = [
    ("Pyrr新", "[N;R].[C;R]([N;R])[a].[C;R][C](=O)[N]"),
    ("Cyclo新", "[C;R]1[C;R][C;R][C;R]1.[a][C;R].[C;R][C](=O)[N]"),
]

print("\n测试新SMARTS:")
for name, smarts in new_smarts:
    pat = Chem.MolFromSmarts(smarts)
    if "Pyrr" in name:
        match = mol_p.HasSubstructMatch(pat) if pat else False
        print(f"  {'✓' if match else '✗'} {name:10s} | {smarts}")
    else:
        match = mol_c.HasSubstructMatch(pat) if pat else False
        print(f"  {'✓' if match else '✗'} {name:10s} | {smarts}")

# 全验证
print("\n完整8种子验证:")
seeds = [
    ("Pyrr", "CC1CC(C(N1)c1ccccc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(F)cc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(Cl)cc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(OC)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccccc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(F)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(Cl)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(OC)cc1)C(N)=O"),
]

PYRR = "[N;R].[C;R]([N;R])[a].[C;R][C](=O)[N]"
CYCLO = "[a][C;R].[C;R][C](=O)[N]"

pyrr_pat = Chem.MolFromSmarts(PYRR)
cyclo_pat = Chem.MolFromSmarts(CYCLO)

for typ, smi in seeds:
    mol = Chem.MolFromSmiles(smi)
    p = mol.HasSubstructMatch(pyrr_pat)
    c = mol.HasSubstructMatch(cyclo_pat)
    ok = (typ=="Pyrr" and p) or (typ=="Cyclo" and c)
    print(f"  {'✓' if ok else '✗'} {typ:5s} | P:{p} C:{c}")
