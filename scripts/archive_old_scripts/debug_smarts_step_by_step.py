from rdkit import Chem

# 测试种子
pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
cyclo_seed = "CC1CC(C1c1ccccc1)C(N)=O"

print("="*80)
print("逐步构建SMARTS")
print("="*80)

# 分析吡咯烷
mol = Chem.MolFromSmiles(pyrr_seed)
print(f"\n吡咯烷: {pyrr_seed}")
print("标注原子索引:")
for atom in mol.GetAtoms():
    print(f"  Idx {atom.GetIdx()}: {atom.GetSymbol()} (in ring: {atom.IsInRing()})")

# 测试各种SMARTS变体
smarts_tests = [
    ("简单5元含N环", "[nR1]1[CR1][CR1][CR1][CR1]1"),
    ("5元环+酰胺", "[nR1]1[CR1][CR1][CR1]([C](=O)[N])[CR1]1"),
    ("N连芳环", "[c][nR1]1[CR1][CR1][CR1][CR1]1"),
    ("完整尝试1", "[c][nR1]1[CR1][CR1][CR1]([C](=O)[N])[CR1]1"),
    ("完整尝试2", "c-[nR1]1-[CR1]-[CR1]-[CR1](-[C](=O)-[N])-[CR1]-1"),
    ("宽松版", "[c,n,o,s][N]1[C][C][C]([C](=O)[N])[C]1"),
]

print("\n吡咯烷SMARTS测试:")
for desc, smarts in smarts_tests:
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol.HasSubstructMatch(pat)
        print(f"  {match} | {desc:20s} | {smarts}")
    else:
        print(f"  ✗ | {desc:20s} | INVALID")

# 分析环丁烷
mol2 = Chem.MolFromSmiles(cyclo_seed)
print(f"\n环丁烷: {cyclo_seed}")

smarts_tests2 = [
    ("简单4元环", "[CR1]1[CR1][CR1][CR1]1"),
    ("4元环+酰胺", "[CR1]1[CR1][CR1]([C](=O)[N])[CR1]1"),
    ("C连芳环", "[c][CR1]1[CR1][CR1][CR1]1"),
    ("完整尝试1", "[c][CR1]1[CR1][CR1]([C](=O)[N])[CR1]1"),
    ("宽松版", "[c,n,o,s][C]1[C][C]([C](=O)[N])[C]1"),
]

print("\n环丁烷SMARTS测试:")
for desc, smarts in smarts_tests2:
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol2.HasSubstructMatch(pat)
        print(f"  {match} | {desc:20s} | {smarts}")
    else:
        print(f"  ✗ | {desc:20s} | INVALID")
