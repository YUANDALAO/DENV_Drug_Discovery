from rdkit import Chem

pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
cyclo_seed = "CC1CC(C1c1ccccc1)C(N)=O"

print("使用递归SMARTS和子结构组合方法:\n")

# 方法1：分离匹配各部分
part_tests = [
    ("5元含N环", "[N;R]1~[C;R]~[C;R]~[C;R]~[C;R]~1", pyrr_seed),
    ("4元C环", "[C;R]1~[C;R]~[C;R]~[C;R]~1", cyclo_seed),
    ("酰胺基", "[C](=[O])[N;H0,H1,H2]", pyrr_seed),
    ("芳环-N键", "[c,n,o,s]~[N;R]", pyrr_seed),
    ("芳环-C键", "[c,n,o,s]~[C;R]", cyclo_seed),
]

for desc, smarts, test_smi in part_tests:
    mol = Chem.MolFromSmiles(test_smi)
    pat = Chem.MolFromSmarts(smarts)
    match = mol.HasSubstructMatch(pat) if pat else False
    print(f"  {'✓' if match else '✗'} {desc:15s} | {smarts}")

# 方法2：使用递归SMARTS $(...)
print("\n递归SMARTS测试:")
recursive_tests = [
    ("Pyrr递归", "[N;R;$(N([c,n,o,s])[C;R][C;R][C;R]([C](=O)[N])[C;R])]", pyrr_seed),
    ("Cyclo递归", "[C;R;$([C]([c,n,o,s])[C;R][C;R]([C](=O)[N])[C;R])]", cyclo_seed),
]

for desc, smarts, test_smi in recursive_tests:
    mol = Chem.MolFromSmiles(test_smi)
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol.HasSubstructMatch(pat)
        print(f"  {'✓' if match else '✗'} {desc:15s}")
    else:
        print(f"  INVALID {desc}")

# 方法3：最简单 - 只要求核心特征存在
print("\n极简SMARTS（仅要求关键特征共存）:")
simple_tests = [
    ("Pyrr简化", "[c,n,o,s][N;R].[C;R][C](=O)[N]", pyrr_seed),
    ("Cyclo简化", "[c,n,o,s][C;R].[C;R][C](=O)[N]", cyclo_seed),
]

for desc, smarts, test_smi in simple_tests:
    mol = Chem.MolFromSmiles(test_smi)
    pat = Chem.MolFromSmarts(smarts)
    match = mol.HasSubstructMatch(pat) if pat else False
    print(f"  {'✓' if match else '✗'} {desc:15s} | {smarts}")

print("\n建议：用最简SMARTS，要求分子中同时包含:")
print("  吡咯烷: (1)芳环连N  (2)环C连酰胺")
print("  环丁烷: (1)芳环连C  (2)环C连酰胺")
