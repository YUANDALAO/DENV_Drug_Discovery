from rdkit import Chem

pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
mol = Chem.MolFromSmiles(pyrr_seed)

print("吡咯烷种子结构详细分析:")
print(f"SMILES: {pyrr_seed}\n")

# 找到N原子和它的邻居
for atom in mol.GetAtoms():
    if atom.GetSymbol() == 'N' and atom.IsInRing():
        print(f"环N原子 (Idx {atom.GetIdx()}):")
        for neighbor in atom.GetNeighbors():
            print(f"  邻居 {neighbor.GetIdx()}: {neighbor.GetSymbol()} aromatic={neighbor.GetIsAromatic()} inRing={neighbor.IsInRing()}")

# 测试不同的芳环-N连接SMARTS
tests = [
    ("芳碳-N", "[c][N;R]"),
    ("任意芳-N", "[a][N;R]"),  # a = 任何芳香原子
    ("递归芳环", "[$(c1ccccc1)][N;R]"),
    ("6元芳环", "[c;r6][N;R]"),
]

print("\n测试芳环-N键的SMARTS:")
for desc, smarts in tests:
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol.HasSubstructMatch(pat)
        print(f"  {'✓' if match else '✗'} {desc:15s} | {smarts}")

# 最终SMARTS组合测试
print("\n完整SMARTS组合测试:")

final_smarts = [
    ("Pyrr-v1", "[a][N;R]", "[C;R][C](=O)[N]"),
    ("Pyrr-v2", "[c;r6][N;R]", "[C;R][C](=O)[N]"),
    ("Cyclo-v1", "[a][C;R]", "[C;R][C](=O)[N]"),
]

for name, part1, part2 in final_smarts:
    smarts = f"{part1}.{part2}"
    pat = Chem.MolFromSmarts(smarts)
    test_mol = Chem.MolFromSmiles(pyrr_seed if "Pyrr" in name else "CC1CC(C1c1ccccc1)C(N)=O")
    if pat:
        match = test_mol.HasSubstructMatch(pat)
        print(f"  {'✓' if match else '✗'} {name:10s} | {smarts}")

print("\n推荐最终方案:")
print("  吡咯烷: [a][N;R].[C;R][C](=O)[N]")
print("  环丁烷: [a][C;R].[C;R][C](=O)[N]")
print("  说明: [a] = 任何芳香原子, . = 两个片段都要存在")
