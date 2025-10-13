from rdkit import Chem

# 仔细分析SMILES
print("SMILES结构拆解:\n")
print("吡咯烷: CC1CC(C(N1)c1ccccc1)C(N)=O")
print("  C-C1-C-C( C(N1)(苯基) )(酰胺) ")
print("  ^^^^环^^  ^^^^^中心C^^^^^")
print("  中心C有3个邻居: 环C, 含N的碳C(N1)(苯), 酰胺")
print("  含N的碳C: 连N和苯环")

print("\n正确理解:")
print("  5元环: C-C-C-C-N")
print("  环外有个C，这个C连接: N(环) + 苯环")
print("  环上某C连: 酰胺 + 那个外部C")

# 简化策略：只要求关键片段存在即可
simple_smarts = [
    ("Pyrr极简", "[N;R]", "[C]([N])[a]", "[C](=O)[N]"),  # 分开3个片段
    ("Cyclo极简", "[C;R]1[C;R][C;R][C;R]1", "[a][C]", "[C](=O)[N]"),
]

print("\n测试极简SMARTS（3个独立片段）:")

seeds = {
    "Pyrr": ["CC1CC(C(N1)c1ccccc1)C(N)=O", "CC1CC(C(N1)c1ccc(F)cc1)C(N)=O"],
    "Cyclo": ["CC1CC(C1c1ccccc1)C(N)=O", "CC1CC(C1c1ccc(F)cc1)C(N)=O"],
}

# 吡咯烷: 要求环N + C(N)(芳) + 酰胺都存在
PYRR = "[N;R].[C]([N])[a].[C](=O)[N]"
CYCLO = "[C;R]1[C;R][C;R][C;R]1.[a][C].[C](=O)[N]"

pyrr_pat = Chem.MolFromSmarts(PYRR)
cyclo_pat = Chem.MolFromSmarts(CYCLO)

for typ in ["Pyrr", "Cyclo"]:
    print(f"\n{typ}测试:")
    for smi in seeds[typ]:
        mol = Chem.MolFromSmiles(smi)
        p = mol.HasSubstructMatch(pyrr_pat)
        c = mol.HasSubstructMatch(cyclo_pat)
        ok = (typ=="Pyrr" and p) or (typ=="Cyclo" and c)
        print(f"  {'✓' if ok else '✗'} P:{p} C:{c} | {smi[:30]}")

print("\n" + "="*80)
print("最终SMARTS (如果全✓):")
print(f"  吡咯烷: {PYRR}")
print(f"  环丁烷: {CYCLO}")
