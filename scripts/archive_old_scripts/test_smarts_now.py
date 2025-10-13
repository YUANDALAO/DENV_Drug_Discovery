from rdkit import Chem

# 当前使用的SMARTS
PYRR = "[#6]~[#6]1~[#6]~[#6](-[#6](=[#8])[#7])~[#7]~[#6]~1[c,n,o,s]"
CYCLO = "[#6]~[#6]1~[#6]~[#6](-[#6](=[#8])[#7])~[#6]~1[c,n,o,s]"

# 种子SMILES
seeds = [
    "CC1CC(C(N1)c1ccccc1)C(N)=O",  # 吡咯烷1
    "CC1CC(C(N1)c1ccc(F)cc1)C(N)=O",  # 吡咯烷2
    "CC1CC(C(N1)c1ccc(Cl)cc1)C(N)=O",  # 吡咯烷3
    "CC1CC(C(N1)c1ccc(OC)cc1)C(N)=O",  # 吡咯烷4
    "CC1CC(C1c1ccccc1)C(N)=O",  # 环丁烷1
    "CC1CC(C1c1ccc(F)cc1)C(N)=O",  # 环丁烷2
    "CC1CC(C1c1ccc(Cl)cc1)C(N)=O",  # 环丁烷3
    "CC1CC(C1c1ccc(OC)cc1)C(N)=O",  # 环丁烷4
]

# Top1生成的分子
generated = "COc1cc2cc(C(=O)N3CCN(C)CC3)n(C(C)C)c2cc1O"

pyrr_pat = Chem.MolFromSmarts(PYRR)
cyclo_pat = Chem.MolFromSmarts(CYCLO)

print("="*80)
print("SMARTS验证测试")
print("="*80)

print("\n种子SMILES匹配情况:")
for i, smi in enumerate(seeds, 1):
    mol = Chem.MolFromSmiles(smi)
    p = mol.HasSubstructMatch(pyrr_pat)
    c = mol.HasSubstructMatch(cyclo_pat)
    expected = "Pyrr" if i <= 4 else "Cyclo"
    status = "✓" if (p and i<=4) or (c and i>4) else "✗"
    print(f"{status} {i}. {expected:5s} | Pyrr:{p} Cyclo:{c} | {smi}")

print("\nTop1生成分子:")
mol = Chem.MolFromSmiles(generated)
p = mol.HasSubstructMatch(pyrr_pat)
c = mol.HasSubstructMatch(cyclo_pat)
print(f"  Pyrr:{p} Cyclo:{c}")
print(f"  {generated}")

# 分析匹配位置
if p:
    matches = mol.GetSubstructMatches(pyrr_pat)
    print(f"\n  吡咯烷SMARTS匹配了 {len(matches)} 个位置")
if c:
    matches = mol.GetSubstructMatches(cyclo_pat)
    print(f"  环丁烷SMARTS匹配了 {len(matches)} 个位置")
