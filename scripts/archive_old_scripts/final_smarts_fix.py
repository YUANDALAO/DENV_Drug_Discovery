from rdkit import Chem

pyrr_seed = "CC1CC(C(N1)c1ccccc1)C(N)=O"
cyclo_seed = "CC1CC(C1c1ccccc1)C(N)=O"

# 正确的SMARTS - 考虑双环系统
smarts_candidates = [
    # 吡咯烷 - 使用[nR]不限定环数
    ("Pyrr-v1", "[c][nR]1[CR][CR][CR]([C](=O)[N])[CR]1"),
    ("Pyrr-v2", "[c,n,o,s][NR]1[CR][CR][CR]([C](=O)[N])[CR]1"),
    ("Pyrr-v3", "[c,n,o,s]~[NR]~1~[CR]~[CR]~[CR](~[C](=[O])~[N])~[CR]~1"),
    
    # 环丁烷
    ("Cyclo-v1", "[c][CR]1[CR][CR]([C](=O)[N])[CR]1"),
    ("Cyclo-v2", "[c,n,o,s][CR]1[CR][CR]([C](=O)[N])[CR]1"),
    ("Cyclo-v3", "[c,n,o,s]~[CR]~1~[CR]~[CR](~[C](=[O])~[N])~[CR]~1"),
]

mol_p = Chem.MolFromSmiles(pyrr_seed)
mol_c = Chem.MolFromSmiles(cyclo_seed)

print("测试吡咯烷 SMARTS:")
for name, smarts in smarts_candidates[:3]:
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol_p.HasSubstructMatch(pat)
        print(f"  {'✓' if match else '✗'} {name:10s} | {smarts}")

print("\n测试环丁烷 SMARTS:")
for name, smarts in smarts_candidates[3:]:
    pat = Chem.MolFromSmarts(smarts)
    if pat:
        match = mol_c.HasSubstructMatch(pat)
        print(f"  {'✓' if match else '✗'} {name:10s} | {smarts}")

print("\n全种子验证:")
seeds = [
    ("Pyrr", "CC1CC(C(N1)c1ccccc1)C(N)=O"),
    ("Pyrr", "CC1CC(C(N1)c1ccc(F)cc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccccc1)C(N)=O"),
    ("Cyclo", "CC1CC(C1c1ccc(F)cc1)C(N)=O"),
]

# 选择最佳SMARTS
PYRR_BEST = "[c,n,o,s]~[NR]~1~[CR]~[CR]~[CR](~[C](=[O])~[N])~[CR]~1"
CYCLO_BEST = "[c,n,o,s]~[CR]~1~[CR]~[CR](~[C](=[O])~[N])~[CR]~1"

pyrr_pat = Chem.MolFromSmarts(PYRR_BEST)
cyclo_pat = Chem.MolFromSmarts(CYCLO_BEST)

for typ, smi in seeds:
    mol = Chem.MolFromSmiles(smi)
    p = mol.HasSubstructMatch(pyrr_pat)
    c = mol.HasSubstructMatch(cyclo_pat)
    ok = (typ=="Pyrr" and p) or (typ=="Cyclo" and c)
    print(f"  {'✓' if ok else '✗'} {typ:5s} | P:{p} C:{c}")

if all((typ=="Pyrr" and Chem.MolFromSmiles(smi).HasSubstructMatch(pyrr_pat)) or 
       (typ=="Cyclo" and Chem.MolFromSmiles(smi).HasSubstructMatch(cyclo_pat)) 
       for typ, smi in seeds):
    print("\n最终SMARTS:")
    print(f"  吡咯烷: {PYRR_BEST}")
    print(f"  环丁烷: {CYCLO_BEST}")
