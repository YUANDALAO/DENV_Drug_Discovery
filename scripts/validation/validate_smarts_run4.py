from rdkit import Chem

# 修正的SMARTS
PYRR = "[c,n,o,s]N1[C,c][C,c][C,c]([#6](=[O])[#7])[C,c]1"
CYCLO = "[c,n,o,s][C,c]1[C,c][C,c]([#6](=[O])[#7])[C,c]1"

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

pyrr_pat = Chem.MolFromSmarts(PYRR)
cyclo_pat = Chem.MolFromSmarts(CYCLO)

if not pyrr_pat:
    print("吡咯烷SMARTS无效！")
    exit(1)
if not cyclo_pat:
    print("环丁烷SMARTS无效！")
    exit(1)

print("="*80)
print("Run4 SMARTS验证")
print("="*80)
print(f"\n吡咯烷: {PYRR}")
print(f"环丁烷: {CYCLO}")

all_pass = True
for expected, smi in seeds:
    mol = Chem.MolFromSmiles(smi)
    p = mol.HasSubstructMatch(pyrr_pat)
    c = mol.HasSubstructMatch(cyclo_pat)
    
    correct = (expected == "Pyrr" and p) or (expected == "Cyclo" and c)
    status = "✓" if correct else "✗"
    
    if not correct:
        all_pass = False
    
    print(f"{status} {expected:5s} | Pyrr:{p} Cyclo:{c} | {smi}")

print("\n" + "="*80)
if all_pass:
    print("✓ 所有种子SMILES验证通过！可以启动Run4")
else:
    print("✗ 验证失败，需要继续调整SMARTS")
print("="*80)
