from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import Counter

# 读取高质量分子
hq = df_all[
    (df_all['QED'] > 0.6) &
    (df_all['RO5_Violations'] <= 1)
].copy()

print(f"分析 {len(hq)} 个高质量分子的骨架...")

# 提取骨架
scaffolds = []
for smi in tqdm(hq['SMILES'], desc="提取骨架"):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        try:
            scf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
            scaffolds.append(scf)
        except:
            scaffolds.append(None)

# 统计
scf_counts = Counter([s for s in scaffolds if s])

print(f"\n独特骨架数: {len(scf_counts)}")
print(f"\n产生高质量分子最多的20个骨架:")
for i, (scf, cnt) in enumerate(scf_counts.most_common(20), 1):
    pct = cnt / len(hq) * 100
    print(f"{i:2d}. {scf[:70]} - {cnt:4d} ({pct:.2f}%)")

# 保存
with open("top_scaffolds_for_quality.txt", 'w') as f:
    for scf, cnt in scf_counts.most_common(100):
        f.write(f"{scf}\t{cnt}\n")

print("\n✓ Top 100骨架已保存: top_scaffolds_for_quality.txt")