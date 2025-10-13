import pandas as pd
import glob
from rdkit import Chem

# 找Run4结果文件
files = glob.glob("experiments/runs/run4/results_*.csv")
if not files:
    print("Run4还未生成结果文件，训练可能未完成或刚启动")
    exit(0)

latest = max(files, key=lambda x: x.split('_')[-1])
df = pd.read_csv(latest, low_memory=False)

print("="*80)
print(f"Run4分析 (共{len(df):,}个分子)")
print("="*80)

# 骨架SMARTS
pyrr_a = Chem.MolFromSmarts("C(N)a")
pyrr_b = Chem.MolFromSmarts("C(=O)N")
cyclo_a = Chem.MolFromSmarts("[C;R]a")
cyclo_b = Chem.MolFromSmarts("C(=O)N")

# 真实骨架匹配检查
matched = []
for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        pyrr_match = (mol.HasSubstructMatch(pyrr_a) and mol.HasSubstructMatch(pyrr_b))
        cyclo_match = (mol.HasSubstructMatch(cyclo_a) and mol.HasSubstructMatch(cyclo_b))
        if pyrr_match or cyclo_match:
            matched.append(row)

print(f"\n真实骨架匹配: {len(matched)} / {len(df)} ({len(matched)/len(df)*100:.1f}%)")

if 'Pyrr_Feature_A (raw)' in df.columns:
    print("\nGroupCount统计分布:")
    print(f"  Pyrr_A平均: {df['Pyrr_Feature_A (raw)'].mean():.2f}")
    print(f"  Pyrr_B平均: {df['Pyrr_Feature_B (raw)'].mean():.2f}")
    print(f"  Cyclo_A平均: {df['Cyclo_Feature_A (raw)'].mean():.2f}")
    print(f"  Cyclo_B平均: {df['Cyclo_Feature_B (raw)'].mean():.2f}")

print(f"\nScore统计:")
print(f"  平均: {df['Score'].mean():.4f}")
print(f"  最大: {df['Score'].max():.4f}")

if 'DENV_Activity_pIC50 (raw)' in df.columns:
    print(f"\npIC50统计:")
    print(f"  平均: {df['DENV_Activity_pIC50 (raw)'].mean():.2f}")
    print(f"  最大: {df['DENV_Activity_pIC50 (raw)'].max():.2f}")
    print(f"  >8.0: {(df['DENV_Activity_pIC50 (raw)'] > 8.0).sum()}")

if matched:
    print("\n匹配骨架的Top10:")
    matched_df = pd.DataFrame(matched)
    top10 = matched_df.nlargest(10, 'Score')[['SMILES', 'Score']].head(10)
    print(top10.to_string(index=False))
