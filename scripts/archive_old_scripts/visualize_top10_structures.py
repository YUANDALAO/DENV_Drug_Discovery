#!/usr/bin/env python3
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import sys
import os

if len(sys.argv) > 1:
    run_folder = sys.argv[1]
else:
    run_folder = "experiments/runs/run3"

# 读取matched文件（包含所有匹配骨架的分子）
matched_file = f"{run_folder}/matched.csv"
if not os.path.exists(matched_file):
    print(f"未找到: {matched_file}")
    sys.exit(1)

df = pd.read_csv(matched_file, low_memory=False)
print(f"匹配骨架的分子: {len(df):,}个")

# 去重后取Top20
df_unique = df.drop_duplicates(subset=['SMILES'])
top20 = df_unique.nlargest(20, 'Score')

print(f"去重后Top20: {len(top20)}个")

mols = []
legends = []

for idx, row in top20.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        # 检查骨架匹配情况
        pyrr = row.get('Pyrrolidine_Scaffold (raw)', 0)
        cyclo = row.get('Cyclobutane_Scaffold (raw)', 0)
        scaffold_type = ""
        if pyrr > 0 and cyclo > 0:
            scaffold_type = "[Both]"
        elif pyrr > 0:
            scaffold_type = "[Pyrr]"
        elif cyclo > 0:
            scaffold_type = "[Cyclo]"
        
        legend = f"{scaffold_type}\nScore: {row['Score']:.4f}\npIC50: {row['DENV_Activity_pIC50 (raw)']:.2f}"
        legends.append(legend)

if mols:
    output_file = f"{run_folder}/top20_scaffolds.png"
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(350, 350), legends=legends, returnPNG=False)
    img.save(output_file)
    print(f"已保存: {output_file}")
