import pandas as pd
from rdkit.Chem import Draw

top500 = pd.read_csv("top500_candidates.csv")

# 生成前50个的结构图
mols = []
legends = []

for idx, row in top500.head(50).iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legends.append(f"QED={row['QED']:.3f}\nMW={row['MW']:.0f}")

img = Draw.MolsToGridImage(
    mols, 
    molsPerRow=5, 
    subImgSize=(300, 300),
    legends=legends,
    returnPNG=False
)

img.save("top50_structures.png")
print("✓ Top 50结构图已保存: top50_structures.png")