cat > visualize_top10.py << 'ENDOFFILE'
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

df = pd.read_csv('run2_final_1.csv')
top10 = df.nlargest(10, 'Score')

mols = []
legends = []

for idx, row in top10.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legend = f"Score: {row['Score']:.3f}\npIC50: {row['Activity (raw)']:.2f}\nQED: {row['QED (raw)']:.2f}\nSA: {row['SA (raw)']:.1f}"
        legends.append(legend)

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(400, 400),
                            legends=legends, returnPNG=False)
img.save('top10_structures.png')

with open('top10_smiles.txt', 'w') as f:
    for i, (idx, row) in enumerate(top10.iterrows(), 1):
        f.write(f"#{i} Score={row['Score']:.3f} pIC50={row['Activity (raw)']:.2f}\n")
        f.write(f"{row['SMILES']}\n\n")

print("完成！")
print("  top10_structures.png")
print("  top10_smiles.txt")
