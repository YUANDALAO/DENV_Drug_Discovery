from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
df = pd.read_csv('scaffold_FIXED_1.csv')
top10 = df.nlargest(10, 'Score')
mols = []
legends = []
for idx, row in top10.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legend = f"Score: {row['Score']:.3f}\npIC50: {row['Activity (raw)']:.2f}"
        legends.append(legend)
Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(400,400), legends=legends).save('top10_structures.png')
print("完成")
