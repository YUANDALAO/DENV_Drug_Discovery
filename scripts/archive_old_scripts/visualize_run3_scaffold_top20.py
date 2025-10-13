from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd

df = pd.read_csv('experiments/runs/run3/scaffold_matched_candidates.csv')
top20 = df.nlargest(20, 'pIC50').drop_duplicates(subset=['SMILES']).head(20)

mols = []
legends = []

for idx, row in top20.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        mols.append(mol)
        legend = f"{row['Type']}\npIC50: {row['pIC50']:.2f}\nQED: {row['QED']:.2f}"
        legends.append(legend)

if mols:
    img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(350, 350), 
                               legends=legends, returnPNG=False)
    img.save('experiments/runs/run3/scaffold_matched_top20.png')
    print(f"Run3骨架匹配Top20已保存")
