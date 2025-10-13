import pandas as pd
from rdkit import Chem

def clean_smiles_data(input_file, output_file):
    """Clean and validate SMILES data for REINVENT4"""
    print(f"Reading data from {input_file}...")
    df = pd.read_csv(input_file)
    
    print("Data shape:", df.shape)
    print("Using 'Smiles' column for SMILES data")
    
    valid_smiles = []
    invalid_count = 0
    
    for idx, smiles in enumerate(df['Smiles']):
        if pd.isna(smiles):
            invalid_count += 1
            continue
            
        try:
            # Validate SMILES
            mol = Chem.MolFromSmiles(str(smiles))
            if mol is not None:
                # Canonicalize SMILES
                canonical_smiles = Chem.MolToSmiles(mol)
                valid_smiles.append(canonical_smiles)
            else:
                invalid_count += 1
        except:
            invalid_count += 1
    
    # Remove duplicates while preserving order
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    # Save as TSV format (required by REINVENT4)
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"\nProcessing complete:")
    print(f"  Original entries: {len(df)}")
    print(f"  Valid SMILES: {len(valid_smiles)}")
    print(f"  Unique SMILES: {len(unique_smiles)}")
    print(f"  Invalid/skipped: {invalid_count}")
    print(f"  Output saved to: {output_file}")

if __name__ == "__main__":
    clean_smiles_data("data/NS3.csv", "data/denv_smiles_clean.tsv")
