import torch
from reinvent.models.model_factory.libinvent_adapter import LibinventAdapter

# 直接加载prior
model_path = "priors/libinvent.prior"
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# 加载checkpoint
checkpoint = torch.load(model_path, map_location=device)

print("=== Checkpoint Keys ===")
print(checkpoint.keys())

print("\n=== Model Type ===")
print(checkpoint.get('model_type', 'NOT FOUND'))

print("\n=== Vocabulary ===")
if 'vocabulary' in checkpoint:
    vocab = checkpoint['vocabulary']
    print("Decoration tokens:", vocab.decoration_vocabulary.tokens())
else:
    print("Vocabulary not in checkpoint")

# 测试tokenization
from rdkit import Chem
test_smiles = ["c1ccccc1*", "c1ccccc1[*]", "CCc1ccccc1*"]
print("\n=== RDKit Parsing ===")
for smi in test_smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol:
        canonical = Chem.MolToSmiles(mol)
        print(f"{smi:20s} -> {canonical}")
    else:
        print(f"{smi:20s} -> INVALID")