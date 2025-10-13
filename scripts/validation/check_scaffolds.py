from rdkit import Chem
import csv

# 目标骨架
scaffold1_smiles = "C1=CC=CC=C1C2CC(C(=O)N)CN2"
scaffold2_smiles = "c1ccccc1C1CC(C1)C(=O)N"

scaffold1 = Chem.MolFromSmiles(scaffold1_smiles)
scaffold2 = Chem.MolFromSmiles(scaffold2_smiles)

print("=== 检查骨架匹配 ===")
print(f"骨架1: {scaffold1_smiles}")
print(f"骨架2: {scaffold2_smiles}\n")

# 读取CSV
matching_s1 = 0
matching_s2 = 0
total_checked = 0
examples_s1 = []
examples_s2 = []

with open('denv_rl_full_results_1.csv', 'r') as f:
    reader = csv.DictReader(f)
    
    for i, row in enumerate(reader):
        if i >= 10000:  # 检查前10000个
            break
            
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            total_checked += 1
            
            # 检查骨架1
            if mol.HasSubstructMatch(scaffold1):
                matching_s1 += 1
                if len(examples_s1) < 5:
                    score = row.get('Score', 'N/A')
                    pic50 = row.get('DENV_Activity_pIC50', 'N/A')
                    examples_s1.append((smiles, score, pic50))
            
            # 检查骨架2
            if mol.HasSubstructMatch(scaffold2):
                matching_s2 += 1
                if len(examples_s2) < 5:
                    score = row.get('Score', 'N/A')
                    pic50 = row.get('DENV_Activity_pIC50', 'N/A')
                    examples_s2.append((smiles, score, pic50))

print(f"总共检查: {total_checked} 个分子")
print(f"\n骨架1 (吡咯烷) 匹配: {matching_s1} ({matching_s1/total_checked*100:.2f}%)")
print(f"骨架2 (环丁烷) 匹配: {matching_s2} ({matching_s2/total_checked*100:.2f}%)\n")

if examples_s1:
    print("=== 骨架1匹配示例 ===")
    for smiles, score, pic50 in examples_s1:
        print(f"  SMILES: {smiles}")
        print(f"  Score: {score}, pIC50: {pic50}\n")
else:
    print("❌ 没有找到骨架1的匹配分子！\n")

if examples_s2:
    print("=== 骨架2匹配示例 ===")
    for smiles, score, pic50 in examples_s2:
        print(f"  SMILES: {smiles}")
        print(f"  Score: {score}, pIC50: {pic50}\n")
else:
    print("❌ 没有找到骨架2的匹配分子！\n")

# 额外：检查最常见的骨架
print("=== 分析最常见的核心结构 ===")
from collections import Counter
from rdkit.Chem.Scaffolds import MurckoScaffold

scaffold_counts = Counter()

with open('denv_rl_full_results_1.csv', 'r') as f:
    reader = csv.DictReader(f)
    for i, row in enumerate(reader):
        if i >= 1000:
            break
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol:
            try:
                core = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_counts[Chem.MolToSmiles(core)] += 1
            except:
                pass

print("Top 10 最常见骨架:")
for scaffold, count in scaffold_counts.most_common(10):
    print(f"  {count:4d}x  {scaffold}")
