import pandas as pd
from rdkit import Chem

# 读取Run3的优质候选
df = pd.read_csv('experiments/runs/run3/promising_candidates.csv', low_memory=False)
print(f"Run3优质候选分子总数: {len(df)}")

# 目标骨架SMARTS
pyrr_a = Chem.MolFromSmarts("C(N)a")
pyrr_b = Chem.MolFromSmarts("C(=O)N")
cyclo_a = Chem.MolFromSmarts("[C;R]a")
cyclo_b = Chem.MolFromSmarts("C(=O)N")

print("\n检查骨架匹配情况...")

results = []
for idx, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    if mol:
        # 检查吡咯烷特征
        pyrr_match = (mol.HasSubstructMatch(pyrr_a) and 
                     mol.HasSubstructMatch(pyrr_b))
        # 检查环丁烷特征
        cyclo_match = (mol.HasSubstructMatch(cyclo_a) and 
                      mol.HasSubstructMatch(cyclo_b))
        
        if pyrr_match or cyclo_match:
            results.append({
                'SMILES': row['SMILES'],
                'Score': row['Score'],
                'pIC50': row['DENV_Activity_pIC50 (raw)'],
                'QED': row['Drug_Likeness (raw)'],
                'SA': row['Synthetic_Accessibility (raw)'],
                'Type': 'Pyrr' if pyrr_match else 'Cyclo'
            })

print(f"\n匹配骨架特征的分子: {len(results)} / {len(df)} ({len(results)/len(df)*100:.1f}%)")

if results:
    result_df = pd.DataFrame(results)
    
    print("\n按pIC50排序的Top20:")
    top20 = result_df.nlargest(20, 'pIC50')
    print(top20.to_string(index=False))
    
    # 保存
    result_df.to_csv('experiments/runs/run3/scaffold_matched_candidates.csv', index=False)
    print(f"\n已保存: experiments/runs/run3/scaffold_matched_candidates.csv")
    
    # 统计
    print("\n统计:")
    print(f"  吡咯烷型: {(result_df['Type']=='Pyrr').sum()}")
    print(f"  环丁烷型: {(result_df['Type']=='Cyclo').sum()}")
    print(f"  平均pIC50: {result_df['pIC50'].mean():.2f}")
    print(f"  最高pIC50: {result_df['pIC50'].max():.2f}")
else:
    print("\n没有找到匹配目标骨架特征的分子。")
    print("\n这说明:")
    print("1. Run3生成的分子虽然满足药性指标")
    print("2. 但不包含您指定的骨架结构")
    print("3. REINVENT自由生成无法精确控制骨架")
