import pandas as pd

df = pd.read_csv('experiments/runs/run3/matched.csv', low_memory=False)
df_unique = df.drop_duplicates(subset=['SMILES'])
top20 = df_unique.nlargest(20, 'Score')

print("="*80)
print("Top20 匹配骨架分子的性质分析")
print("="*80)

cols = ['Score', 'DENV_Activity_pIC50 (raw)', 'Drug_Likeness (raw)', 
        'Synthetic_Accessibility (raw)', 'Molecular_Weight (raw)', 'LogP (raw)']

stats = top20[cols].describe()
print("\n统计摘要:")
print(stats.to_string())

print("\n\n详细列表:")
display_cols = ['SMILES', 'Score', 'DENV_Activity_pIC50 (raw)', 
                'Drug_Likeness (raw)', 'Synthetic_Accessibility (raw)', 
                'Pyrrolidine_Scaffold (raw)', 'Cyclobutane_Scaffold (raw)']
print(top20[display_cols].to_string(index=False))

print("\n骨架类型统计:")
pyrr_only = ((top20['Pyrrolidine_Scaffold (raw)'] > 0) & (top20['Cyclobutane_Scaffold (raw)'] == 0)).sum()
cyclo_only = ((top20['Cyclobutane_Scaffold (raw)'] > 0) & (top20['Pyrrolidine_Scaffold (raw)'] == 0)).sum()
both = ((top20['Pyrrolidine_Scaffold (raw)'] > 0) & (top20['Cyclobutane_Scaffold (raw)'] > 0)).sum()
print(f"  仅吡咯烷: {pyrr_only}")
print(f"  仅环丁烷: {cyclo_only}")
print(f"  两者都有: {both}")
