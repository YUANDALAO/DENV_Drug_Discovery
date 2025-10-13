import pandas as pd

print("Run1 vs Run2 对比")
print("="*70)

# Run1
df1 = pd.read_csv('runs/run1_scaffold_original/results.csv')
pyr1 = (df1['Pyrrolidine (raw)'] > 0.5).sum()
cyc1 = (df1['Cyclobutane (raw)'] > 0.5).sum()
both1 = ((df1['Pyrrolidine (raw)'] > 0.5) & (df1['Cyclobutane (raw)'] > 0.5)).sum()

print("Run1 (原始SMARTS + geometric_mean):")
print(f"  最高Score: {df1['Score'].max():.3f}")
print(f"  平均pIC50: {df1['Activity (raw)'].mean():.2f}")
print(f"  吡咯烷匹配: {pyr1} ({pyr1/len(df1)*100:.1f}%)")
print(f"  环丁烷匹配: {cyc1} ({cyc1/len(df1)*100:.1f}%)")
print(f"  同时匹配: {both1} ({both1/len(df1)*100:.1f}%) ← 问题！")

print("\n" + "-"*70 + "\n")

# Run2
df2 = pd.read_csv('runs/run2_corrected_smarts/results.csv')
pyr2 = (df2['Pyrrolidine_Scaffold (raw)'] > 0.5).sum()
cyc2 = (df2['Cyclobutane_Scaffold (raw)'] > 0.5).sum()
both2 = ((df2['Pyrrolidine_Scaffold (raw)'] > 0.5) & 
         (df2['Cyclobutane_Scaffold (raw)'] > 0.5)).sum()

print("Run2 (修正SMARTS + arithmetic_mean):")
print(f"  最高Score: {df2['Score'].max():.3f}")
print(f"  平均pIC50: {df2['Activity (raw)'].mean():.2f}")
print(f"  吡咯烷匹配: {pyr2} ({pyr2/len(df2)*100:.1f}%)")
print(f"  环丁烷匹配: {cyc2} ({cyc2/len(df2)*100:.1f}%)")
print(f"  同时匹配: {both2} ({both2/len(df2)*100:.1f}%) ← 期望接近0")

print("\n改进:")
print(f"  Score提升: {(df2['Score'].max() - df1['Score'].max()):.3f}")
print(f"  活性提升: {(df2['Activity (raw)'].mean() - df1['Activity (raw)'].mean()):.2f}")
