from Bio.PDB import PDBParser
import numpy as np

# 读取你的AlphaFold预测结构
parser = PDBParser()
structure = parser.get_structure('DENV2_AF', 'DENV2_alphafold.pdb')

# 提取pLDDT分数（存在B-factor列）
plddt_scores = []
residue_numbers = []

for model in structure:
    for chain in model:
        for residue in chain:
            try:
                ca = residue['CA']
                plddt_scores.append(ca.bfactor)
                residue_numbers.append(residue.id[1])
            except:
                continue

plddt_array = np.array(plddt_scores)

print("="*60)
print("AlphaFold 质量评估")
print("="*60)
print(f"平均 pLDDT: {plddt_array.mean():.2f}")
print(f"最低 pLDDT: {plddt_array.min():.2f}")
print(f"最高 pLDDT: {plddt_array.max():.2f}")
print(f"\npLDDT分布:")
print(f"  >90 (优秀):  {(plddt_array > 90).sum()}/{len(plddt_array)} ({(plddt_array > 90).sum()/len(plddt_array)*100:.1f}%)")
print(f"  70-90 (良好): {((plddt_array > 70) & (plddt_array <= 90)).sum()}/{len(plddt_array)} ({((plddt_array > 70) & (plddt_array <= 90)).sum()/len(plddt_array)*100:.1f}%)")
print(f"  50-70 (一般): {((plddt_array > 50) & (plddt_array <= 70)).sum()}/{len(plddt_array)} ({((plddt_array > 50) & (plddt_array <= 70)).sum()/len(plddt_array)*100:.1f}%)")
print(f"  <50 (差):    {(plddt_array <= 50).sum()}/{len(plddt_array)} ({(plddt_array <= 50).sum()/len(plddt_array)*100:.1f}%)")

# 关键：检查活性位点区域的pLDDT
# 假设催化三联体在残基50-140区域
active_site_region = (np.array(residue_numbers) >= 50) & (np.array(residue_numbers) <= 140)
active_site_plddt = plddt_array[active_site_region]

print(f"\n活性位点区域 (aa 50-140) pLDDT: {active_site_plddt.mean():.2f}")
if active_site_plddt.mean() > 80:
    print("✓ 活性位点质量优秀，可用于对接")
elif active_site_plddt.mean() > 70:
    print("⚠ 活性位点质量良好，建议与4M9K对比")
else:
    print("✗ 活性位点质量不佳，建议用4M9K")