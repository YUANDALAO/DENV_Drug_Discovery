"""
使用REINVENT4的评分系统重新评分恢复的分子
"""

import sys
from pathlib import Path
import pandas as pd

def create_scoring_config(smiles_csv, output_csv, config_toml):
    """创建评分配置文件"""
    
    config = f"""
run_type = "scoring"
device = "cuda:0"
json_out_config = "{Path(config_toml).parent / '_scoring.json'}"

[parameters]
smiles_file = "{smiles_csv}"
output_csv = "{output_csv}"

[scoring]
type = "arithmetic_mean"

# 使用与训练时完全相同的评分组件
[[scoring.component]]
[scoring.component.QSARScorer]
[[scoring.component.QSARScorer.endpoint]]
name = "DENV_Activity"
weight = 1.5
[scoring.component.QSARScorer.endpoint.params]
model_path = "/mnt/c/Users/ucsaheu/python_projects/DENV_Drug_Discovery/02_ML_Modeling/models/random_forest_champion.joblib"
[scoring.component.QSARScorer.endpoint.transform]
type = "double_sigmoid"
low = 5.0
high = 9.0
coef_div = 9.0
coef_si = 10.0
coef_se = 10.0

[[scoring.component]]
[scoring.component.NumAtomStereoCenters]
[[scoring.component.NumAtomStereoCenters.endpoint]]
name = "Stereo_Centers"
weight = 1.0
[scoring.component.NumAtomStereoCenters.endpoint.transform]
type = "reverse_sigmoid"
low = 0
high = 2
k = 1.5

[[scoring.component]]
[scoring.component.CustomAlerts]
[[scoring.component.CustomAlerts.endpoint]]
name = "Stability_Alerts"
weight = 1.0
[scoring.component.CustomAlerts.endpoint.params]
smarts = [
    "[*;r8]", "[*;r9]", "[*;r10]",
    "[#8][#8]", "[#6;+]", "C#C", "[NX3][NX3]", "[SH]",
    "C(=O)OC(=O)", "[N+](=O)[O-]", "S(=O)(=O)Cl",
    "[F,Cl,Br,I][C,c][F,Cl,Br,I]",
    "c1ccc2c(c1)C(=O)c1c(C2=O)cccc1",
    "C1=CC(=O)C=CC1=O",
    "[C;R1]1[C;R1][C;R1][N;R1][C;R1]1([!#1])([!#1])([!#1])",
    "[C;R1]1[C;R1]([!#1])([!#1])[C;R1]([!#1])([!#1])[N;R1][C;R1]1",
    "C(=O)N([!#1])C(=O)",
]

[[scoring.component]]
[scoring.component.Qed]
[[scoring.component.Qed.endpoint]]
name = "QED"
weight = 0.5

[[scoring.component]]
[scoring.component.SAScore]
[[scoring.component.SAScore.endpoint]]
name = "SA"
weight = 0.6
[scoring.component.SAScore.endpoint.transform]
type = "reverse_sigmoid"
low = 1.0
high = 4.5
k = 1.0

[[scoring.component]]
[scoring.component.MolecularWeight]
[[scoring.component.MolecularWeight.endpoint]]
name = "MW"
weight = 0.4
[scoring.component.MolecularWeight.endpoint.transform]
type = "double_sigmoid"
low = 300.0
high = 500.0
coef_div = 500.0
coef_si = 20.0
coef_se = 20.0

[[scoring.component]]
[scoring.component.NumHeavyAtoms]
[[scoring.component.NumHeavyAtoms.endpoint]]
name = "Heavy_Atoms"
weight = 0.3
[scoring.component.NumHeavyAtoms.endpoint.transform]
type = "double_sigmoid"
low = 20
high = 35
coef_div = 35.0
coef_si = 5.0
coef_se = 5.0

[[scoring.component]]
[scoring.component.SlogP]
[[scoring.component.SlogP.endpoint]]
name = "LogP"
weight = 0.4
[scoring.component.SlogP.endpoint.transform]
type = "double_sigmoid"
low = 1.0
high = 4.0
coef_div = 4.0
coef_si = 10.0
coef_se = 10.0

[[scoring.component]]
[scoring.component.TPSA]
[[scoring.component.TPSA.endpoint]]
name = "TPSA"
weight = 0.3
[scoring.component.TPSA.endpoint.transform]
type = "double_sigmoid"
low = 40.0
high = 100.0
coef_div = 100.0
coef_si = 20.0
coef_se = 20.0

[[scoring.component]]
[scoring.component.HBondAcceptors]
[[scoring.component.HBondAcceptors.endpoint]]
name = "HBA"
weight = 0.4
[scoring.component.HBondAcceptors.endpoint.transform]
type = "reverse_sigmoid"
low = 2
high = 7
k = 0.5

[[scoring.component]]
[scoring.component.HBondDonors]
[[scoring.component.HBondDonors.endpoint]]
name = "HBD"
weight = 0.3
[scoring.component.HBondDonors.endpoint.transform]
type = "reverse_sigmoid"
low = 0
high = 3
k = 0.5

[[scoring.component]]
[scoring.component.NumRotBond]
[[scoring.component.NumRotBond.endpoint]]
name = "Rotatable_Bonds"
weight = 0.3
[scoring.component.NumRotBond.endpoint.transform]
type = "reverse_sigmoid"
low = 0
high = 7
k = 0.5

[[scoring.component]]
[scoring.component.NumRings]
[[scoring.component.NumRings.endpoint]]
name = "Ring_Count"
weight = 0.3
[scoring.component.NumRings.endpoint.transform]
type = "reverse_sigmoid"
low = 2
high = 4
k = 0.5

[[scoring.component]]
[scoring.component.NumAromaticRings]
[[scoring.component.NumAromaticRings.endpoint]]
name = "Aromatic_Rings"
weight = 0.2
[scoring.component.NumAromaticRings.endpoint.transform]
type = "reverse_sigmoid"
low = 1
high = 3
k = 0.5
"""
    
    with open(config_toml, 'w') as f:
        f.write(config)
    
    print(f"✓ 评分配置已保存到: {config_toml}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python rescore_molecules.py <实验文件夹路径>")
        print("示例: python rescore_molecules.py experiments/runs/run11_run7_moresteps")
        sys.exit(1)
    
    run_dir = Path(sys.argv[1])
    
    # 1. 准备SMILES文件（只要第一列）
    recovered_csv = run_dir / "recovered_all_results.csv"
    smiles_only = run_dir / "smiles_for_scoring.smi"
    
    if not recovered_csv.exists():
        print(f"❌ 错误: 未找到恢复的CSV文件: {recovered_csv}")
        print("请先运行: python recover_checkpoint.py " + str(run_dir))
        sys.exit(1)
    
    # 读取SMILES
    df = pd.read_csv(recovered_csv)
    print(f"✓ 读取了 {len(df)} 个分子")
    
    # 保存为.smi格式（只有SMILES列）
    df[['SMILES']].to_csv(smiles_only, index=False, header=False)
    print(f"✓ SMILES已保存到: {smiles_only}")
    
    # 2. 创建评分配置
    scoring_config = run_dir / "scoring_config.toml"
    rescored_csv = run_dir / "rescored_results.csv"
    
    create_scoring_config(
        smiles_csv=str(smiles_only),
        output_csv=str(rescored_csv),
        config_toml=str(scoring_config)
    )
    
    print("\n" + "="*60)
    print("准备完成！现在运行:")
    print("="*60)
    print(f"reinvent {scoring_config}")
    print("\n评分完成后，结果将保存在:")
    print(f"  {rescored_csv}")
    print("="*60)