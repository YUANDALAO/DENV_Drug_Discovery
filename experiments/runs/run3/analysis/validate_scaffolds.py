#!/usr/bin/env python3
"""
Run3 骨架验证脚本
验证生成的分子是否符合严格的骨架定义
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import sys
from pathlib import Path

# SMARTS模式定义
SCAFFOLD_PATTERNS = {
    'Pyrrolidine_Strict': 'C1N([!H])C(c2ccccc2)C(C(=O)N)C1',
    'Cyclobutane_Strict': 'C1C(c2ccccc2)C([!H])C1C(=O)N'
}

def validate_molecule(smiles, patterns):
    """验证分子是否匹配任一骨架"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES"
    
    matches = {}
    for name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        
        has_match = mol.HasSubstructMatch(pattern)
        matches[name] = has_match
    
    return matches, None

def main():
    # 查找最新的CSV文件
    csv_dir = Path("results/csv")
    csv_files = list(csv_dir.glob("run3_scaffold_strict*.csv"))
    
    if not csv_files:
        print("❌ 未找到结果CSV文件")
        sys.exit(1)
    
    latest_csv = max(csv_files, key=lambda p: p.stat().st_mtime)
    print(f"📂 分析文件: {latest_csv.name}")
    print("=" * 60)
    
    # 读取数据
    df = pd.read_csv(latest_csv)
    print(f"📊 总分子数: {len(df)}")
    
    # 验证每个分子
    results = []
    for idx, row in df.iterrows():
        smiles = row.get('SMILES', '')
        matches, error = validate_molecule(smiles, SCAFFOLD_PATTERNS)
        
        if error:
            results.append({
                'SMILES': smiles,
                'Valid': False,
                'Error': error,
                'Pyrrolidine': False,
                'Cyclobutane': False
            })
        else:
            results.append({
                'SMILES': smiles,
                'Valid': True,
                'Error': None,
                'Pyrrolidine': matches.get('Pyrrolidine_Strict', False),
                'Cyclobutane': matches.get('Cyclobutane_Strict', False)
            })
    
    # 统计结果
    results_df = pd.DataFrame(results)
    
    print("\n🔍 骨架匹配统计:")
    print("-" * 60)
    print(f"有效分子: {results_df['Valid'].sum()}")
    print(f"吡咯烷骨架: {results_df['Pyrrolidine'].sum()}")
    print(f"环丁烷骨架: {results_df['Cyclobutane'].sum()}")
    print(f"任一骨架: {(results_df['Pyrrolidine'] | results_df['Cyclobutane']).sum()}")
    print(f"无骨架匹配: {(~results_df['Pyrrolidine'] & ~results_df['Cyclobutane'] & results_df['Valid']).sum()}")
    
    # 保存验证结果
    output_file = latest_csv.parent / f"{latest_csv.stem}_validated.csv"
    results_df.to_csv(output_file, index=False)
    print(f"\n💾 验证结果已保存: {output_file.name}")
    
    # 绘制示例分子
    print("\n🎨 绘制骨架示例...")
    pyrrolidine_mols = results_df[results_df['Pyrrolidine']].head(5)
    cyclobutane_mols = results_df[results_df['Cyclobutane']].head(5)
    
    if len(pyrrolidine_mols) > 0:
        mols = [Chem.MolFromSmiles(s) for s in pyrrolidine_mols['SMILES'] if Chem.MolFromSmiles(s)]
        if mols:
            img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 300),
                                         legends=['Pyrrolidine'] * len(mols))
            img.save('results/images/pyrrolidine_examples.png')
            print(f"  ✅ 吡咯烷示例: results/images/pyrrolidine_examples.png")
    
    if len(cyclobutane_mols) > 0:
        mols = [Chem.MolFromSmiles(s) for s in cyclobutane_mols['SMILES'] if Chem.MolFromSmiles(s)]
        if mols:
            img = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(300, 300),
                                         legends=['Cyclobutane'] * len(mols))
            img.save('results/images/cyclobutane_examples.png')
            print(f"  ✅ 环丁烷示例: results/images/cyclobutane_examples.png")
    
    print("\n" + "=" * 60)
    print("✅ 验证完成")

if __name__ == "__main__":
    main()
