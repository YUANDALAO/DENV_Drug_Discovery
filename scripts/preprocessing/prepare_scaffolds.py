from rdkit import Chem

# 定义两个核心骨架（用[*]标记可变位置）
scaffolds = [
    # 骨架1：吡咯烷骨架
    "O=C([*:1])C1C([*:2])NC([*:3])C1",
    
    # 骨架2：环丁烷骨架  
    "O=C([*:1])C1C([*:2])C([*:3])C1",
]

# 验证并保存
with open('data/denv_scaffolds.smi', 'w') as f:
    for i, scaffold in enumerate(scaffolds, 1):
        mol = Chem.MolFromSmiles(scaffold)
        if mol:
            canonical = Chem.MolToSmiles(mol)
            f.write(f"{canonical}\tScaffold_{i}\n")
            print(f"骨架 {i}: {canonical}")
        else:
            print(f"骨架 {i} 解析失败")

print(f"\n已保存 {len(scaffolds)} 个骨架到 data/denv_scaffolds.smi")
