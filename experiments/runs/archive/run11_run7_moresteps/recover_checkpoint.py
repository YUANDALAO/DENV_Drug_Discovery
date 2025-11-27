"""
从REINVENT4 checkpoint恢复完整的生成结果
适用于 IdenticalMurckoScaffold diversity filter
"""

import torch
import pandas as pd
from pathlib import Path

def recover_from_checkpoint(checkpoint_path, output_csv):
    """
    从checkpoint文件中提取所有生成的分子
    
    Args:
        checkpoint_path: checkpoint文件路径
        output_csv: 输出CSV文件路径
    """
    print(f"正在加载 checkpoint: {checkpoint_path}")
    
    # 加载checkpoint (使用CPU以节省GPU内存)
    ckpt = torch.load(checkpoint_path, map_location='cpu')
    
    print("\n=== Checkpoint 结构 ===")
    for key in ckpt.keys():
        print(f"- {key}: {type(ckpt[key])}")
    
    # 检查staged_learning中的内容
    if 'staged_learning' in ckpt:
        print("\n=== Staged Learning 内容 ===")
        sl_data = ckpt['staged_learning']
        for key in sl_data.keys():
            print(f"- {key}: {type(sl_data[key])}")
        
        # diversity filter应该在这里
        if 'diversity_filter' in sl_data:
            df_data = sl_data['diversity_filter']
            print(f"\n=== Diversity Filter 内容 ===")
            for key in df_data.keys():
                print(f"- {key}: {type(df_data[key])}")
        else:
            print("\n❌ 错误: staged_learning中未找到 'diversity_filter'")
            return
    elif 'diversity_filter' in ckpt:
        df_data = ckpt['diversity_filter']
        print(f"\n=== Diversity Filter 内容 ===")
        for key in df_data.keys():
            print(f"- {key}: {type(df_data[key])}")
    else:
        print("\n❌ 错误: checkpoint中未找到 'diversity_filter'")
        print("尝试查找其他可能的位置...")
        
        # 递归查找diversity_filter
        def find_diversity_filter(d, path=""):
            if isinstance(d, dict):
                if 'diversity_filter' in d:
                    print(f"✓ 找到 diversity_filter 在: {path}")
                    return d['diversity_filter']
                for k, v in d.items():
                    result = find_diversity_filter(v, f"{path}.{k}" if path else k)
                    if result is not None:
                        return result
            return None
        
        df_data = find_diversity_filter(ckpt)
        if df_data is None:
            print("❌ 无法找到 diversity_filter")
            return
    
    # IdenticalMurckoScaffold 使用 'scaffolds' 存储
    if 'scaffolds' in df_data:
        scaffolds_dict = df_data['scaffolds']
        print(f"\n✓ 找到 {len(scaffolds_dict)} 个不同的 Murcko scaffolds")
        
        # 收集所有SMILES和分数
        all_molecules = []
        
        for scaffold_hash, molecules_list in scaffolds_dict.items():
            # molecules_list 是该scaffold下的分子列表
            for mol_data in molecules_list:
                # 提取SMILES和分数
                if isinstance(mol_data, dict):
                    smiles = mol_data.get('smiles', '')
                    score = mol_data.get('score', 0.0)
                elif isinstance(mol_data, tuple):
                    smiles = mol_data[0]
                    score = mol_data[1] if len(mol_data) > 1 else 0.0
                else:
                    smiles = str(mol_data)
                    score = 0.0
                
                all_molecules.append({
                    'SMILES': smiles,
                    'Score': score,
                    'Scaffold_Hash': scaffold_hash
                })
        
        print(f"✓ 总共恢复 {len(all_molecules)} 个分子")
        
        # 创建DataFrame
        df = pd.DataFrame(all_molecules)
        
        # 按分数降序排列
        df = df.sort_values('Score', ascending=False)
        
        # 保存结果
        df.to_csv(output_csv, index=False)
        print(f"\n✓ 结果已保存到: {output_csv}")
        
        # 统计信息
        print("\n=== 统计信息 ===")
        print(f"总分子数: {len(df)}")
        print(f"唯一SMILES数: {df['SMILES'].nunique()}")
        print(f"分数范围: {df['Score'].min():.4f} - {df['Score'].max():.4f}")
        print(f"平均分数: {df['Score'].mean():.4f}")
        print(f"中位数分数: {df['Score'].median():.4f}")
        print(f"\n分数分布:")
        print(df['Score'].describe())
        
        # 显示top 10
        print("\n=== Top 10 分子 ===")
        print(df.head(10).to_string(index=False))
        
    elif 'memory' in df_data:
        # 如果是其他类型的diversity filter
        print("\n使用 'memory' 字段提取数据...")
        memory = df_data['memory']
        
        molecules = []
        for smiles, score in memory.items():
            molecules.append({
                'SMILES': smiles,
                'Score': score
            })
        
        df = pd.DataFrame(molecules)
        df = df.sort_values('Score', ascending=False)
        df.to_csv(output_csv, index=False)
        
        print(f"✓ 恢复 {len(df)} 个分子")
        print(f"✓ 结果已保存到: {output_csv}")
    
    else:
        print("\n❌ 未识别的 diversity filter 数据结构")
        print("可用字段:", list(df_data.keys()))

    # 额外检查：是否有step信息
    if 'step' in ckpt:
        print(f"\n训练步数: {ckpt['step']}")


if __name__ == "__main__":
    import sys
    
    # 检查命令行参数
    if len(sys.argv) < 2:
        print("用法: python recover_checkpoint.py <实验文件夹路径>")
        print("示例: python recover_checkpoint.py experiments/runs/run11_run7_moresteps")
        sys.exit(1)
    
    # 获取实验文件夹路径
    run_dir = Path(sys.argv[1])
    
    if not run_dir.exists():
        print(f"❌ 错误: 文件夹不存在: {run_dir}")
        sys.exit(1)
    
    # 设置输入输出路径
    checkpoint_file = run_dir / "checkpoint.chkpt"
    output_file = run_dir / "recovered_all_results.csv"
    
    if not checkpoint_file.exists():
        print(f"❌ 错误: checkpoint文件不存在: {checkpoint_file}")
        sys.exit(1)
    
    print(f"实验文件夹: {run_dir}")
    print(f"Checkpoint文件: {checkpoint_file}")
    print(f"输出文件: {output_file}")
    print("")
    
    # 恢复数据
    recover_from_checkpoint(str(checkpoint_file), str(output_file))
    
    print("\n" + "="*50)
    print("恢复完成！")
    print(f"结果保存在: {output_file}")
    print("="*50)