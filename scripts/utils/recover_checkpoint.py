"""
从REINVENT4 checkpoint恢复完整的生成结果
适用于 IdenticalMurckoScaffold diversity filter
"""

import torch
import pandas as pd
from pathlib import Path
import sys


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
    if 'staged_learning' not in ckpt:
        print("\n❌ 错误: checkpoint中未找到 'staged_learning'")
        return
    
    print("\n=== Staged Learning 内容 ===")
    sl_data = ckpt['staged_learning']
    for key in sl_data.keys():
        print(f"- {key}: {type(sl_data[key])}")
    
    # diversity filter应该在这里
    if 'diversity_filter' not in sl_data:
        print("\n❌ 错误: staged_learning中未找到 'diversity_filter'")
        return
    
    df_obj = sl_data['diversity_filter']
    print(f"\n=== Diversity Filter 对象 ===")
    print(f"类型: {type(df_obj)}")
    
    # smiles_memory是一个set，只包含SMILES字符串
    if not hasattr(df_obj, 'smiles_memory'):
        print("\n❌ 错误: diversity_filter对象没有 'smiles_memory' 属性")
        return
    
    smiles_memory = df_obj.smiles_memory
    print(f"\n✓ 找到 smiles_memory")
    print(f"类型: {type(smiles_memory)}")
    
    if not isinstance(smiles_memory, set):
        print(f"⚠️  smiles_memory不是set类型: {type(smiles_memory)}")
        return
    
    print(f"分子数量: {len(smiles_memory)}")
    
    # 提取所有SMILES
    all_molecules = []
    for smiles in smiles_memory:
        all_molecules.append({
            'SMILES': smiles,
            'Note': 'Score needs to be recalculated'
        })
    
    print(f"\n✓ 提取了 {len(all_molecules)} 个唯一SMILES")
    print("\n⚠️  注意: Checkpoint中未保存分数信息")
    print("   需要使用REINVENT4的scoring模式重新评分")
    
    # 保存SMILES
    df = pd.DataFrame(all_molecules)
    df.to_csv(output_csv, index=False)
    
    print(f"\n✓ SMILES已保存到: {output_csv}")
    
    # 显示统计
    print("\n=== 统计信息 ===")
    print(f"唯一分子数: {len(df)}")
    print(f"\n前10个分子:")
    for i, smiles in enumerate(list(smiles_memory)[:10]):
        print(f"  {i+1}. {smiles}")
    
    # 额外检查：是否有step信息
    if 'metadata' in ckpt and 'step' in ckpt['metadata']:
        print(f"\n训练步数: {ckpt['metadata']['step']}")


if __name__ == "__main__":
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
    output_file = run_dir / "recovered_smiles.csv"
    
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
    print(f"SMILES已保存在: {output_file}")
    print("\n下一步: 运行重新评分")
    print(f"  python rescore_molecules.py {run_dir}")
    print("="*50)