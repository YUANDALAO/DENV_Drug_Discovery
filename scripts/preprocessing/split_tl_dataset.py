#!/usr/bin/env python3
"""
分割LibInvent训练集和验证集
用于Transfer Learning
"""

import random
from pathlib import Path

def split_dataset(input_file, train_ratio=0.8, random_seed=42):
    """
    分割数据集为训练集和验证集
    
    Args:
        input_file: 输入SMILES文件
        train_ratio: 训练集比例（默认0.8）
        random_seed: 随机种子
    """
    # 读取数据
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    print(f"📂 读取文件: {input_file}")
    print(f"   总行数: {len(lines)}")
    
    # 去除空行
    lines = [line for line in lines if line.strip()]
    print(f"   有效行数: {len(lines)}")
    
    # 随机打乱（保证可重复）
    random.seed(random_seed)
    random.shuffle(lines)
    
    # 分割
    split_idx = int(len(lines) * train_ratio)
    train_lines = lines[:split_idx]
    valid_lines = lines[split_idx:]
    
    print(f"\n📊 数据分割:")
    print(f"   训练集: {len(train_lines)} ({len(train_lines)/len(lines)*100:.1f}%)")
    print(f"   验证集: {len(valid_lines)} ({len(valid_lines)/len(lines)*100:.1f}%)")
    
    # 保存训练集
    input_path = Path(input_file)
    train_file = input_path.parent / f"{input_path.stem}_train{input_path.suffix}"
    valid_file = input_path.parent / f"{input_path.stem}_valid{input_path.suffix}"
    
    with open(train_file, 'w') as f:
        f.writelines(train_lines)
    
    with open(valid_file, 'w') as f:
        f.writelines(valid_lines)
    
    print(f"\n✅ 保存完成:")
    print(f"   训练集: {train_file}")
    print(f"   验证集: {valid_file}")
    
    # 显示前几行样例
    print(f"\n📋 训练集样例（前3行）:")
    for line in train_lines[:3]:
        print(f"   {line.strip()}")
    
    print(f"\n📋 验证集样例（前3行）:")
    for line in valid_lines[:3]:
        print(f"   {line.strip()}")
    
    return train_file, valid_file

def validate_libinvent_format(file_path):
    """验证LibInvent格式是否正确"""
    print(f"\n🔍 验证LibInvent格式: {file_path}")
    
    errors = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            # LibInvent格式: scaffold R1|R2|R3...
            parts = line.split()
            if len(parts) < 2:
                errors.append(f"  行{i}: 缺少R-groups - {line[:50]}")
                continue
            
            scaffold = parts[0]
            rgroups = parts[1]
            
            # 检查scaffold是否有attachment points [*]
            if '[*]' not in scaffold:
                errors.append(f"  行{i}: scaffold缺少[*] - {scaffold[:50]}")
            
            # 检查R-groups格式
            if '|' not in rgroups and len(rgroups) > 5:
                errors.append(f"  行{i}: R-groups格式可能错误 - {rgroups[:30]}")
            
            if i <= 3:
                print(f"  ✓ 行{i}: {line[:60]}...")
    
    if errors:
        print(f"\n⚠️  发现 {len(errors)} 个格式问题:")
        for err in errors[:5]:  # 只显示前5个
            print(err)
        if len(errors) > 5:
            print(f"  ... 还有 {len(errors)-5} 个问题")
    else:
        print(f"  ✅ 格式检查通过！")
    
    return len(errors) == 0

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("用法: python split_tl_dataset.py <input_file>")
        print("示例: python split_tl_dataset.py data/libinvent_final_format.smi")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # 验证格式
    if not validate_libinvent_format(input_file):
        print("\n❌ 格式验证失败！请修复后再分割")
        sys.exit(1)
    
    # 分割数据集
    train_file, valid_file = split_dataset(input_file)
    
    print("\n" + "="*60)
    print("下一步: 使用改进的TL配置重新训练")
    print("="*60)
    print("\n修改配置文件:")
    print(f'  smiles_file = "{train_file}"')
    print(f'  validation_smiles_file = "{valid_file}"')
    print("\n然后运行:")
    print("  reinvent libinvent_transfer_learning_v2.toml")