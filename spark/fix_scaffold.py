#!/usr/bin/env python3
"""
修复Scaffold SMILES文件中的[*]问题
将 [*] 替换为 * 以匹配LibInvent词汇表
"""

import sys
from pathlib import Path

def fix_scaffold_smiles(input_file, output_file=None):
    """
    修复scaffold SMILES文件
    
    Args:
        input_file: 输入文件路径
        output_file: 输出文件路径（如果为None，则覆盖原文件）
    """
    input_path = Path(input_file)
    
    if not input_path.exists():
        print(f"❌ 文件不存在: {input_file}")
        return False
    
    # 读取文件
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    print(f"📖 读取文件: {input_file}")
    print(f"   行数: {len(lines)}")
    
    # 修复SMILES
    fixed_lines = []
    changes = 0
    
    for i, line in enumerate(lines):
        original = line.strip()
        
        # 替换 [*] 为 *
        fixed = original.replace('[*]', '*')
        
        if fixed != original:
            changes += 1
            print(f"   第{i+1}行: {original[:50]} → {fixed[:50]}")
        
        fixed_lines.append(fixed + '\n')
    
    print(f"\n✏️  修改统计:")
    print(f"   总行数: {len(lines)}")
    print(f"   修改行数: {changes}")
    
    # 写入文件
    if output_file is None:
        # 备份原文件
        backup_path = input_path.with_suffix('.smi.backup')
        input_path.rename(backup_path)
        print(f"\n💾 备份原文件: {backup_path}")
        output_file = input_file
    
    output_path = Path(output_file)
    with open(output_path, 'w') as f:
        f.writelines(fixed_lines)
    
    print(f"✅ 修复完成: {output_file}")
    
    # 验证
    print(f"\n🔍 验证修复后的文件:")
    with open(output_path, 'r') as f:
        sample = f.readlines()[:5]
    
    for i, line in enumerate(sample, 1):
        has_bracket_star = '[*]' in line
        status = '❌' if has_bracket_star else '✅'
        print(f"   {status} 第{i}行: {line.strip()[:60]}")
    
    return True


if __name__ == "__main__":
    # 默认处理标准位置的文件
    scaffold_file = "data/pyrrolidine_dual_aryl.smi"
    
    if len(sys.argv) > 1:
        scaffold_file = sys.argv[1]
    
    print("="*80)
    print("Scaffold SMILES修复工具")
    print("="*80)
    print()
    
    success = fix_scaffold_smiles(scaffold_file)
    
    if success:
        print()
        print("="*80)
        print("✅ 修复成功！")
        print("="*80)
        print()
        print("现在可以运行REINVENT4了：")
        print("  reinvent experiments/runs/spark_run1/config.toml")
    else:
        print()
        print("="*80)
        print("❌ 修复失败")
        print("="*80)
        sys.exit(1)
