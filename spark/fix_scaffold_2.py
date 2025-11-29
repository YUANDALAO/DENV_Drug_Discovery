#!/usr/bin/env python3
"""
修复scaffold文件中的attachment point格式
将 [*] 转换为 * （LibInvent vocabulary兼容格式）
"""

import sys
from pathlib import Path

def fix_scaffold_file(input_path, output_path=None):
    """修复scaffold文件中的[*]为*"""
    
    if output_path is None:
        output_path = input_path
    
    input_path = Path(input_path)
    
    if not input_path.exists():
        print(f"❌ 文件不存在: {input_path}")
        return False
    
    # 读取原始内容
    with open(input_path, 'r') as f:
        original = f.read()
    
    # 统计[*]出现次数
    bracket_count = original.count('[*]')
    
    if bracket_count == 0:
        print(f"✅ 文件已经是正确格式（无[*]）: {input_path}")
        print(f"   裸*出现次数: {original.count('*')}")
        return True
    
    # 替换 [*] 为 *
    fixed = original.replace('[*]', '*')
    
    # 写入修复后的内容
    with open(output_path, 'w') as f:
        f.write(fixed)
    
    print(f"✅ 已修复 {bracket_count} 处 [*] -> *")
    print(f"   输入: {input_path}")
    print(f"   输出: {output_path}")
    
    # 显示前几行
    lines = fixed.strip().split('\n')[:5]
    print(f"\n前{len(lines)}行内容:")
    for i, line in enumerate(lines, 1):
        print(f"  {i}: {line}")
    
    return True

def check_vocabulary_compatibility(smiles):
    """检查SMILES是否与LibInvent vocabulary兼容"""
    
    # LibInvent支持的tokens
    supported = {
        '<pad>', '$', '^', '#', '(', ')', '*', '-', 
        '1', '2', '3', '4', '5', '6', '=', 
        'Br', 'C', 'Cl', 'F', 'N', 'O', 'S',
        '[N+]', '[N-]', '[N]', '[O-]', '[O]', '[S+]', 
        '[n+]', '[nH]', '[s+]', 
        'c', 'n', 'o', 's', '|'
    }
    
    # 检查是否包含不支持的[*]
    if '[*]' in smiles:
        return False, "[*] 不支持，应使用裸 *"
    
    return True, "格式正确"

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python fix_scaffold.py <scaffold_file.smi> [output_file.smi]")
        print("\n示例:")
        print("  python fix_scaffold.py data/pyrrolidine_dual_aryl.smi")
        print("  python fix_scaffold.py data/scaffold.smi data/scaffold_fixed.smi")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    fix_scaffold_file(input_file, output_file)
