"""
诊断REINVENT4 checkpoint的完整数据结构
"""

import torch
import sys
from pathlib import Path

def diagnose_checkpoint(checkpoint_path):
    """详细诊断checkpoint结构"""
    
    print(f"正在加载: {checkpoint_path}\n")
    ckpt = torch.load(checkpoint_path, map_location='cpu')
    
    # 1. 获取diversity filter对象
    df_obj = ckpt['staged_learning']['diversity_filter']
    
    print("="*60)
    print("DIVERSITY FILTER 对象分析")
    print("="*60)
    print(f"类型: {type(df_obj)}")
    print(f"\n可用属性:")
    for attr in dir(df_obj):
        if not attr.startswith('_'):
            print(f"  - {attr}")
    
    # 2. 检查 smiles_memory
    print("\n" + "="*60)
    print("SMILES MEMORY 详细分析")
    print("="*60)
    
    if hasattr(df_obj, 'smiles_memory'):
        smiles_mem = df_obj.smiles_memory
        print(f"类型: {type(smiles_mem)}")
        print(f"属性: {dir(smiles_mem)}")
        
        # 如果是dict
        if isinstance(smiles_mem, dict):
            print(f"\n字典长度: {len(smiles_mem)}")
            print("\n前5个条目:")
            for i, (k, v) in enumerate(list(smiles_mem.items())[:5]):
                print(f"  {k[:60]}: {v} (type: {type(v)})")
        
        # 如果是对象，检查其属性
        elif hasattr(smiles_mem, '__dict__'):
            print(f"\n对象属性: {smiles_mem.__dict__.keys()}")
            for key, value in list(smiles_mem.__dict__.items())[:5]:
                print(f"  {key}: {type(value)}")
                if isinstance(value, dict) and len(value) > 0:
                    first_k = list(value.keys())[0]
                    first_v = value[first_k]
                    print(f"    示例: {first_k[:40]}: {first_v}")
        
        # 如果有内部计数器
        if hasattr(smiles_mem, 'counter'):
            print(f"\nCounter属性: {type(smiles_mem.counter)}")
            if isinstance(smiles_mem.counter, dict):
                print(f"Counter长度: {len(smiles_mem.counter)}")
                print("\n前5个条目:")
                for i, (k, v) in enumerate(list(smiles_mem.counter.items())[:5]):
                    print(f"  {k[:60]}: {v}")
        
        # 如果有buckets
        if hasattr(smiles_mem, 'buckets'):
            print(f"\nBuckets属性: {type(smiles_mem.buckets)}")
            if isinstance(smiles_mem.buckets, dict):
                print(f"Buckets长度: {len(smiles_mem.buckets)}")
    
    # 3. 检查 scaffold_memory
    print("\n" + "="*60)
    print("SCAFFOLD MEMORY 详细分析")
    print("="*60)
    
    if hasattr(df_obj, 'scaffold_memory'):
        scaffold_mem = df_obj.scaffold_memory
        print(f"类型: {type(scaffold_mem)}")
        
        # 检查BucketCounter的内部结构
        if hasattr(scaffold_mem, '__dict__'):
            print(f"\n对象属性:")
            for key in scaffold_mem.__dict__.keys():
                attr_value = getattr(scaffold_mem, key)
                print(f"  {key}: {type(attr_value)}")
                
                # 如果是字典，显示详情
                if isinstance(attr_value, dict) and len(attr_value) > 0:
                    print(f"    长度: {len(attr_value)}")
                    first_k = list(attr_value.keys())[0]
                    first_v = attr_value[first_k]
                    print(f"    示例 - Key: {first_k}")
                    print(f"    示例 - Value: {first_v} (type: {type(first_v)})")
                    
                    # 如果value是list或dict，进一步检查
                    if isinstance(first_v, (list, tuple)) and len(first_v) > 0:
                        print(f"    第一个元素: {first_v[0]}")
                    elif isinstance(first_v, dict) and len(first_v) > 0:
                        sub_k = list(first_v.keys())[0]
                        print(f"    子Key: {sub_k}, 子Value: {first_v[sub_k]}")
    
    # 4. 检查是否有其他可能存储SMILES和分数的地方
    print("\n" + "="*60)
    print("查找所有可能包含分数的属性")
    print("="*60)
    
    def find_dict_with_smiles(obj, path="", max_depth=3):
        """递归查找包含SMILES和分数的字典"""
        if max_depth == 0:
            return
        
        if isinstance(obj, dict):
            # 检查是否看起来像SMILES: score映射
            if len(obj) > 0:
                first_k = list(obj.keys())[0]
                first_v = obj[first_k]
                if isinstance(first_k, str) and 'c' in first_k.lower() and isinstance(first_v, (int, float)):
                    print(f"\n✓ 找到可能的SMILES字典: {path}")
                    print(f"  长度: {len(obj)}")
                    print(f"  示例: {first_k[:50]} -> {first_v}")
            
            for k, v in list(obj.items())[:10]:  # 只检查前10个
                find_dict_with_smiles(v, f"{path}.{k}", max_depth-1)
        
        elif hasattr(obj, '__dict__'):
            for k, v in list(obj.__dict__.items())[:10]:
                find_dict_with_smiles(v, f"{path}.{k}", max_depth-1)
    
    find_dict_with_smiles(df_obj, "diversity_filter")
    
    print("\n" + "="*60)
    print("诊断完成")
    print("="*60)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python diagnose_checkpoint.py <checkpoint文件路径>")
        sys.exit(1)
    
    checkpoint_file = sys.argv[1]
    
    if not Path(checkpoint_file).exists():
        print(f"错误: 文件不存在: {checkpoint_file}")
        sys.exit(1)
    
    diagnose_checkpoint(checkpoint_file)