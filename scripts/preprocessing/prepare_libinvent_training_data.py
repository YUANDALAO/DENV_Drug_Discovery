import pandas as pd
from rdkit import Chem
from rdkit.Chem import BRICS

def decompose_molecule_for_libinvent(smiles):
    """将分子分解为骨架和装饰基团"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    
    # 使用BRICS分解
    try:
        # 尝试找到环系统作为骨架
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            # 简化版：将最大的环系统作为骨架
            # 其他部分作为装饰基团
            
            # 这里使用简化策略：
            # 找到分子中的核心片段
            fragments = list(BRICS.BRICSDecompose(mol))
            
            if len(fragments) >= 2:
                # 第一个片段作为骨架，其他作为装饰
                scaffold = fragments[0]
                decorations = '.'.join(fragments[1:])
                return f"{scaffold}\t{decorations}"
        
        # 如果无法分解，跳过
        return None
        
    except:
        return None

def create_libinvent_dataset(input_file, output_file):
    """创建LibInvent格式的训练数据"""
    df = pd.read_csv(input_file)
    
    valid_entries = []
    skipped = 0
    
    print("正在分解分子...")
    for i, smiles in enumerate(df['Smiles']):
        if pd.isna(smiles):
            continue
        
        decomposed = decompose_molecule_for_libinvent(str(smiles))
        if decomposed:
            valid_entries.append(decomposed)
        else:
            skipped += 1
        
        if (i + 1) % 100 == 0:
            print(f"  处理进度: {i+1}/{len(df)}")
    
    # 保存
    with open(output_file, 'w') as f:
        for entry in valid_entries:
            f.write(entry + '\n')
    
    print(f"\nLibInvent数据集创建完成:")
    print(f"  成功分解: {len(valid_entries)}")
    print(f"  跳过: {skipped}")
    print(f"  保存到: {output_file}")
    
    return len(valid_entries)

# 创建LibInvent训练数据
count = create_libinvent_dataset(
    'data/NS3.csv',
    'data/libinvent_training.smi'
)
