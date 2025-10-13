from rdkit import Chem

def prepare_correct_format():
    """使用LibInvent实际支持的格式"""
    
    # 骨架中使用[*]
    scaffolds = [
        "O=C([*])C1C([*])NC([*])C1",
        "O=C([*])C1C([*])C([*])C1",
    ]
    
    # 装饰基团中只用普通原子，不用*标记
    # 装饰基团会自动连接到骨架的[*]位置
    decorations = [
        "C", "CC", "CCC", "C(C)C",  # 烷基
        "c1ccccc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1", "c1ccc(C)cc1",  # 芳香基
        "c1ccncc1", "c1ccoc1",  # 杂环
        "O", "OC", "OCC",  # 含氧
        "N", "NC", "NCC",  # 含氮
        "F", "Cl", "Br",  # 卤素
        "C(F)(F)F", "CN", "C#N",  # 其他
    ]
    
    valid_decorations = []
    for dec in decorations:
        mol = Chem.MolFromSmiles(dec)
        if mol:
            valid_decorations.append(dec)
    
    print(f"有效装饰基团: {len(valid_decorations)}")
    
    # LibInvent格式: scaffold  decoration1|decoration2|decoration3
    entries = []
    for scaffold in scaffolds:
        for r1 in valid_decorations[:10]:
            for r2 in valid_decorations[:8]:
                for r3 in valid_decorations[:6]:
                    entry = f"{scaffold}\t{r1}|{r2}|{r3}"
                    entries.append(entry)
                    if len(entries) >= 800:
                        break
                if len(entries) >= 800:
                    break
            if len(entries) >= 800:
                break
    
    with open('data/libinvent_final_format.smi', 'w') as f:
        for entry in entries:
            f.write(entry + '\n')
    
    print(f"创建了 {len(entries)} 个训练样本")
    for i in range(3):
        print(f"  {entries[i]}")

prepare_correct_format()
