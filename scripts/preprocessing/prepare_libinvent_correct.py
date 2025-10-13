from rdkit import Chem

def prepare_correct_libinvent_data():
    """正确的LibInvent格式数据准备"""
    
    # 您的骨架（用单个连接点标记）
    scaffolds = [
        "O=C([*:1])C1C([*:2])NC([*:3])C1",  # 吡咯烷
        "O=C([*:1])C1C([*:2])C([*:3])C1",   # 环丁烷
    ]
    
    # 正确的装饰基团（每个只有一个连接点）
    decorations = [
        # 简单烷基
        "[*:1]C", "[*:1]CC", "[*:1]CCC", "[*:1]C(C)C",
        # 芳香基团
        "[*:1]c1ccccc1", "[*:1]c1ccc(F)cc1", "[*:1]c1ccc(Cl)cc1",
        "[*:1]c1ccc(C)cc1", "[*:1]c1ccc(OC)cc1",
        # 杂环
        "[*:1]c1ccncc1", "[*:1]c1ccoc1", "[*:1]c1ccsc1",
        # 含氧基团
        "[*:1]O", "[*:1]OC", "[*:1]OCC", "[*:1]C(=O)C",
        # 含氮基团
        "[*:1]N", "[*:1]NC", "[*:1]NCC", "[*:1]N(C)C",
        # 卤素
        "[*:1]F", "[*:1]Cl", "[*:1]Br",
        # 其他
        "[*:1]CF3", "[*:1]CN", "[*:1]C#N",
    ]
    
    # 验证装饰基团
    valid_decorations = []
    for dec in decorations:
        mol = Chem.MolFromSmiles(dec)
        if mol:
            valid_decorations.append(dec)
    
    print(f"有效装饰基团: {len(valid_decorations)}")
    
    # 创建LibInvent格式数据
    # 格式: scaffold  decoration1.decoration2.decoration3
    libinvent_entries = []
    
    for scaffold in scaffolds:
        # 为每个骨架创建不同的R基团组合
        for i, r1 in enumerate(valid_decorations[:10]):
            for r2 in valid_decorations[:8]:
                for r3 in valid_decorations[:6]:
                    # 组合装饰基团
                    decorations_str = f"{r1}|{r2}|{r3}"
                    entry = f"{scaffold}\t{decorations_str}"
                    libinvent_entries.append(entry)
                    
                    if len(libinvent_entries) >= 800:
                        break
                if len(libinvent_entries) >= 800:
                    break
            if len(libinvent_entries) >= 800:
                break
    
    # 保存
    with open('data/libinvent_training_correct.smi', 'w') as f:
        for entry in libinvent_entries:
            f.write(entry + '\n')
    
    print(f"\n创建了 {len(libinvent_entries)} 个LibInvent训练样本")
    print("前3个示例:")
    for i in range(min(3, len(libinvent_entries))):
        print(f"  {libinvent_entries[i]}")
    
    return len(libinvent_entries)

prepare_correct_libinvent_data()
