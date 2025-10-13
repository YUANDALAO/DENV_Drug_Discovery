from rdkit import Chem

def prepare_final_libinvent_data():
    """使用LibInvent支持的正确格式"""
    
    # 骨架：使用简单的[*]标记连接点
    scaffolds = [
        "O=C([*])C1C([*])NC([*])C1",  # 吡咯烷 - 3个连接点
        "O=C([*])C1C([*])C([*])C1",   # 环丁烷 - 3个连接点
    ]
    
    # 装饰基团：使用[*]标记
    decorations = [
        "[*]C", "[*]CC", "[*]CCC", "[*]C(C)C",
        "[*]c1ccccc1", "[*]c1ccc(F)cc1", "[*]c1ccc(Cl)cc1",
        "[*]c1ccc(C)cc1", "[*]c1ccc(OC)cc1",
        "[*]c1ccncc1", "[*]c1ccoc1",
        "[*]O", "[*]OC", "[*]OCC",
        "[*]N", "[*]NC", "[*]NCC",
        "[*]F", "[*]Cl", "[*]Br",
        "[*]CF3", "[*]CN", "[*]C#N",
    ]
    
    # 验证
    valid_decorations = []
    for dec in decorations:
        mol = Chem.MolFromSmiles(dec)
        if mol:
            valid_decorations.append(dec)
    
    print(f"有效装饰基团: {len(valid_decorations)}")
    
    # LibInvent格式: scaffold  decoration1|decoration2|decoration3
    libinvent_entries = []
    
    for scaffold in scaffolds:
        for i, r1 in enumerate(valid_decorations[:10]):
            for r2 in valid_decorations[:8]:
                for r3 in valid_decorations[:6]:
                    decorations_str = f"{r1}|{r2}|{r3}"
                    entry = f"{scaffold}\t{decorations_str}"
                    libinvent_entries.append(entry)
                    
                    if len(libinvent_entries) >= 800:
                        break
                if len(libinvent_entries) >= 800:
                    break
            if len(libinvent_entries) >= 800:
                break
    
    with open('data/libinvent_training_final.smi', 'w') as f:
        for entry in libinvent_entries:
            f.write(entry + '\n')
    
    print(f"创建了 {len(libinvent_entries)} 个LibInvent训练样本")
    print("示例:")
    for entry in libinvent_entries[:3]:
        print(f"  {entry}")

prepare_final_libinvent_data()
