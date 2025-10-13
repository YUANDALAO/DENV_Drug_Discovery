print("开始测试...")

try:
    import torch
    print("1. PyTorch导入成功")
    
    from reinvent.models import ModelAdapter
    print("2. ModelAdapter导入成功")
    
    print("3. 开始加载模型...")
    model_adapter = ModelAdapter.load_from_file("priors/denv_ultra_clean_model.model", "cpu")
    print("4. 模型加载成功")
    
    print("5. 开始生成分子...")
    sampled = model_adapter.sample(3)
    print("6. 生成完成")
    
    for i, smiles in enumerate(sampled, 1):
        print(f"   {i}. {smiles}")
        
except Exception as e:
    print(f"错误: {e}")
    import traceback
    traceback.print_exc()
