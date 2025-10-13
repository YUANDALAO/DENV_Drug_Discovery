print("开始正确的测试...")

try:
    import torch
    print("1. PyTorch导入成功")
    
    # 使用REINVENT的配置方式进行采样
    import tempfile
    import os
    
    # 创建采样配置文件
    sampling_config = """
run_type = "sampling"
device = "cuda:0"

[parameters]
model_file = "priors/denv_ultra_clean_model.model"
num_smiles = 10
"""
    
    # 写入临时配置文件
    with tempfile.NamedTemporaryFile(mode='w', suffix='.toml', delete=False) as f:
        f.write(sampling_config)
        config_file = f.name
    
    print("2. 配置文件创建成功")
    
    # 使用reinvent命令行工具进行采样
    import subprocess
    
    print("3. 开始生成分子...")
    result = subprocess.run(
        ['reinvent', config_file], 
        capture_output=True, 
        text=True,
        timeout=60
    )
    
    print("4. 生成完成")
    print("输出:", result.stdout)
    if result.stderr:
        print("错误:", result.stderr)
    
    # 清理临时文件
    os.unlink(config_file)
        
except Exception as e:
    print(f"错误: {e}")
    import traceback
    traceback.print_exc()
