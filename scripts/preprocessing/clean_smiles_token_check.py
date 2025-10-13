import pandas as pd
from rdkit import Chem
import re

def extract_tokens_from_smiles(smiles):
    """从SMILES中提取所有token"""
    # 更精确的token匹配模式
    pattern = r'\[[^\]]+\]|Br|Cl|[BCNOSPFIHbcnops()]|\d|[=#\-+\\\/]|%\d{2}|\^|\$'
    tokens = re.findall(pattern, smiles)
    return set(tokens)

def is_token_compatible(smiles):
    """检查SMILES是否只包含REINVENT支持的token"""
    # REINVENT精确的允许token集合
    allowed_tokens = {
        '4', 'S', 's', '[S+]', '2', '=', '(', '6', ')', '7', '[nH]', 'Br', '5', 
        'n', '-', '%10', '^', 'O', '[N+]', '3', 'C', 'Cl', '[n+]', 'F', '9', 
        'o', '1', '[O-]', '[N-]', '8', '#', 'N', '$', 'c'
    }
    
    # 提取SMILES中的所有token
    smiles_tokens = extract_tokens_from_smiles(smiles)
    
    # 检查是否有不支持的token
    unsupported = smiles_tokens - allowed_tokens
    
    if unsupported:
        return False, unsupported
    return True, set()

def ultra_clean_smiles(input_file, output_file):
    """超严格的SMILES清理，基于token验证"""
    df = pd.read_csv(input_file)
    
    print(f"开始ultra清理，原始数据: {len(df)} 个化合物")
    
    valid_smiles = []
    rejected_examples = []
    
    for i, smiles in enumerate(df['Smiles']):
        if pd.isna(smiles):
            continue
            
        smiles_str = str(smiles).strip()
        
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is not None:
                canonical_smiles = Chem.MolToSmiles(mol)
                
                # 基本过滤
                if '.' in canonical_smiles or '|' in canonical_smiles:
                    continue
                
                if mol.GetNumAtoms() > 60:
                    continue
                
                # Token兼容性检查
                is_compatible, unsupported_tokens = is_token_compatible(canonical_smiles)
                
                if is_compatible:
                    valid_smiles.append(canonical_smiles)
                else:
                    if len(rejected_examples) < 3:
                        rejected_examples.append((canonical_smiles[:50], unsupported_tokens))
                        print(f"拒绝: {canonical_smiles[:50]}... 不支持token: {unsupported_tokens}")
        except Exception as e:
            continue
    
    # 去重
    unique_smiles = list(dict.fromkeys(valid_smiles))
    
    # 保存
    with open(output_file, 'w') as f:
        for smiles in unique_smiles:
            f.write(f"{smiles}\n")
    
    print(f"\nUltra清理结果:")
    print(f"  Token兼容的SMILES: {len(unique_smiles)}")
    print(f"  被拒绝: {len(df) - len(unique_smiles)}")
    print(f"  最终成功率: {len(unique_smiles)/len(df)*100:.1f}%")
    print(f"  保存到: {output_file}")

if __name__ == "__main__":
    ultra_clean_smiles("data/NS3.csv", "data/denv_smiles_ultra.tsv")
