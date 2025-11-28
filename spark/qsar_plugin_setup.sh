#!/bin/bash
# QSAR插件安装脚本 for NVIDIA DGX Spark

echo "=========================================="
echo "REINVENT4 QSAR插件安装向导"
echo "=========================================="
echo ""

# 1. 设置路径
REINVENT_DIR="$HOME/aidd/REINVENT4"
PLUGIN_DIR="$REINVENT_DIR/reinvent/reinvent_plugins/components"

echo "步骤 1: 检查REINVENT4安装..."
if [ ! -d "$REINVENT_DIR" ]; then
    echo "❌ 错误: 找不到REINVENT4目录: $REINVENT_DIR"
    exit 1
fi
echo "✅ REINVENT4目录存在"

# 2. 创建插件目录（如果不存在）
echo ""
echo "步骤 2: 检查插件目录..."
if [ ! -d "$PLUGIN_DIR" ]; then
    echo "创建插件目录: $PLUGIN_DIR"
    mkdir -p "$PLUGIN_DIR"
fi
echo "✅ 插件目录: $PLUGIN_DIR"

# 3. 复制QSAR scorer插件
echo ""
echo "步骤 3: 安装QSAR scorer插件..."
cat > "$PLUGIN_DIR/comp_qsar_scorer.py" << 'PLUGIN_EOF'
"""
QSAR Scorer Plugin for REINVENT4
Predicts DENV NS2B-NS3 protease inhibitor activity (pIC50)
Optimized for NVIDIA DGX Spark
"""

import joblib
import numpy as np
from typing import List, Union
from pydantic.dataclasses import dataclass
from pydantic import field_validator

from rdkit import Chem
from rdkit.Chem import AllChem

from .component_results import ComponentResults
from .add_tag import add_tag


@add_tag("__parameters")
@dataclass
class Parameters:
    """Parameters for QSAR scoring component"""
    model_path: Union[str, List[str]]
    
    @field_validator('model_path')
    @classmethod
    def convert_list_to_str(cls, v):
        """Handle both string and list inputs from TOML parser"""
        if isinstance(v, list):
            return v[0] if v else ""
        return v


@add_tag("__component")
class QSARScorer:
    """
    Custom QSAR scorer for DENV inhibitor activity prediction.
    Optimized for large-scale generation on DGX Spark.
    """
    
    def __init__(self, params: Parameters):
        self.params = params
        self.call_count = 0
        
        # Extract model path
        model_path = params.model_path
        if isinstance(model_path, list):
            model_path = model_path[0]
        
        # Load model
        try:
            self.model = joblib.load(model_path)
            print(f"[QSAR Scorer] ✓ Model loaded from {model_path}")
        except Exception as e:
            raise RuntimeError(f"[QSAR Scorer] ✗ Failed to load model: {e}")
    
    def __call__(self, smiles_list: List[str]) -> ComponentResults:
        """
        Calculate QSAR scores for given SMILES strings.
        
        Args:
            smiles_list: List of SMILES strings
            
        Returns:
            ComponentResults with predicted pIC50 values
        """
        self.call_count += 1
        scores = []
        
        if self.call_count == 1:
            print(f"[QSAR Scorer] Processing batch of {len(smiles_list)} SMILES")
        
        for i, smi in enumerate(smiles_list):
            score = 5.0  # Default medium score
            
            try:
                # Parse SMILES to molecule
                mol = Chem.MolFromSmiles(smi)
                
                if mol is not None:
                    # Generate Morgan fingerprint (optimized)
                    fp = AllChem.GetMorganFingerprintAsBitVect(
                        mol, radius=2, nBits=4096
                    )
                    fp_array = np.array(list(fp), dtype=np.float32).reshape(1, -1)
                    
                    # Predict pIC50
                    predicted_pic50 = float(self.model.predict(fp_array)[0])
                    score = predicted_pic50
                    
                    # Log first few predictions
                    if self.call_count == 1 and i < 3:
                        print(f"[QSAR Scorer]   {smi[:50]:50s} -> pIC50={score:.4f}")
                        
            except Exception as e:
                if self.call_count == 1 and i < 2:
                    print(f"[QSAR Scorer]   ERROR on '{smi[:30]}': {e}")
            
            scores.append(score)
        
        scores_array = np.array(scores, dtype=np.float32)
        
        if self.call_count == 1:
            print(f"[QSAR Scorer] Stats: min={scores_array.min():.2f}, "
                  f"max={scores_array.max():.2f}, mean={scores_array.mean():.2f}")
        
        return ComponentResults(scores=[scores_array])
PLUGIN_EOF

echo "✅ QSAR插件已安装: $PLUGIN_DIR/comp_qsar_scorer.py"

# 4. 验证安装
echo ""
echo "步骤 4: 验证安装..."
if [ -f "$PLUGIN_DIR/comp_qsar_scorer.py" ]; then
    echo "✅ 插件文件存在"
    # 检查文件大小
    SIZE=$(wc -c < "$PLUGIN_DIR/comp_qsar_scorer.py")
    echo "   文件大小: $SIZE bytes"
else
    echo "❌ 警告: 插件文件未找到"
fi

# 5. 创建模型目录
echo ""
echo "步骤 5: 创建模型目录..."
MODEL_DIR="$REINVENT_DIR/models"
mkdir -p "$MODEL_DIR"
echo "✅ 模型目录: $MODEL_DIR"
echo "   请将 random_forest_champion.joblib 复制到此目录"

# 6. 创建实验目录
echo ""
echo "步骤 6: 创建实验目录..."
EXP_DIR="$REINVENT_DIR/experiments/runs/spark_run1"
mkdir -p "$EXP_DIR"
echo "✅ 实验目录: $EXP_DIR"

echo ""
echo "=========================================="
echo "安装完成！"
echo "=========================================="
echo ""
echo "接下来需要手动完成："
echo "1. 复制 random_forest_champion.joblib 到 $MODEL_DIR/"
echo "2. 确保 priors/ 目录包含:"
echo "   - libinvent.prior"
echo "   - denv_libinvent_model.model"
echo "3. 确保 data/ 目录包含:"
echo "   - pyrrolidine_dual_aryl.smi"
echo ""
echo "然后运行配置文件验证插件加载"
echo ""

