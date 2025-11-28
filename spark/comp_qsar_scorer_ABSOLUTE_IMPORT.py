"""
QSAR Scorer Plugin for REINVENT4
Predicts DENV NS2B-NS3 protease inhibitor activity (pIC50)
Optimized for NVIDIA DGX Spark

使用绝对导入以兼容REINVENT4插件系统
"""

import joblib
import numpy as np
from typing import List
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem

# 修改：使用绝对导入而不是相对导入
from reinvent.scoring.component_results import ComponentResults
from reinvent.runmodes.utils.add_tag import add_tag


@add_tag("__parameters")
@dataclass
class Parameters:
    """Parameters for QSAR scoring component"""
    model_path: str


@add_tag("__component")
class QSARScorer:
    """
    Custom QSAR scorer for DENV inhibitor activity prediction.
    Optimized for large-scale generation on DGX Spark.
    """
    
    def __init__(self, params: Parameters):
        self.params = params
        self.call_count = 0
        
        # Load model
        try:
            self.model = joblib.load(params.model_path)
            print(f"[QSAR Scorer] ✓ Model loaded from {params.model_path}")
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
            score = 5.0  # Default medium score for invalid molecules
            
            try:
                # Parse SMILES to molecule
                mol = Chem.MolFromSmiles(smi)
                
                if mol is not None:
                    # Generate Morgan fingerprint (radius=2, 4096 bits)
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
