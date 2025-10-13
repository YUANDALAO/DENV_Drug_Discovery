"""
QSAR Scorer Plugin for REINVENT4
Predicts DENV NS2B-NS3 protease inhibitor activity (pIC50)
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
            print(f"[QSAR Scorer] Processing {len(smiles_list)} SMILES strings")
        
        for i, smi in enumerate(smiles_list):
            score = 5.0  # Default medium score
            
            try:
                # Parse SMILES to molecule
                mol = Chem.MolFromSmiles(smi)
                
                if mol is not None:
                    # Generate Morgan fingerprint
                    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=4096)
                    fp_array = np.array(list(fp), dtype=float).reshape(1, -1)
                    
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
            print(f"[QSAR Scorer] Predictions: min={scores_array.min():.2f}, max={scores_array.max():.2f}, mean={scores_array.mean():.2f}")
        
        return ComponentResults(scores=[scores_array])
