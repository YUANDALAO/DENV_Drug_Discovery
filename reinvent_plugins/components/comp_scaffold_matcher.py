"""
Custom Scaffold Matcher - 修复版
奖励包含目标骨架的分子
"""

from __future__ import annotations
__all__ = ["ScaffoldMatcher"]

from typing import List
from rdkit import Chem
import numpy as np
from pydantic.dataclasses import dataclass

from .component_results import ComponentResults
from reinvent_plugins.mol_cache import molcache
from .add_tag import add_tag


@add_tag("__parameters")
@dataclass
class Parameters:
    smarts: List[str]
    use_chirality: List[bool]


@add_tag("__component")  # 移除"penalty"标记
class ScaffoldMatcher:
    def __init__(self, params: Parameters):
        self.patterns = []
        self.use_chirality = params.use_chirality

        for smarts in params.smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                self.patterns.append(pattern)

        if not self.patterns:
            raise ValueError("No valid SMARTS patterns")

    @molcache
    def __call__(self, mols: List[Chem.Mol]) -> ComponentResults:
        scores = []

        for pattern, use_chirality in zip(self.patterns, self.use_chirality):
            # 修复：直接返回0.0或1.0
            match = np.array([
                1.0 if mol.HasSubstructMatch(pattern, useChirality=use_chirality) else 0.0
                for mol in mols
            ], dtype=np.float32)
            
            scores.append(match)

        return ComponentResults(scores)
