# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
import numpy as numpy
__all__: list[str] = ['TanimotoSimilarity', 'numpy']
def TanimotoSimilarity(arr1, arr2):
    ...

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdReducedGraphs import GenerateErGFingerprintForReducedGraph as GenerateErGFingerprintForReducedGraph
from rdkit.Chem.rdReducedGraphs import GenerateMolExtendedReducedGraph as GenerateMolExtendedReducedGraph
from rdkit.Chem.rdReducedGraphs import GetErGFingerprint as GetErGFingerprint
