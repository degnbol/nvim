"""
Module containing functions to generate and work with reduced graphs
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['GenerateErGFingerprintForReducedGraph', 'GenerateMolExtendedReducedGraph', 'GetErGFingerprint']
def GenerateErGFingerprintForReducedGraph(mol: Mol, atomTypes: typing.Any = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> typing.Any:
    """
        Returns the ErG fingerprint vector for a reduced graph
    
    """
def GenerateMolExtendedReducedGraph(mol: Mol, atomTypes: typing.Any = 0) -> rdkit.Chem.Mol:
    """
        Returns the reduced graph for a molecule
    
    """
def GetErGFingerprint(mol: Mol, atomTypes: typing.Any = 0, fuzzIncrement: float = 0.3, minPath: int = 1, maxPath: int = 15) -> typing.Any:
    """
        Returns the ErG fingerprint vector for a molecule
    
    """
