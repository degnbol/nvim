# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing a C++ implementation of code for doing MMPA
"""
from __future__ import annotations
import typing
__all__: list[str] = ['FragmentMol']
@typing.overload
def FragmentMol(*args, **kwargs) -> tuple:
    """
        Does the fragmentation necessary for an MMPA analysis
    
    """
@typing.overload
def FragmentMol(*args, **kwargs) -> tuple:
    """
        Does the fragmentation necessary for an MMPA analysis
    
    """
@typing.overload
def FragmentMol(mol: Mol, bondsToCut: typing.Any, minCuts: int = 1, maxCuts: int = 3, resultsAsMols: bool = True) -> tuple:
    """
        Does the fragmentation necessary for an MMPA analysis
    
    """
