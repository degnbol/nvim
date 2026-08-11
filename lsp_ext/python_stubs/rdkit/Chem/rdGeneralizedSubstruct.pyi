# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing functions for generalized substructure searching
"""
from __future__ import annotations
import typing
__all__: list[str] = ['CreateExtendedQueryMol', 'ExtendedQueryMol', 'MolGetSubstructMatch', 'MolGetSubstructMatches', 'MolHasSubstructMatch', 'PatternFingerprintTarget']
class ExtendedQueryMol(Boost.Python.instance):
    """
    Extended query molecule for use in generalized substructure searching.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def InitFromBinary(self, pkl: str) -> None:
        """
        """
    def InitFromJSON(self, text: str) -> None:
        """
        """
    def PatternFingerprintQuery(self, fingerprintSize: int = 2048) -> typing.Any:
        """
            C++ signature :
                std::__1::unique_ptr<ExplicitBitVect, std::__1::default_delete<ExplicitBitVect>> PatternFingerprintQuery(RDKit::GeneralizedSubstruct::ExtendedQueryMol {lvalue} [,unsigned int=2048])
        """
    def ToBinary(self) -> typing.Any:
        """
        """
    def ToJSON(self) -> str:
        """
        """
    def __init__(self, text: str, isJSON: bool = False) -> None:
        ...
def CreateExtendedQueryMol(mol: Mol, doEnumeration: bool = True, doTautomers: bool = True, adjustQueryProperties: bool = False, adjustQueryParameters: AdjustQueryParameters = None) -> ExtendedQueryMol:
    """
        Creates an ExtendedQueryMol from the input molecule
        
          This takes a query molecule and, conceptually, performs the following steps to
          produce an ExtendedQueryMol:
        
            1. Enumerates features like Link Nodes and SRUs
            2. Converts everything into TautomerQueries
            3. Runs adjustQueryProperties()
        
          Each step is optional
        
    
    """
def MolGetSubstructMatch(mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> typing.Any:
    """
        returns first match (if any) of a molecule to a generalized substructure query
    
    """
def MolGetSubstructMatches(mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> typing.Any:
    """
        returns all matches (if any) of a molecule to a generalized substructure query
    
    """
def MolHasSubstructMatch(mol: Mol, query: ExtendedQueryMol, params: SubstructMatchParameters = None) -> bool:
    """
        determines whether or not a molecule is a match to a generalized substructure query
    
    """
def PatternFingerprintTarget(target: Mol, fingerprintSize: int = 2048) -> typing.Any:
    """
        Creates a pattern fingerprint for a target molecule that is compatible with an extended query
    
        C++ signature :
            std::__1::unique_ptr<ExplicitBitVect, std::__1::default_delete<ExplicitBitVect>> PatternFingerprintTarget(RDKit::ROMol [,unsigned int=2048])
    """
