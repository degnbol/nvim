# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing functions for interchange of molecules.
Note that this should be considered beta and that the format
  and API will very likely change in future releases.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['JSONParseParameters', 'JSONToMols', 'JSONWriteParameters', 'MolToJSON', 'MolsToJSON']
class JSONParseParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON parser
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def parseConformers(self) -> bool:
        """parse conformers in the JSON (default: True)"""
    @parseConformers.setter
    def parseConformers(self, value: bool) -> None: ...
    @property
    def parseProperties(self) -> bool:
        """parse molecular properties in the JSON (default: True)"""
    @parseProperties.setter
    def parseProperties(self, value: bool) -> None: ...
    @property
    def setAromaticBonds(self) -> bool:
        """set bond types to aromatic for bonds flagged aromatic (default: True)"""
    @setAromaticBonds.setter
    def setAromaticBonds(self, value: bool) -> None: ...
    @property
    def strictValenceCheck(self) -> bool:
        """be strict when checking atom valences (default: False)"""
    @strictValenceCheck.setter
    def strictValenceCheck(self, value: bool) -> None: ...
    @property
    def useHCounts(self) -> bool:
        """use atomic H counts from the JSON. You may want to set this to False when parsing queries. (default: True)"""
    @useHCounts.setter
    def useHCounts(self, value: bool) -> None: ...
class JSONWriteParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON writer
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def useRDKitExtensions(self) -> bool:
        """use RDKit extensions to the commonchem format (default: True)"""
    @useRDKitExtensions.setter
    def useRDKitExtensions(self, value: bool) -> None: ...
def JSONToMols(jsonBlock: str, params: typing.Any = None) -> tuple:
    """
        Convert JSON to a tuple of molecules
        
            ARGUMENTS:
              - jsonBlock: the molecule to work with
              - params: (optional) JSONParseParameters controlling the JSON parsing
            RETURNS:
              a tuple of Mols
        
    
    """
def MolToJSON(mol: Mol, params: typing.Any = None) -> str:
    """
        Convert a single molecule to JSON
        
            ARGUMENTS:
              - mol: the molecule to work with
            RETURNS:
              a string
        
    
    """
def MolsToJSON(mols: typing.Any, params: typing.Any = None) -> str:
    """
        Convert a set of molecules to JSON
        
            ARGUMENTS:
              - mols: the molecules to work with
            RETURNS:
              a string
        
    
    """
