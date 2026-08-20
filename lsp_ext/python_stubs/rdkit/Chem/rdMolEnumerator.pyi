# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing classes and functions for enumerating molecules
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['Enumerate', 'EnumeratorType', 'MolEnumeratorParams']
class EnumeratorType(Boost.Python.enum):
    LinkNode: typing.ClassVar[EnumeratorType]  # value = rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode
    PositionVariation: typing.ClassVar[EnumeratorType]  # value = rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation
    RepeatUnit: typing.ClassVar[EnumeratorType]  # value = rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'LinkNode': rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode, 'PositionVariation': rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation, 'RepeatUnit': rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdMolEnumerator.EnumeratorType.LinkNode, 1: rdkit.Chem.rdMolEnumerator.EnumeratorType.PositionVariation, 2: rdkit.Chem.rdMolEnumerator.EnumeratorType.RepeatUnit}
class MolEnumeratorParams(Boost.Python.instance):
    """
    Molecular enumerator parameters
    """
    __instance_size__: typing.ClassVar[int] = 64
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def SetEnumerationOperator(self, typ: EnumeratorType) -> None:
        """
            set the operator to be used for enumeration
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg1: EnumeratorType) -> typing.Any:
        ...
    @property
    def doRandom(self) -> bool:
        """do random enumeration (not yet implemented (default: False)"""
    @doRandom.setter
    def doRandom(self, value: bool) -> None: ...
    @property
    def maxToEnumerate(self) -> int:
        """maximum number of molecules to enumerate (default: 1000)"""
    @maxToEnumerate.setter
    def maxToEnumerate(self, value: int) -> None: ...
    @property
    def randomSeed(self) -> int:
        """seed for the random enumeration (not yet implemented (default: -1)"""
    @randomSeed.setter
    def randomSeed(self, value: int) -> None: ...
    @property
    def sanitize(self) -> bool:
        """sanitize molecules after enumeration (default: False)"""
    @sanitize.setter
    def sanitize(self, value: bool) -> None: ...
@typing.overload
def Enumerate(mol: Mol, maxPerOperation: int = 0) -> rdkit.Chem.MolBundle:
    """
        do an enumeration and return a MolBundle.
          If maxPerOperation is >0 that will be used as the maximum number of molecules which 
            can be returned by any given operation.
        Limitations:
          - the current implementation does not support molecules which include both
            SRUs and LINKNODEs
          - Overlapping SRUs, i.e. where one monomer is contained within another, are
            not supported
    
    """
@typing.overload
def Enumerate(mol: Mol, enumParams: MolEnumeratorParams) -> rdkit.Chem.MolBundle:
    """
        do an enumeration for the supplied parameter type and return a MolBundle
        Limitations:
          - Overlapping SRUs, i.e. where one monomer is contained within another, are
            not supported
    
    """
