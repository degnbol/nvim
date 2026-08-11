# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing functions for working with groups of molecules
"""
from __future__ import annotations
import typing
__all__: list[str] = ['GetFingerprintsForMolsInFile', 'SupplierOptions']
class SupplierOptions(Boost.Python.instance):
    """
    Supplier Options
    """
    __instance_size__: typing.ClassVar[int] = 112
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
    def confId2D(self) -> int:
        """used for TDT files (default: -1)"""
    @confId2D.setter
    def confId2D(self, value: int) -> None: ...
    @property
    def confId3D(self) -> int:
        """used for TDT files (default: 0)"""
    @confId3D.setter
    def confId3D(self, value: int) -> None: ...
    @property
    def delimiter(self) -> str:
        """used for SMILES files (default: '\t')"""
    @delimiter.setter
    def delimiter(self, value: str) -> None: ...
    @property
    def nameColumn(self) -> int:
        """used for SMILES files (default: 1)"""
    @nameColumn.setter
    def nameColumn(self, value: int) -> None: ...
    @property
    def nameRecord(self) -> str:
        """used for TDT files (default: '')"""
    @nameRecord.setter
    def nameRecord(self, value: str) -> None: ...
    @property
    def numThreads(self) -> int:
        """the number of threads to use while working (default: 0)"""
    @numThreads.setter
    def numThreads(self, value: int) -> None: ...
    @property
    def removeHs(self) -> bool:
        """default: True"""
    @removeHs.setter
    def removeHs(self, value: bool) -> None: ...
    @property
    def sanitize(self) -> bool:
        """default: True"""
    @sanitize.setter
    def sanitize(self, value: bool) -> None: ...
    @property
    def smilesColumn(self) -> int:
        """used for SMILES files (default: 0)"""
    @smilesColumn.setter
    def smilesColumn(self, value: int) -> None: ...
    @property
    def strictParsing(self) -> bool:
        """default: True"""
    @strictParsing.setter
    def strictParsing(self, value: bool) -> None: ...
    @property
    def titleLine(self) -> bool:
        """used for SMILES files (default: True)"""
    @titleLine.setter
    def titleLine(self, value: bool) -> None: ...
@typing.overload
def GetFingerprintsForMolsInFile(filename: str, generator: typing.Any = None, options: SupplierOptions = ...) -> tuple:
    """
        returns the fingerprints for the molecules in a file (32 bit version)
    
    """
@typing.overload
def GetFingerprintsForMolsInFile(filename: str, generator: typing.Any = None, options: SupplierOptions = ...) -> tuple:
    """
        returns the fingerprints for the molecules in a file (64 bit version)
    
    """
