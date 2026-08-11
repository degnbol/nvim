# fix_pybind_stubs: rdkit 2026.3.5
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['PatternFingerprintTautomerTarget', 'TautomerQuery', 'TautomerQueryCanSerialize']
class TautomerQuery(Boost.Python.instance):
    """
    The Tautomer Query Class.
      Creates a query that enables structure search accounting for matching of
      Tautomeric forms
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetModifiedAtoms(self) -> UnsignedLong_Vect:
        """
        """
    def GetModifiedBonds(self) -> UnsignedLong_Vect:
        """
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatch(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatches(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatches(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatchesWithTautomers(self, target: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
        """
    def GetTautomers(self) -> typing.Any:
        """
        """
    def GetTemplateMolecule(self) -> rdkit.Chem.Mol:
        """
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        """
    @typing.overload
    def IsSubstructOf(self, target: Mol, params: SubstructMatchParameters) -> bool:
        """
        """
    def PatternFingerprintTemplate(self, fingerprintSize: int = 2048) -> ExplicitBitVect:
        """
        """
    def ToBinary(self) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self, arg1: Mol) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, arg1: Mol, arg2: str) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, pickle: str) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
def PatternFingerprintTautomerTarget(target: Mol, fingerprintSize: int = 2048) -> ExplicitBitVect:
    """
    """
def TautomerQueryCanSerialize() -> bool:
    """
        Returns True if the TautomerQuery is serializable (requires that the RDKit was built with boost::serialization)
    
    """
