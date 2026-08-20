# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['FragCatGenerator', 'FragCatParams', 'FragCatalog', 'FragFPGenerator']
class FragCatGenerator(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddFragsFromMol(self, mol: Mol, fcat: FragCatalog) -> int:
        """
        """
    def __init__(self) -> None:
        ...
class FragCatParams(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetFuncGroup(self, fid: int) -> rdkit.Chem.Mol:
        """
        """
    def GetLowerFragLength(self) -> int:
        """
        """
    def GetNumFuncGroups(self) -> int:
        """
        """
    def GetTolerance(self) -> float:
        """
        """
    def GetTypeString(self) -> str:
        """
        """
    def GetUpperFragLength(self) -> int:
        """
        """
    def Serialize(self) -> str:
        """
        """
    def __init__(self, lLen: int, uLen: int, fgroupFilename: str, tol: float = 1e-08) -> None:
        ...
class FragCatalog(Boost.Python.instance):
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 128
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetBitDescription(self, idx: int) -> str:
        """
        """
    def GetBitDiscrims(self, idx: int) -> typing.Sequence[double]:
        """
        """
    def GetBitEntryId(self, idx: int) -> int:
        """
        """
    def GetBitFuncGroupIds(self, idx: int) -> typing.Sequence[int]:
        """
        """
    def GetBitOrder(self, idx: int) -> int:
        """
        """
    def GetCatalogParams(self) -> FragCatParams:
        """
        """
    def GetEntryBitId(self, idx: int) -> int:
        """
        """
    def GetEntryDescription(self, idx: int) -> str:
        """
        """
    def GetEntryDownIds(self, idx: int) -> typing.Sequence[int]:
        """
        """
    def GetEntryFuncGroupIds(self, idx: int) -> typing.Sequence[int]:
        """
        """
    def GetEntryOrder(self, idx: int) -> int:
        """
        """
    def GetFPLength(self) -> int:
        """
        """
    def GetNumEntries(self) -> int:
        """
        """
    def Serialize(self) -> str:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self, params: FragCatParams) -> None:
        ...
    @typing.overload
    def __init__(self, pickle: str) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
class FragFPGenerator(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetFPForMol(self, mol: Mol, fcat: FragCatalog) -> ExplicitBitVect:
        """
        """
    def __init__(self) -> None:
        ...
