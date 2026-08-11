from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['CreateMolCatalog', 'MolCatalog', 'MolCatalogEntry']
class MolCatalog(Boost.Python.instance):
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 128
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddEdge(self, id1: int, id2: int) -> None:
        """
        """
    def AddEntry(self, entry: MolCatalogEntry) -> int:
        """
        """
    def GetBitDescription(self, idx: int) -> str:
        """
        """
    def GetBitEntryId(self, idx: int) -> int:
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
    def __init__(self, pickle: str) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
class MolCatalogEntry(Boost.Python.instance):
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 88
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetDescription(self) -> str:
        """
        """
    def GetMol(self) -> rdkit.Chem.Mol:
        """
        """
    def GetOrder(self) -> int:
        """
        """
    def SetDescription(self, val: str) -> None:
        """
        """
    def SetMol(self, mol: Mol) -> None:
        """
        """
    def SetOrder(self, order: int) -> None:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, pickle: str) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
def CreateMolCatalog() -> MolCatalog:
    """
    """
