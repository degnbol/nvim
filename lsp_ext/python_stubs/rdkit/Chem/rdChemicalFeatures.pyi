# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing free chemical feature functionality
     These are features that are not associated with molecules. They are 
     typically derived from pharmacophores and site-maps.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['FreeChemicalFeature']
class FreeChemicalFeature(Boost.Python.instance):
    """
    Class to represent free chemical features.
        These chemical features are not associated with a molecule, though they can be matched 
        to molecular features
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 120
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetFamily(self) -> str:
        """
            Get the family of the feature
        
        """
    def GetId(self) -> int:
        """
            Get the id of the feature
        
        """
    def GetPos(self) -> Point3D:
        """
            Get the position of the feature
        
        """
    def GetType(self) -> str:
        """
            Get the sepcific type for the feature
        
        """
    def SetFamily(self, family: str) -> None:
        """
            Set the family of the feature
        
        """
    def SetId(self, id: int) -> None:
        """
            Set the id of the feature
        
        """
    def SetPos(self, loc: Point3D) -> None:
        """
            Set the feature position
        
        """
    def SetType(self, type: str) -> None:
        """
            Set the sepcific type for the feature
        
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self, pickle: str) -> None:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, family: str, type: str, loc: Point3D, id: int = -1) -> None:
        ...
    @typing.overload
    def __init__(self, family: str, loc: Point3D) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
