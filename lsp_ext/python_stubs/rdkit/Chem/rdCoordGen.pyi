"""
Module containing interface to the CoordGen library.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['AddCoords', 'CoordGenParams', 'SetDefaultTemplateFileDir']
class CoordGenParams(Boost.Python.instance):
    """
    Parameters controlling coordinate generation
    """
    __instance_size__: typing.ClassVar[int] = 112
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def SetCoordMap(self, coordMap: dict) -> None:
        """
            expects a dictionary of Point2D objects with template coordinates
        
        """
    def SetTemplateMol(self, templ: Mol) -> None:
        """
            sets a molecule to be used as the template
        
        """
    def __init__(self) -> None:
        ...
    @property
    def coordgenScaling(self) -> float:
        """scaling factor for a single bond (default: 50.0)"""
    @coordgenScaling.setter
    def coordgenScaling(self, value: float) -> None: ...
    @property
    def dbg_useConstrained(self) -> bool:
        """for debugging use (default: True)"""
    @dbg_useConstrained.setter
    def dbg_useConstrained(self, value: bool) -> None: ...
    @property
    def dbg_useFixed(self) -> bool:
        """for debugging use (default: False)"""
    @dbg_useFixed.setter
    def dbg_useFixed(self, value: bool) -> None: ...
    @property
    def minimizerPrecision(self) -> float:
        """controls sketcher precision (default: 0.009999999776482582)"""
    @minimizerPrecision.setter
    def minimizerPrecision(self, value: float) -> None: ...
    @property
    def sketcherBestPrecision(self) -> float:
        """highest quality (and slowest) precision setting (default: 3.0)"""
    @property
    def sketcherCoarsePrecision(self) -> float:
        """"coarse" (fastest) precision setting, produces good-quality coordinates most of the time, this is the default setting for the RDKit (default: 0.009999999776482582)"""
    @property
    def sketcherQuickPrecision(self) -> float:
        """faster precision setting (default: 0.20000000298023224)"""
    @property
    def sketcherStandardPrecision(self) -> float:
        """standard quality precision setting, the default for the coordgen project (default: 1.0)"""
    @property
    def templateFileDir(self) -> str:
        """directory containing the templates.mae file (default: '')"""
    @templateFileDir.setter
    def templateFileDir(self, value: str) -> None: ...
    @property
    def treatNonterminalBondsToMetalAsZOBs(self) -> bool:
        """default: False"""
    @treatNonterminalBondsToMetalAsZOBs.setter
    def treatNonterminalBondsToMetalAsZOBs(self, value: bool) -> None: ...
def AddCoords(mol: Mol, params: typing.Any = None) -> None:
    """
        Add 2D coordinates.
        ARGUMENTS:
           - mol: molecule to modify
           - params: (optional) parameters controlling the coordinate generation
        
        
    
    """
def SetDefaultTemplateFileDir(dir: str) -> None:
    """
    """
