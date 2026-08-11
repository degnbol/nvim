# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing functions to enumerate stereoisomers of a molecule.  Chiral centers and double bonds will be enumerated if unassigned, or, if the appropriate option is set, if assigned.  Atropisomers will only be enumerated if assigned.  There is, as yet, no means of finding  unassigned atropisomers.
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['StereoEnumerationOptions', 'StereoisomerEnumerator']
class StereoEnumerationOptions(Boost.Python.instance):
    """
    EnumerateSteroisomers options.
    """
    __instance_size__: typing.ClassVar[int] = 48
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
    def maxIsomers(self) -> int:
        """The maximum number of isomers to yield.  If the number of possible isomers is greater than maxIsomers, a random subset will be yielded.  If 0, there is no maximum.  Since every additional stereocenter doubles the number of results (and execution time) it's important to keep an eye on this. (default: 0)"""
    @maxIsomers.setter
    def maxIsomers(self, value: int) -> None: ...
    @property
    def onlyStereoGroups(self) -> bool:
        """If true, only find stereoisomers that differ at the StereoGroups associated with the molecule.  Default=False. (default: False)"""
    @onlyStereoGroups.setter
    def onlyStereoGroups(self, value: bool) -> None: ...
    @property
    def onlyUnassigned(self) -> bool:
        """If true, stereocenters which have a specified stereochemistry will not be perturbed unless they are part of a relative stereo group.  Default=True. (default: True)"""
    @onlyUnassigned.setter
    def onlyUnassigned(self, value: bool) -> None: ...
    @property
    def randomSeed(self) -> int:
        """Seed for random number generator.  Default=-1 means no seed. (default: -1)"""
    @randomSeed.setter
    def randomSeed(self, value: int) -> None: ...
    @property
    def tryEmbedding(self) -> bool:
        """If true, the process attempts to generate a standard RDKit distance geometry conformation for the stereoisomer.  If this fails, we assume that the stereoisomer is non-physical and don't return it.  NOTE that this is computationally expensive and is just a heuristic that could result in stereoisomers being lost.  Default=False (default: False)"""
    @tryEmbedding.setter
    def tryEmbedding(self, value: bool) -> None: ...
    @property
    def unique(self) -> bool:
        """If true, only stereoisomers that differ in canonical CXSmiles will be returned.  Default=True. (default: True)"""
    @unique.setter
    def unique(self, value: bool) -> None: ...
class StereoisomerEnumerator(Boost.Python.instance):
    """
    Stereoisomer enumerator.
    """
    @staticmethod
    def GetStereoisomerCount(arg1: StereoisomerEnumerator) -> int:
        """
            Get the number of stereoisomers.
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, arg1: typing.Any, options: typing.Any = None, verbose: bool = True) -> None:
        ...
    def next(self) -> rdkit.Chem.Mol:
        """
            Get next isomer in the sequence, or None if at the end.
        
        """
