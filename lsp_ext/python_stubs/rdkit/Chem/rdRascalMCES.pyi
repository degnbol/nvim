# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing implementation of RASCAL Maximum Common Edge Substructure algorithm.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['FindMCES', 'RascalButinaCluster', 'RascalCluster', 'RascalClusterOptions', 'RascalOptions', 'RascalResult']
class RascalClusterOptions(Boost.Python.instance):
    """
    RASCAL Cluster Options.  Most of these pertain to RascalCluster calculations.  Only similarityCutoff is used by RascalButinaCluster.
    """
    __instance_size__: typing.ClassVar[int] = 80
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
    def a(self) -> float:
        """The penalty score for each unconnected component in the MCES. Default=0.05. (default: 0.05)"""
    @a.setter
    def a(self, value: float) -> None: ...
    @property
    def b(self) -> float:
        """The weight of matched bonds over matched atoms. Default=2. (default: 0.05)"""
    @b.setter
    def b(self, value: float) -> None: ...
    @property
    def clusterMergeSim(self) -> float:
        """Two clusters are merged if the fraction of molecules they have in common is greater than this.  Default=0.6. (default: 0.6)"""
    @clusterMergeSim.setter
    def clusterMergeSim(self, value: float) -> None: ...
    @property
    def maxNumFrags(self) -> int:
        """The maximum number of fragments allowed in the MCES for each pair of molecules. Default=2.  So that the MCES isn't a lot of small fragments scattered around the molecules giving an inflated estimate of similarity. (default: 2)"""
    @maxNumFrags.setter
    def maxNumFrags(self, value: int) -> None: ...
    @property
    def minFragSize(self) -> int:
        """The minimum number of atoms in a fragment for it to be included in the MCES.  Default=3. (default: 3)"""
    @minFragSize.setter
    def minFragSize(self, value: int) -> None: ...
    @property
    def minIntraClusterSim(self) -> float:
        """Two pairs of molecules are included in the same cluster if the similarity between their MCESs is greater than this.  Default=0.9. (default: 0.9)"""
    @minIntraClusterSim.setter
    def minIntraClusterSim(self, value: float) -> None: ...
    @property
    def numThreads(self) -> int:
        """Number of threads to use during clustering.  Default=-1 means all the hardware threads less one. (default: -1)"""
    @numThreads.setter
    def numThreads(self, value: int) -> None: ...
    @property
    def similarityCutoff(self) -> float:
        """Similarity cutoff for molecules to be in the same cluster.  Between 0.0 and 1.0, default=0.7. (default: 0.7)"""
    @similarityCutoff.setter
    def similarityCutoff(self, value: float) -> None: ...
class RascalOptions(Boost.Python.instance):
    """
    RASCAL Options
    """
    __instance_size__: typing.ClassVar[int] = 104
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
    def allBestMCESs(self) -> bool:
        """If True, reports all MCESs found of the same maximum size.  Default False means just report the first found. (default: False)"""
    @allBestMCESs.setter
    def allBestMCESs(self, value: bool) -> None: ...
    @property
    def completeAromaticRings(self) -> bool:
        """If True (default), partial aromatic rings won't be returned. (default: True)"""
    @completeAromaticRings.setter
    def completeAromaticRings(self, value: bool) -> None: ...
    @property
    def completeSmallestRings(self) -> bool:
        """If True (default is False), only complete rings present in both input molecule's RingInfo will be returned. Implies completeAromaticRings and ringMatchesRingOnly. (default: False)"""
    @completeSmallestRings.setter
    def completeSmallestRings(self, value: bool) -> None: ...
    @property
    def equivalentAtoms(self) -> str:
        """SMARTS strings defining atoms that shouldbe considered equivalent. e.g.[F,Cl,Br,I] so all halogens will match each other.Space-separated list allowing more than 1class of equivalent atoms. (default: '')"""
    @equivalentAtoms.setter
    def equivalentAtoms(self, value: str) -> None: ...
    @property
    def exactConnectionsMatch(self) -> bool:
        """If True (default is False), atoms will only match atoms if they have the same
         number of explicit connections.  E.g. the central atom of
         C(C)(C) won't match either atom in CC (default: False)"""
    @exactConnectionsMatch.setter
    def exactConnectionsMatch(self, value: bool) -> None: ...
    @property
    def ignoreAtomAromaticity(self) -> bool:
        """If True, matches atoms solely on atomic number.  If False, will treat aromatic and aliphatic atoms as different.  Default=True. (default: True)"""
    @ignoreAtomAromaticity.setter
    def ignoreAtomAromaticity(self, value: bool) -> None: ...
    @property
    def ignoreBondOrders(self) -> bool:
        """If True, will treat all bonds as the same, irrespective of order.  Default=False. (default: False)"""
    @ignoreBondOrders.setter
    def ignoreBondOrders(self, value: bool) -> None: ...
    @property
    def maxBestMCESs(self) -> int:
        """Some pathological cases produce huge numbers of equivalent solutions that can crash the program due to memory depletion.  This caps the number of such solutions to prevent this happening.  Default=10000. (default: 10000)"""
    @maxBestMCESs.setter
    def maxBestMCESs(self, value: int) -> None: ...
    @property
    def maxBondMatchPairs(self) -> int:
        """Too many matching bond (vertex) pairs can cause the process to run out of memory.  The default of 1000 is fairly safe.  Increase with caution, as memory use increases with the square of this number. (default: 1000)"""
    @maxBondMatchPairs.setter
    def maxBondMatchPairs(self, value: int) -> None: ...
    @property
    def maxFragSeparation(self) -> int:
        """Maximum number of bonds between fragments in the MCES for both to be reported.  Default -1 means no maximum.  If exceeded, the smaller fragment will be removed. (default: -1)"""
    @maxFragSeparation.setter
    def maxFragSeparation(self, value: int) -> None: ...
    @property
    def minCliqueSize(self) -> int:
        """Normally, the minimum clique size is specified via the similarityThreshold.  Sometimes it's more convenient to specify it directly.  If this is > 0, it will over-ride the similarityThreshold.  Note that this refers to the minimum number of BONDS in the MCES. Default=0. (default: 0)"""
    @minCliqueSize.setter
    def minCliqueSize(self, value: int) -> None: ...
    @property
    def minFragSize(self) -> int:
        """Imposes a minimum on the number of atoms in a fragment that may be part of the MCES.  Default -1 means no minimum. (default: -1)"""
    @minFragSize.setter
    def minFragSize(self, value: int) -> None: ...
    @property
    def returnEmptyMCES(self) -> bool:
        """If the estimated similarity between the 2 molecules doesn't meet the similarityThreshold, no results are returned.  If you want to know what the estimates were, set this to True, and examine the tier1Sim and tier2Sim properties of the result then returned. (default: False)"""
    @returnEmptyMCES.setter
    def returnEmptyMCES(self, value: bool) -> None: ...
    @property
    def ringMatchesRingOnly(self) -> bool:
        """If True (default is False), ring bonds won't match non-ring bonds. (default: False)"""
    @ringMatchesRingOnly.setter
    def ringMatchesRingOnly(self, value: bool) -> None: ...
    @property
    def similarityThreshold(self) -> float:
        """Threshold below which MCES won't be run.  Between 0.0 and 1.0, default=0.7. (default: 0.7)"""
    @similarityThreshold.setter
    def similarityThreshold(self, value: float) -> None: ...
    @property
    def singleLargestFrag(self) -> bool:
        """Return the just single largest fragment of the MCES. It is equivalent to running with allBestMCEs=True, finding the result with the largest largestFragmentSize, and calling its largestFragmentOnly method.  This option may not produce the largest possible single fragment that the molecules have in common. If you definitely want that you may be better off using rdFMCS. (default: False)"""
    @singleLargestFrag.setter
    def singleLargestFrag(self, value: bool) -> None: ...
    @property
    def timeout(self) -> int:
        """Maximum time (in seconds) to spend on an individual MCESs determination.  Default 60, -1 means no limit. (default: 60)"""
    @timeout.setter
    def timeout(self, value: int) -> None: ...
class RascalResult(Boost.Python.instance):
    """
    Used to return RASCAL MCES results.
    """
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def atomMatches(self) -> list:
        """
            Likewise for atoms.
        
        """
    def bondMatches(self) -> list:
        """
            A function returning a list of list of tuples, each inner list containing the matching bonds in the MCES as tuples of bond indices from mol1 and mol2
        
        """
    def largestFragmentOnly(self) -> None:
        """
            Function that cuts the MCES down to the single largest frag.  This cannot be undone.
        
        """
    @property
    def largestFragmentSize(*args, **kwargs):
        """
        Number of atoms in largest fragment.
        """
    @property
    def numFragments(*args, **kwargs):
        """
        Number of fragments in MCES.
        """
    @property
    def similarity(*args, **kwargs):
        """
        Johnson similarity between 2 molecules.
        """
    @property
    def smartsString(*args, **kwargs):
        """
        SMARTS string defining the MCES.
        """
    @property
    def tier1Sim(*args, **kwargs):
        """
        The tier 1 similarity estimate.
        """
    @property
    def tier2Sim(*args, **kwargs):
        """
        The tier 2 similarity estimate.
        """
    @property
    def timedOut(*args, **kwargs):
        """
        Whether it timed out.
        """
def FindMCES(mol1: Mol, mol2: Mol, opts: typing.Any = None) -> list:
    """
        Find one or more MCESs between the 2 molecules given.  Returns a list of RascalResult objects.- mol1- mol2 The two molecules for which to find the MCES- opts Optional RascalOptions object changing the default run mode.
    
    """
def RascalButinaCluster(mols: typing.Any, opts: typing.Any = None) -> list:
    """
        Use the RASCAL MCES similarity metric to do Butina clustering (Butina JCICS 39 747-750 (1999)).  Returns a list of lists of molecules, each inner list being a cluster.  The last cluster is all the molecules that didn't fit into another cluster (the singletons).- mols List of molecules to be clustered- opts Optional RascalOptions object changing the default run mode.
    
    """
def RascalCluster(mols: typing.Any, opts: typing.Any = None) -> list:
    """
        Use the RASCAL MCES similarity metric to do fuzzy clustering.  Returns a list of lists of molecules, each inner list being a cluster.  The last cluster is all the molecules that didn't fit into another cluster (the singletons).- mols List of molecules to be clustered- opts Optional RascalOptions object changing the default run mode.
    
    """
