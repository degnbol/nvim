"""
Module containing implementation of SynthonSpace search of Synthon-based chemical libraries such as Enamine REAL.  NOTE: This functionality is experimental and the API and/or results may change in future releases.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['ConvertTextToDBFile', 'FormattedIntegerString', 'SubstructureResult', 'SynthonSpace', 'SynthonSpaceSearchParams']
class SubstructureResult(Boost.Python.instance):
    """
    Used to return results of SynthonSpace searches.
    """
    @staticmethod
    def GetCancelled(arg1: SubstructureResult) -> bool:
        """
            Returns whether the search was cancelled or not.
        
        """
    @staticmethod
    def GetMaxNumResults(arg1: SubstructureResult) -> int:
        """
            The upper bound on number of results possible.  There may be fewer than this in practice for several reasons such as duplicate reagent sets being removed or the final product not matching the query even though the synthons suggested they would.
        
        """
    @staticmethod
    def GetTimedOut(arg1: SubstructureResult) -> bool:
        """
            Returns whether the search timed out or not.
        
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
    def GetHitMolecules(self) -> list:
        """
            A function returning hits from the search
        
        """
class SynthonSpace(Boost.Python.instance):
    """
    SynthonSpaceSearch object.
    """
    __instance_size__: typing.ClassVar[int] = 152
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def BuildSynthonFingerprints(self, fingerprintGenerator: typing.Any) -> None:
        """
            Build the synthon fingerprints ready for similarity searching.  This is done automatically when the first similarity search is done, but if converting a text file to binary format it might need to be done explicitly.
        
        """
    def FingerprintSearch(self, query: Mol, fingerprintGenerator: typing.Any, params: typing.Any = None) -> SubstructureResult:
        """
            Does a fingerprint search in the SynthonSpace using the FingerprintGenerator passed in.
        
        """
    def FingerprintSearchIncremental(self, query: Mol, fingerprintGenerator: typing.Any, callback: typing.Any, params: typing.Any = None) -> None:
        """
            Does a fingerprint search in the SynthonSpace using the FingerprintGenerator passed in, returning results the callback.
        
        """
    def GetNumProducts(self) -> int:
        """
            Returns number of products in the SynthonSpace, with multiple counting of any duplicates.
        
        """
    def GetNumReactions(self) -> int:
        """
            Returns number of reactions in the SynthonSpace.
        
        """
    def GetSynthonFingerprintType(self) -> str:
        """
            Returns the information string for the fingerprint generator used to create this space.
        
        """
    def RascalSearch(self, query: Mol, rascalOptions: typing.Any, params: typing.Any = None) -> SubstructureResult:
        """
            Does a search using the Rascal similarity score.  The similarity threshold used is provided by rascalOptions, and the one in params is ignored.
        
        """
    def RascalSearchIncremental(self, query: Mol, rascalOptions: typing.Any, callback: typing.Any, params: typing.Any = None) -> None:
        """
            Does a search using the Rascal similarity score.  The similarity threshold used is provided by rascalOptions, and the one in params is ignored.  Returns results iteratively in the callback.
        
        """
    def ReadDBFile(self, inFile: str, numThreads: int = 1) -> None:
        """
            Reads binary database file.  Takes optional number of threads,default=1.
        
        """
    def ReadTextFile(self, inFile: str) -> None:
        """
            Reads text file of the sort used by ChemSpace/Enamine.
        
        """
    @typing.overload
    def SubstructureSearch(self, query: Mol, substructMatchParams: typing.Any = None, params: typing.Any = None) -> SubstructureResult:
        """
            Does a substructure search in the SynthonSpace.
        
        """
    @typing.overload
    def SubstructureSearch(self, query: typing.Any, substructMatchParams: typing.Any = None, params: typing.Any = None) -> SubstructureResult:
        """
            Does a substructure search in the SynthonSpace using an extended query.
        
        """
    def SubstructureSearchIncremental(self, query: Mol, callback: typing.Any, substructMatchParams: typing.Any = None, params: typing.Any = None) -> None:
        """
            Does a substructure search in the SynthonSpace returning results in the callback.
        
        """
    def Summarise(self) -> None:
        """
            Writes a summary of the SynthonSpace to stdout.
        
        """
    def WriteDBFile(self, outFile: str) -> None:
        """
            Writes binary database file.
        
        """
    def WriteEnumeratedFile(self, outFile: str) -> None:
        """
            Writes enumerated library to file.
        
        """
    def __init__(self) -> None:
        ...
class SynthonSpaceSearchParams(Boost.Python.instance):
    """
    SynthonSpaceSearch parameters.
    """
    __instance_size__: typing.ClassVar[int] = 144
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
    def approxSimilarityAdjuster(self) -> float:
        """The fingerprint search uses an approximate similarity method before building a product and doing a final check.  The similarityCutoff is reduced by this value for the approximate check.  A lower value will give faster run times at the risk of missing some hits.  The value you use should have a positive correlation with your FOMO.  The default of 0.1 is appropriate for Morgan fingerprints.  With RDKit fingerprints, 0.05 is adequate, and higher than that has been seen to produce long run times. (default: 0.1)"""
    @approxSimilarityAdjuster.setter
    def approxSimilarityAdjuster(self, value: float) -> None: ...
    @property
    def buildHits(self) -> bool:
        """If false, reports the maximum number of hits that the search could produce, but doesn't return them. (default: True)"""
    @buildHits.setter
    def buildHits(self, value: bool) -> None: ...
    @property
    def fragSimilarityAdjuster(self) -> float:
        """Similarities of fragments are generally low due to low bit densities.  For the fragment matching, reduce the similarity cutoff off by this amount.  Default=0.1. (default: 0.1)"""
    @fragSimilarityAdjuster.setter
    def fragSimilarityAdjuster(self, value: float) -> None: ...
    @property
    def hitStart(self) -> int:
        """The sequence number of the hit to start from.  So that you can return the next N hits of a search having already obtained N-1.  Default=0 (default: 0)"""
    @hitStart.setter
    def hitStart(self, value: int) -> None: ...
    @property
    def maxHitChiralAtoms(self) -> int:
        """Maximum number of chiral atoms in a hit.  Default=-1 means no maximum. (default: -1)"""
    @maxHitChiralAtoms.setter
    def maxHitChiralAtoms(self, value: int) -> None: ...
    @property
    def maxHitHeavyAtoms(self) -> int:
        """Maximum number of heavy atoms in a hit.  Default=-1 means no maximum. (default: -1)"""
    @maxHitHeavyAtoms.setter
    def maxHitHeavyAtoms(self, value: int) -> None: ...
    @property
    def maxHitMolWt(self) -> float:
        """Maximum molecular weight for a hit.  Default=0.0 mean no maximum. (default: 0.0)"""
    @maxHitMolWt.setter
    def maxHitMolWt(self, value: float) -> None: ...
    @property
    def maxHits(self) -> int:
        """The maximum number of hits to return.  Default=1000.Use -1 for no maximum. (default: 1000)"""
    @maxHits.setter
    def maxHits(self, value: int) -> None: ...
    @property
    def maxNumFrags(self) -> int:
        """The maximum number of fragments the query can be broken into.  Big molecules will create huge numbers of fragments that may cause excessive memory use.  If the number of fragments hits this number, fragmentation stops and the search results will likely be incomplete.  Default=100000. (default: 100000)"""
    @maxNumFrags.setter
    def maxNumFrags(self, value: int) -> None: ...
    @property
    def minHitChiralAtoms(self) -> int:
        """Minimum number of chiral atoms in a hit.  Default=0. (default: 0)"""
    @minHitChiralAtoms.setter
    def minHitChiralAtoms(self, value: int) -> None: ...
    @property
    def minHitHeavyAtoms(self) -> int:
        """Minimum number of heavy atoms in a hit.  Default=0. (default: 0)"""
    @minHitHeavyAtoms.setter
    def minHitHeavyAtoms(self, value: int) -> None: ...
    @property
    def minHitMolWt(self) -> float:
        """Minimum molecular weight for a hit.  Default=0.0. (default: 0.0)"""
    @minHitMolWt.setter
    def minHitMolWt(self, value: float) -> None: ...
    @property
    def numRandomSweeps(self) -> int:
        """The random sampling doesn't always produce the required number of hits in 1 go.  This parameter controls how many loops it makes to try and get the hits before giving up.  Default=10. (default: 10)"""
    @numRandomSweeps.setter
    def numRandomSweeps(self, value: int) -> None: ...
    @property
    def numThreads(self) -> int:
        """The number of threads to use for search.  If > 0, will use that number.  If <= 0, will use the number of hardware threads plus this number.  So if the number of hardware threads is 8, and numThreads is -1, it will use 7 threads.  Default=1. (default: 1)"""
    @numThreads.setter
    def numThreads(self, value: int) -> None: ...
    @property
    def randomSample(self) -> bool:
        """If True, returns a random sample of the hits, up to maxHits in number.  Default=False. (default: False)"""
    @randomSample.setter
    def randomSample(self, value: bool) -> None: ...
    @property
    def randomSeed(self) -> int:
        """If using randomSample, this seeds the random number generator so as to give reproducible results.  Default=-1 means use a random seed. (default: -1)"""
    @randomSeed.setter
    def randomSeed(self, value: int) -> None: ...
    @property
    def similarityCutoff(self) -> float:
        """Similarity cutoff for returning hits by fingerprint similarity.  At present the fp is hard-coded to be Morgan, bits, radius=2.  Default=0.5. (default: 0.5)"""
    @similarityCutoff.setter
    def similarityCutoff(self, value: float) -> None: ...
    @property
    def timeOut(self) -> int:
        """Time limit for search, in seconds.  Default is 600s, 0 means no timeout.  Requires an integer (default: 600)"""
    @timeOut.setter
    def timeOut(self, value: int) -> None: ...
    @property
    def toTryChunkSize(self) -> int:
        """Process possible hits using the given chunk size (default: 2500000)"""
    @toTryChunkSize.setter
    def toTryChunkSize(self, value: int) -> None: ...
def ConvertTextToDBFile(inFilename: str, outFilename: str, fpGen: typing.Any = None) -> None:
    """
        Convert the text file into the binary DB file in our format.  Assumes that all synthons from a reaction are contiguous in the input file.  This uses a lot less memory than using ReadTextFile() followed by  WriteDBFile().- inFilename the name of the text file- outFilename the name of the binary file- optional fingerprint generator
    
    """
def FormattedIntegerString(value: int) -> str:
    """
        Format an integer with spaces every 3 digits for ease of reading
    
    """
