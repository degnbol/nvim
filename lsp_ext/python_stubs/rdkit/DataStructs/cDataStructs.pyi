# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing an assortment of functionality for basic data structures.

At the moment the data structures defined are:
  Bit Vector classes (for storing signatures, fingerprints and the like:
    - ExplicitBitVect: class for relatively small (10s of thousands of bits) or
                       dense bit vectors.
    - SparseBitVect:   class for large, sparse bit vectors
  DiscreteValueVect:   class for storing vectors of integers
  SparseIntVect:       class for storing sparse vectors of integers
"""
from __future__ import annotations
import typing
__all__: list[str] = ['AllBitSimilarity', 'AllProbeBitsMatch', 'AsymmetricSimilarity', 'AsymmetricSimilarityNeighbors', 'AsymmetricSimilarityNeighbors_sparse', 'BitVectToBinaryText', 'BitVectToFPSText', 'BitVectToText', 'BraunBlanquetSimilarity', 'BraunBlanquetSimilarityNeighbors', 'BraunBlanquetSimilarityNeighbors_sparse', 'BulkAllBitSimilarity', 'BulkAsymmetricSimilarity', 'BulkBraunBlanquetSimilarity', 'BulkCosineSimilarity', 'BulkDiceSimilarity', 'BulkKulczynskiSimilarity', 'BulkMcConnaugheySimilarity', 'BulkOnBitSimilarity', 'BulkRogotGoldbergSimilarity', 'BulkRusselSimilarity', 'BulkSokalSimilarity', 'BulkTanimotoSimilarity', 'BulkTverskySimilarity', 'ComputeL1Norm', 'ConvertToExplicit', 'ConvertToNumpyArray', 'CosineSimilarity', 'CosineSimilarityNeighbors', 'CosineSimilarityNeighbors_sparse', 'CreateFromBinaryText', 'CreateFromBitString', 'CreateFromFPSText', 'DiceSimilarity', 'DiceSimilarityNeighbors', 'DiceSimilarityNeighbors_sparse', 'DiscreteValueType', 'DiscreteValueVect', 'EIGHTBITVALUE', 'ExplicitBitVect', 'FOURBITVALUE', 'FPBReader', 'FoldFingerprint', 'InitFromDaylightString', 'IntSparseIntVect', 'KulczynskiSimilarity', 'KulczynskiSimilarityNeighbors', 'KulczynskiSimilarityNeighbors_sparse', 'LongSparseIntVect', 'McConnaugheySimilarity', 'McConnaugheySimilarityNeighbors', 'McConnaugheySimilarityNeighbors_sparse', 'MultiFPBReader', 'NumBitsInCommon', 'ONEBITVALUE', 'OffBitProjSimilarity', 'OffBitsInCommon', 'OnBitProjSimilarity', 'OnBitSimilarity', 'OnBitsInCommon', 'RealValueVect', 'RogotGoldbergSimilarity', 'RogotGoldbergSimilarityNeighbors', 'RogotGoldbergSimilarityNeighbors_sparse', 'RusselSimilarity', 'RusselSimilarityNeighbors', 'RusselSimilarityNeighbors_sparse', 'SIXTEENBITVALUE', 'SokalSimilarity', 'SokalSimilarityNeighbors', 'SokalSimilarityNeighbors_sparse', 'SparseBitVect', 'TWOBITVALUE', 'TanimotoSimilarity', 'TanimotoSimilarityNeighbors', 'TanimotoSimilarityNeighbors_sparse', 'TverskySimilarity', 'UIntSparseIntVect', 'ULongSparseIntVect']
class DiscreteValueType(Boost.Python.enum):
    EIGHTBITVALUE: typing.ClassVar[DiscreteValueType]  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE
    FOURBITVALUE: typing.ClassVar[DiscreteValueType]  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE
    ONEBITVALUE: typing.ClassVar[DiscreteValueType]  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE
    SIXTEENBITVALUE: typing.ClassVar[DiscreteValueType]  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE
    TWOBITVALUE: typing.ClassVar[DiscreteValueType]  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'ONEBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE, 'TWOBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, 'FOURBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE, 'EIGHTBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE, 'SIXTEENBITVALUE': rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE}
    values: typing.ClassVar[dict]  # value = {0: rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE, 1: rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, 2: rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE, 3: rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE, 4: rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE}
class DiscreteValueVect(Boost.Python.instance):
    """
    A container class for storing unsigned integer
    values within a particular range.
    
    The length of the vector and type of its elements (determines the maximum value
    that can be stored) are both set at construction time.
    
    As you would expect, _DiscreteValueVects_ support a set of binary operations
    so you can do things like:
      dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry
      dvv3 = dvv1 | dvv2  the result contains the largest value in each entry
      dvv1 += dvv2     values are truncated when necessary
      dvv3 = dvv1 + dvv2    values are truncated when necessary
      dvv1 -= dvv3    would-be negative values are set to zero
      dvv3 = dvv1 - dvv2    would-be negative values are set to zero
    
    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 64
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetTotalVal(self) -> int:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    def GetValueType(self) -> DiscreteValueType:
        """
            Get the type of value stored in the vector
        
        """
    def __add__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
    def __and__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, i: int) -> int:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, valType: DiscreteValueType, length: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    def __isub__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
            Get the number of entries in the vector
        
        """
    def __or__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
    def __setitem__(self, i: int, val: int) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: DiscreteValueVect) -> typing.Any:
        """
        """
class ExplicitBitVect(Boost.Python.instance):
    """
    A class to store explicit bit vectors.
    
    This class is most useful for situations where the size of the vector
    is relatively small (tens of thousands or smaller).
    
    For larger vectors, use the _SparseBitVect_ class instead.
    
    As you would expect, _ExplicitBitVects_ support a set of binary operations
    so you can do things like:
      bv3 = bv1 & bv2  (bitwise and)
      bv3 = bv1 | bv2  (bitwise or)
      bv3 = bv1 ^ bv2  (bitwise xor)
      bv3 = ~bv1       (bitwise negation)
    
    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def ToBitString(*args, **kwargs):
        """
        
        BitVectToText( (SparseBitVect)bv1) -> str :
        
            C++ signature :
                std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> BitVectToText(SparseBitVect)
        
        BitVectToText( (ExplicitBitVect)bv1) -> str :
            Returns a string of zeros and ones representing the bit vector.
        
            C++ signature :
                std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> BitVectToText(ExplicitBitVect)
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def FromBase64(self, inD: str) -> None:
        """
            Initializes the vector from a base64 encoded binary string.
            
        
        """
    def GetBit(self, which: int) -> bool:
        """
            Returns the value of a bit.
            
        
        """
    def GetNumBits(self) -> int:
        """
            Returns the number of bits in the vector (the vector's size).
            
        
        """
    def GetNumOffBits(self) -> int:
        """
            Returns the number of off bits.
            
        
        """
    def GetNumOnBits(self) -> int:
        """
            Returns the number of on bits.
            
        
        """
    def GetOnBits(self) -> typing.Sequence[int]:
        """
            Returns a tuple containing IDs of the on bits.
            
        
        """
    def SetBit(self, which: int) -> bool:
        """
            Turns on a particular bit.  Returns the original state of the bit.
            
        
        """
    def SetBitsFromList(self, onBitList: typing.Any) -> None:
        """
            Turns on a set of bits.  The argument should be a tuple or list of bit ids.
            
        
        """
    def ToBase64(self) -> str:
        """
            Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
            
        
        """
    def ToBinary(self) -> typing.Any:
        """
            Returns an internal binary representation of the vector.
            
        
        """
    def ToList(self) -> list:
        """
            Return the Bitvector as a python list (faster than list(vect))
        
        """
    def UnSetBit(self, which: int) -> bool:
        """
            Turns off a particular bit.  Returns the original state of the bit.
            
        
        """
    def UnSetBitsFromList(self, offBitList: typing.Any) -> None:
        """
            Turns off a set of bits.  The argument should be a tuple or list of bit ids.
            
        
        """
    def __add__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    def __and__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    def __eq__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, which: int) -> int:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, size: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    @typing.overload
    def __init__(self, size: int, bitsSet: bool) -> None:
        ...
    def __invert__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __ne__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    def __or__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
    def __setitem__(self, which: int, val: int) -> int:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __xor__(self, other: ExplicitBitVect) -> typing.Any:
        """
        """
class FPBReader(Boost.Python.instance):
    """
    A class for reading and searching FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetBytes(self, which: int) -> typing.Any:
        """
            returns a particular fingerprint as bytes
        
        """
    def GetContainingNeighbors(self, bv: str) -> tuple:
        """
            returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
        
        """
    def GetFP(self, idx: int) -> ExplicitBitVect:
        """
            returns a particular fingerprint as an ExplicitBitVect
        
        """
    def GetId(self, idx: int) -> str:
        """
            returns the id of a particular fingerprint
        
        """
    def GetNumBits(self) -> int:
        """
            returns the number of bits in a fingerprint
        
        """
    def GetTanimoto(self, which: int, bytes: str) -> float:
        """
            return the tanimoto similarity of a particular fingerprint to the bytes provided
        
        """
    def GetTanimotoNeighbors(self, bv: str, threshold: float = 0.7) -> tuple:
        """
            returns tanimoto similarities to and indices of all neighbors above the specified threshold
        
        """
    def GetTversky(self, which: int, bytes: str, ca: float, cb: float) -> float:
        """
            return the Tverksy similarity of a particular fingerprint to the bytes provided
        
        """
    def GetTverskyNeighbors(self, bv: str, ca: float, cb: float, threshold: float = 0.7) -> tuple:
        """
            returns Tversky similarities to and indices of all neighbors above the specified threshold
        
        """
    def Init(self) -> None:
        """
            Read the fingerprints from the file. This can take a while.
            
        
        """
    def __getitem__(self, which: int) -> tuple:
        """
        """
    def __init__(self, filename: str, lazy: bool = False) -> None:
        ...
    def __len__(self) -> int:
        """
        """
class IntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    
    The length of the vector is set at construction time.
    
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry
    
    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetLength(self) -> int:
        """
            Returns the length of the vector
        
        """
    def GetNonzeroElements(self) -> dict:
        """
            returns a dictionary of the nonzero elements
        
        """
    def GetTotalVal(self, useAbs: bool = False) -> int:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    def ToBinary(self) -> typing.Any:
        """
            returns a binary (pickle) representation of the vector
        
        """
    def ToList(self) -> list:
        """
            Return the SparseIntVect as a python list
        
        """
    def UpdateFromSequence(self, seq: typing.Any) -> None:
        """
            update the vector based on the values in the list or tuple
        
        """
    def __add__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    def __and__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    def __eq__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, item: int) -> int:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __iadd__(self, other: int) -> typing.Any:
        """
        """
    def __idiv__(self, other: int) -> typing.Any:
        """
        """
    def __imul__(self, other: int) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    @typing.overload
    def __isub__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __isub__(self, other: int) -> typing.Any:
        """
        """
    def __ne__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    def __or__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
    def __setitem__(self, item: int, value: int) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: IntSparseIntVect) -> typing.Any:
        """
        """
class LongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    
    The length of the vector is set at construction time.
    
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry
    
    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetLength(self) -> int:
        """
            Returns the length of the vector
        
        """
    def GetNonzeroElements(self) -> dict:
        """
            returns a dictionary of the nonzero elements
        
        """
    def GetTotalVal(self, useAbs: bool = False) -> int:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    def ToBinary(self) -> typing.Any:
        """
            returns a binary (pickle) representation of the vector
        
        """
    def ToList(self) -> list:
        """
            Return the SparseIntVect as a python list
        
        """
    def UpdateFromSequence(self, seq: typing.Any) -> None:
        """
            update the vector based on the values in the list or tuple
        
        """
    def __add__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    def __and__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    def __eq__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, item: int) -> int:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __iadd__(self, other: int) -> typing.Any:
        """
        """
    def __idiv__(self, other: int) -> typing.Any:
        """
        """
    def __imul__(self, other: int) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    @typing.overload
    def __isub__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __isub__(self, other: int) -> typing.Any:
        """
        """
    def __ne__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    def __or__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
    def __setitem__(self, item: int, value: int) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: LongSparseIntVect) -> typing.Any:
        """
        """
class MultiFPBReader(Boost.Python.instance):
    """
    A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddReader(self, rdr: FPBReader) -> int:
        """
            adds an FPBReader to our set of readers
        
        """
    def GetContainingNeighbors(self, bv: str, numThreads: int = 1) -> tuple:
        """
            returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
        
        """
    def GetNumBits(self) -> int:
        """
            returns the number of bits in a fingerprint
        
        """
    def GetReader(self, which: int) -> FPBReader:
        """
            returns one of our readers
        
        """
    def GetTanimotoNeighbors(self, bv: str, threshold: float = 0.7, numThreads: int = 1) -> tuple:
        """
            returns tanimoto similarities to and indices of all neighbors above the specified threshold
        
        """
    def GetTverskyNeighbors(self, bv: str, ca: float, cb: float, threshold: float = 0.7, numThreads: int = 1) -> tuple:
        """
            returns Tversky similarities to and indices of all neighbors above the specified threshold
        
        """
    def Init(self) -> None:
        """
            Call Init() on each of our children. This can take a while.
            
        
        """
    def __init__(self, initOnSearch: bool = False) -> None:
        ...
    def __len__(self) -> int:
        """
        """
class RealValueVect(Boost.Python.instance):
    """
    A container class for storing real
    values.
    
    The length of the vector is set at construction time.
    
    As you would expect, _RealValueVects_ support a set of binary operations
    so you can do things like:
      rvv3 = rvv1 & rvv2  the result contains the smallest value in each entry
      rvv3 = rvv1 | rvv2  the result contains the largest value in each entry
      rvv1 += rvv2     
      rvv3 = rvv1 + rvv2    
      rvv1 -= rvv3    
      rvv3 = rvv1 - rvv2    
    
    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 56
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def GetTotalVal(arg1: RealValueVect) -> float:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __add__(self, other: RealValueVect) -> typing.Any:
        """
        """
    def __and__(self, other: RealValueVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, item: int) -> float:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: RealValueVect) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, arg1: str) -> None:
        ...
    def __isub__(self, other: RealValueVect) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
            Get the number of entries in the vector
        
        """
    def __or__(self, other: RealValueVect) -> typing.Any:
        """
        """
    def __setitem__(self, item: int, value: float) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: RealValueVect) -> typing.Any:
        """
        """
class SparseBitVect(Boost.Python.instance):
    """
    A class to store sparse bit vectors.
    
    This class is most useful for situations where the size of the vector
    is large and relatively few bits are set
    
    For smaller or denser vectors, the _ExplicitBitVect_ class is much faster.
    
    As you would expect, _SparseBitVects_ support a set of binary operations
    so you can do things like:
      bv3 = bv1 & bv2  (bitwise and)
      bv3 = bv1 | bv2  (bitwise or)
      bv3 = bv1 ^ bv2  (bitwise xor)
      bv3 = ~bv1       (bitwise negation) NOTE: this operation is likely
                        to be VERY slow and inefficient.
    
    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def ToBitString(*args, **kwargs):
        """
        
        BitVectToText( (SparseBitVect)bv1) -> str :
        
            C++ signature :
                std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> BitVectToText(SparseBitVect)
        
        BitVectToText( (ExplicitBitVect)bv1) -> str :
            Returns a string of zeros and ones representing the bit vector.
        
            C++ signature :
                std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char>> BitVectToText(ExplicitBitVect)
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def FromBase64(self, inD: str) -> None:
        """
            Initializes the vector from a base64 encoded binary string.
            
        
        """
    def GetBit(self, which: int) -> bool:
        """
            Returns the value of a bit.
            
        
        """
    def GetNumBits(self) -> int:
        """
            Returns the number of bits in the vector (the vector's size).
            
        
        """
    def GetNumOffBits(self) -> int:
        """
            Returns the number of off bits.
            
        
        """
    def GetNumOnBits(self) -> int:
        """
            Returns the number of on bits.
            
        
        """
    def GetOnBits(self) -> typing.Sequence[int]:
        """
            Returns a tuple containing IDs of the on bits.
            
        
        """
    def SetBit(self, which: int) -> bool:
        """
            Turns on a particular bit.  Returns the original state of the bit.
            
        
        """
    def SetBitsFromList(self, onBitList: typing.Any) -> None:
        """
            Turns on a set of bits.  The argument should be a tuple or list of bit ids.
            
        
        """
    def ToBase64(self) -> str:
        """
            Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
            
        
        """
    def ToBinary(self) -> typing.Any:
        """
            Returns an internal binary representation of the vector.
            
        
        """
    def ToList(self) -> list:
        """
            Return the BitVector as a python list
        
        """
    def UnSetBit(self, which: int) -> bool:
        """
            Turns off a particular bit.  Returns the original state of the bit.
            
        
        """
    def UnSetBitsFromList(self, offBitList: typing.Any) -> None:
        """
            Turns off a set of bits.  The argument should be a tuple or list of bit ids.
            
        
        """
    def __and__(self, other: SparseBitVect) -> typing.Any:
        """
        """
    def __eq__(self, other: SparseBitVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, which: int) -> int:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self, size: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    def __invert__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __ne__(self, other: SparseBitVect) -> typing.Any:
        """
        """
    def __or__(self, other: SparseBitVect) -> typing.Any:
        """
        """
    def __setitem__(self, which: int, val: int) -> int:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __xor__(self, other: SparseBitVect) -> typing.Any:
        """
        """
class UIntSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    
    The length of the vector is set at construction time.
    
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry
    
    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetLength(self) -> int:
        """
            Returns the length of the vector
        
        """
    def GetNonzeroElements(self) -> dict:
        """
            returns a dictionary of the nonzero elements
        
        """
    def GetTotalVal(self, useAbs: bool = False) -> int:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    def ToBinary(self) -> typing.Any:
        """
            returns a binary (pickle) representation of the vector
        
        """
    def ToList(self) -> list:
        """
            Return the SparseIntVect as a python list
        
        """
    def UpdateFromSequence(self, seq: typing.Any) -> None:
        """
            update the vector based on the values in the list or tuple
        
        """
    def __add__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    def __and__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    def __eq__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, item: int) -> int:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __iadd__(self, other: int) -> typing.Any:
        """
        """
    def __idiv__(self, other: int) -> typing.Any:
        """
        """
    def __imul__(self, other: int) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    @typing.overload
    def __isub__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __isub__(self, other: int) -> typing.Any:
        """
        """
    def __ne__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    def __or__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
    def __setitem__(self, item: int, value: int) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: UIntSparseIntVect) -> typing.Any:
        """
        """
class ULongSparseIntVect(Boost.Python.instance):
    """
    A container class for storing integer
    values within a particular range.
    
    The length of the vector is set at construction time.
    
    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry
    
    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetLength(self) -> int:
        """
            Returns the length of the vector
        
        """
    def GetNonzeroElements(self) -> dict:
        """
            returns a dictionary of the nonzero elements
        
        """
    def GetTotalVal(self, useAbs: bool = False) -> int:
        """
            Get the sum of the values in the vector, basically L1 norm
        
        """
    def ToBinary(self) -> typing.Any:
        """
            returns a binary (pickle) representation of the vector
        
        """
    def ToList(self) -> list:
        """
            Return the SparseIntVect as a python list
        
        """
    def UpdateFromSequence(self, seq: typing.Any) -> None:
        """
            update the vector based on the values in the list or tuple
        
        """
    def __add__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    def __and__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    def __eq__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, item: int) -> int:
        """
            Get the value at a specified location
        
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __iadd__(self, other: int) -> typing.Any:
        """
        """
    def __idiv__(self, other: int) -> typing.Any:
        """
        """
    def __imul__(self, other: int) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    @typing.overload
    def __isub__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    @typing.overload
    def __isub__(self, other: int) -> typing.Any:
        """
        """
    def __ne__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    def __or__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
    def __setitem__(self, item: int, value: int) -> None:
        """
            Set the value at a specified location
        
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: ULongSparseIntVect) -> typing.Any:
        """
        """
@typing.overload
def AllBitSimilarity(v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    """
@typing.overload
def AllBitSimilarity(v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    """
        (B(bv1) - B(bv1^bv2)) / B(bv1)
    
    """
@typing.overload
def AllProbeBitsMatch(probe: SparseBitVect, ref: SparseBitVect) -> bool:
    """
    """
@typing.overload
def AllProbeBitsMatch(probe: ExplicitBitVect, ref: ExplicitBitVect) -> bool:
    """
    """
@typing.overload
def AllProbeBitsMatch(probe: SparseBitVect, ref: str) -> bool:
    """
    """
@typing.overload
def AllProbeBitsMatch(probe: ExplicitBitVect, ref: str) -> bool:
    """
        Returns True if all bits in the first argument match all bits in the 
          vector defined by the pickle in the second argument.
        
    
    """
@typing.overload
def AsymmetricSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def AsymmetricSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / min(B(bv1),B(bv2))
    
    """
@typing.overload
def AsymmetricSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def AsymmetricSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / min(B(bv1),B(bv2))
    
    """
def AsymmetricSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / min(B(bv1),B(bv2))
    
    """
def AsymmetricSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / min(B(bv1),B(bv2))
    
    """
@typing.overload
def BitVectToBinaryText(bv: SparseBitVect) -> typing.Any:
    """
    """
@typing.overload
def BitVectToBinaryText(bv: ExplicitBitVect) -> typing.Any:
    """
        Returns a binary string (byte array) representing the bit vector.
    
    """
@typing.overload
def BitVectToFPSText(bv1: SparseBitVect) -> str:
    """
    """
@typing.overload
def BitVectToFPSText(bv1: ExplicitBitVect) -> str:
    """
        Returns an FPS string representing the bit vector.
    
    """
@typing.overload
def BitVectToText(bv1: SparseBitVect) -> str:
    """
    """
@typing.overload
def BitVectToText(bv1: ExplicitBitVect) -> str:
    """
        Returns a string of zeros and ones representing the bit vector.
    
    """
@typing.overload
def BraunBlanquetSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def BraunBlanquetSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / max(B(bv1),B(bv2))
    
    """
@typing.overload
def BraunBlanquetSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def BraunBlanquetSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / max(B(bv1),B(bv2))
    
    """
def BraunBlanquetSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / max(B(bv1),B(bv2))
    
    """
def BraunBlanquetSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / max(B(bv1),B(bv2))
    
    """
@typing.overload
def BulkAllBitSimilarity(v1: ExplicitBitVect, v2: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkAllBitSimilarity(v1: ExplicitBitVect, v2: typing.Any, returnDistance: bool = 0) -> list:
    """
        (B(bv1) - B(bv1^bv2)) / B(bv1)
    
    """
@typing.overload
def BulkAsymmetricSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkAsymmetricSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / min(B(bv1),B(bv2))
    
    """
@typing.overload
def BulkBraunBlanquetSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkBraunBlanquetSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / max(B(bv1),B(bv2))
    
    """
@typing.overload
def BulkCosineSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkCosineSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))
    
    """
@typing.overload
def BulkDiceSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkDiceSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        2*B(bv1&bv2) / (B(bv1) + B(bv2))
    
    """
@typing.overload
def BulkDiceSimilarity(v1: IntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Dice similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkDiceSimilarity(v1: LongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Dice similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkDiceSimilarity(v1: UIntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Dice similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkDiceSimilarity(v1: ULongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Dice similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkKulczynskiSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkKulczynskiSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))
    
    """
@typing.overload
def BulkMcConnaugheySimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkMcConnaugheySimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))
    
    """
@typing.overload
def BulkOnBitSimilarity(v1: ExplicitBitVect, v2: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkOnBitSimilarity(v1: ExplicitBitVect, v2: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / B(bv1|bv2)
    
    """
@typing.overload
def BulkRogotGoldbergSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkRogotGoldbergSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def BulkRusselSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkRusselSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def BulkSokalSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkSokalSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))
    
    """
@typing.overload
def BulkTanimotoSimilarity(bv1: SparseBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkTanimotoSimilarity(bv1: ExplicitBitVect, bvList: typing.Any, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))
    
    """
@typing.overload
def BulkTanimotoSimilarity(v1: IntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Tanimoto similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTanimotoSimilarity(v1: LongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Tanimoto similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTanimotoSimilarity(v1: UIntSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Tanimoto similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTanimotoSimilarity(v1: ULongSparseIntVect, v2: list, returnDistance: bool = False) -> list:
    """
        return the Tanimoto similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTverskySimilarity(bv1: SparseBitVect, bvList: typing.Any, a: float, b: float, returnDistance: bool = 0) -> list:
    """
    """
@typing.overload
def BulkTverskySimilarity(bv1: ExplicitBitVect, bvList: typing.Any, a: float, b: float, returnDistance: bool = 0) -> list:
    """
        B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)
    
    """
@typing.overload
def BulkTverskySimilarity(v1: IntSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    """
        return the Tversky similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTverskySimilarity(v1: LongSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    """
        return the Tversky similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTverskySimilarity(v1: UIntSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    """
        return the Tversky similarities between one vector and a sequence of others
    
    """
@typing.overload
def BulkTverskySimilarity(v1: ULongSparseIntVect, v2: list, a: float, b: float, returnDistance: bool = False) -> list:
    """
        return the Tversky similarities between one vector and a sequence of others
    
    """
@typing.overload
def ComputeL1Norm(v1: DiscreteValueVect, v2: DiscreteValueVect) -> int:
    """
        Compute the distance between two discrete vector values
        
    
    """
@typing.overload
def ComputeL1Norm(arg1: RealValueVect, arg2: RealValueVect) -> float:
    """
        Compute the distance between two real vector values
        
    
    """
def ConvertToExplicit(sbv: SparseBitVect) -> ExplicitBitVect:
    """
        Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector
    
    """
@typing.overload
def ConvertToNumpyArray(bv: ExplicitBitVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(bv: DiscreteValueVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(bv: IntSparseIntVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(bv: LongSparseIntVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(bv: UIntSparseIntVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(bv: ULongSparseIntVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def ConvertToNumpyArray(rvv: RealValueVect, destArray: typing.Any) -> None:
    """
    """
@typing.overload
def CosineSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def CosineSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))
    
    """
@typing.overload
def CosineSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def CosineSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))
    
    """
def CosineSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))
    
    """
def CosineSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / sqrt(B(bv1) * B(bv2))
    
    """
def CreateFromBinaryText(fps: str) -> ExplicitBitVect:
    """
        Creates an ExplicitBitVect from a binary string (byte array).
    
    """
def CreateFromBitString(bits: str) -> ExplicitBitVect:
    """
        Creates an ExplicitBitVect from a bit string (string of 0s and 1s).
    
    """
def CreateFromFPSText(fps: str) -> ExplicitBitVect:
    """
        Creates an ExplicitBitVect from an FPS string.
    
    """
@typing.overload
def DiceSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def DiceSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        2*B(bv1&bv2) / (B(bv1) + B(bv2))
    
    """
@typing.overload
def DiceSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def DiceSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        2*B(bv1&bv2) / (B(bv1) + B(bv2))
    
    """
@typing.overload
def DiceSimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Dice similarity between two vectors
    
    """
@typing.overload
def DiceSimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Dice similarity between two vectors
    
    """
@typing.overload
def DiceSimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Dice similarity between two vectors
    
    """
@typing.overload
def DiceSimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Dice similarity between two vectors
    
    """
def DiceSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        2*B(bv1&bv2) / (B(bv1) + B(bv2))
    
    """
def DiceSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        2*B(bv1&bv2) / (B(bv1) + B(bv2))
    
    """
@typing.overload
def FoldFingerprint(bv: SparseBitVect, foldFactor: int = 2) -> SparseBitVect:
    """
    """
@typing.overload
def FoldFingerprint(bv: ExplicitBitVect, foldFactor: int = 2) -> ExplicitBitVect:
    """
        Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.
    
    """
@typing.overload
def InitFromDaylightString(sbv: SparseBitVect, s: str) -> None:
    """
    """
@typing.overload
def InitFromDaylightString(sbv: ExplicitBitVect, s: str) -> None:
    """
        Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.
        
           **Arguments**
             - bv: either a _SparseBitVect_ or an _ExplicitBitVect_
             - txt: a string with the Daylight encoding (this is the text that
                    the Daylight tools put in the FP field of a TDT)
        
        
    
    """
@typing.overload
def KulczynskiSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def KulczynskiSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))
    
    """
@typing.overload
def KulczynskiSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def KulczynskiSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))
    
    """
def KulczynskiSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))
    
    """
def KulczynskiSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))
    
    """
@typing.overload
def McConnaugheySimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def McConnaugheySimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))
    
    """
@typing.overload
def McConnaugheySimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def McConnaugheySimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))
    
    """
def McConnaugheySimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))
    
    """
def McConnaugheySimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        (B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))
    
    """
@typing.overload
def NumBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> int:
    """
    """
@typing.overload
def NumBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> int:
    """
        Returns the total number of bits in common between the two bit vectors
    
    """
@typing.overload
def OffBitProjSimilarity(bv1: SparseBitVect, bv2: SparseBitVect) -> typing.Sequence[double]:
    """
    """
@typing.overload
def OffBitProjSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> typing.Sequence[double]:
    """
    """
@typing.overload
def OffBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> typing.Sequence[int]:
    """
    """
@typing.overload
def OffBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> typing.Sequence[int]:
    """
        Returns the number of off bits in common between the two bit vectors
    
    """
@typing.overload
def OnBitProjSimilarity(bv1: SparseBitVect, bv2: SparseBitVect) -> typing.Sequence[double]:
    """
    """
@typing.overload
def OnBitProjSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> typing.Sequence[double]:
    """
        Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))
    
    """
@typing.overload
def OnBitSimilarity(v1: SparseBitVect, v2: SparseBitVect) -> float:
    """
    """
@typing.overload
def OnBitSimilarity(v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    """
        B(bv1&bv2) / B(bv1|bv2)
    
    """
@typing.overload
def OnBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> typing.Sequence[int]:
    """
    """
@typing.overload
def OnBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> typing.Sequence[int]:
    """
        Returns the number of on bits in common between the two bit vectors
    
    """
@typing.overload
def RogotGoldbergSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def RogotGoldbergSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def RogotGoldbergSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def RogotGoldbergSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / B(bv1)
    
    """
def RogotGoldbergSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
def RogotGoldbergSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def RusselSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def RusselSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def RusselSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def RusselSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / B(bv1)
    
    """
def RusselSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
def RusselSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / B(bv1)
    
    """
@typing.overload
def SokalSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def SokalSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))
    
    """
@typing.overload
def SokalSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def SokalSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))
    
    """
def SokalSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))
    
    """
def SokalSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))
    
    """
@typing.overload
def TanimotoSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def TanimotoSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))
    
    """
@typing.overload
def TanimotoSimilarity(bv1: SparseBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def TanimotoSimilarity(bv1: ExplicitBitVect, pkl: str, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))
    
    """
@typing.overload
def TanimotoSimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tanimoto similarity between two vectors
    
    """
@typing.overload
def TanimotoSimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tanimoto similarity between two vectors
    
    """
@typing.overload
def TanimotoSimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tanimoto similarity between two vectors
    
    """
@typing.overload
def TanimotoSimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tanimoto similarity between two vectors
    
    """
def TanimotoSimilarityNeighbors(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))
    
    """
def TanimotoSimilarityNeighbors_sparse(bvqueries: typing.Any, bvList: typing.Any) -> list:
    """
        B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))
    
    """
@typing.overload
def TverskySimilarity(bv1: SparseBitVect, bv2: SparseBitVect, a: float, b: float, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def TverskySimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, a: float, b: float, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)
    
    """
@typing.overload
def TverskySimilarity(bv1: SparseBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0) -> float:
    """
    """
@typing.overload
def TverskySimilarity(bv1: ExplicitBitVect, pkl: str, a: float, b: float, returnDistance: bool = 0) -> float:
    """
        B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2)
    
    """
@typing.overload
def TverskySimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tversky similarity between two vectors
    
    """
@typing.overload
def TverskySimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tversky similarity between two vectors
    
    """
@typing.overload
def TverskySimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tversky similarity between two vectors
    
    """
@typing.overload
def TverskySimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """
        return the Tversky similarity between two vectors
    
    """
EIGHTBITVALUE: DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE
FOURBITVALUE: DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE
ONEBITVALUE: DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE
SIXTEENBITVALUE: DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE
TWOBITVALUE: DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE
