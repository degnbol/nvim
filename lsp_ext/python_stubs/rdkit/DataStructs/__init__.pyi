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
import math as math
from rdkit.DataStructs.cDataStructs import DiscreteValueType
from rdkit.DataStructs.cDataStructs import DiscreteValueVect
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.DataStructs.cDataStructs import FPBReader
from rdkit.DataStructs.cDataStructs import IntSparseIntVect
from rdkit.DataStructs.cDataStructs import LongSparseIntVect
from rdkit.DataStructs.cDataStructs import MultiFPBReader
from rdkit.DataStructs.cDataStructs import RealValueVect
from rdkit.DataStructs.cDataStructs import SparseBitVect
from rdkit.DataStructs.cDataStructs import UIntSparseIntVect
from rdkit.DataStructs.cDataStructs import ULongSparseIntVect
from rdkit import rdBase
from .cDataStructs import *
__all__: list[str] = ['DiscreteValueType', 'DiscreteValueVect', 'EIGHTBITVALUE', 'ExplicitBitVect', 'FOURBITVALUE', 'FPBReader', 'FingerprintSimilarity', 'FoldToTargetDensity', 'IntSparseIntVect', 'LongSparseIntVect', 'MultiFPBReader', 'ONEBITVALUE', 'RealValueVect', 'SIXTEENBITVALUE', 'SparseBitVect', 'TWOBITVALUE', 'UIntSparseIntVect', 'ULongSparseIntVect', 'cDataStructs', 'getElementFromFlatMatrix', 'getNForFlatMatrix', 'math', 'rdBase', 'similarityFunctions']
def FingerprintSimilarity(fp1, fp2, metric = ...):
    """
     returns the calculated similarity between two fingerprints,
          handles any folding that may need to be done to ensure that they
          are compatible
    
        
    """
def FoldToTargetDensity(fp, density = 0.3, minLength = 64):
    ...
def getElementFromFlatMatrix(matrix, i, j):
    """
    Return element (i,j); diagonal is 0; lower side mirrors upper.
    """
def getNForFlatMatrix(matrix):
    """
    Get n for a strict upper- (or lower-) triangular matrix.
    """
EIGHTBITVALUE: cDataStructs.DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.EIGHTBITVALUE
FOURBITVALUE: cDataStructs.DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.FOURBITVALUE
ONEBITVALUE: cDataStructs.DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.ONEBITVALUE
SIXTEENBITVALUE: cDataStructs.DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.SIXTEENBITVALUE
TWOBITVALUE: cDataStructs.DiscreteValueType  # value = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE
similarityFunctions: list  # value = [('Tanimoto', <Boost.Python.function object>, ''), ('Dice', <Boost.Python.function object>, ''), ('Cosine', <Boost.Python.function object>, ''), ('Sokal', <Boost.Python.function object>, ''), ('Russel', <Boost.Python.function object>, ''), ('RogotGoldberg', <Boost.Python.function object>, ''), ('AllBit', <Boost.Python.function object>, ''), ('Kulczynski', <Boost.Python.function object>, ''), ('McConnaughey', <Boost.Python.function object>, ''), ('Asymmetric', <Boost.Python.function object>, ''), ('BraunBlanquet', <Boost.Python.function object>, '')]

# present at runtime, absent from the generated stub:
from rdkit.DataStructs.cDataStructs import AllBitSimilarity as AllBitSimilarity
from rdkit.DataStructs.cDataStructs import AllProbeBitsMatch as AllProbeBitsMatch
from rdkit.DataStructs.cDataStructs import AsymmetricSimilarity as AsymmetricSimilarity
from rdkit.DataStructs.cDataStructs import AsymmetricSimilarityNeighbors as AsymmetricSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import AsymmetricSimilarityNeighbors_sparse as AsymmetricSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import BitVectToBinaryText as BitVectToBinaryText
from rdkit.DataStructs.cDataStructs import BitVectToFPSText as BitVectToFPSText
from rdkit.DataStructs.cDataStructs import BitVectToText as BitVectToText
from rdkit.DataStructs.cDataStructs import BraunBlanquetSimilarity as BraunBlanquetSimilarity
from rdkit.DataStructs.cDataStructs import BraunBlanquetSimilarityNeighbors as BraunBlanquetSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import BraunBlanquetSimilarityNeighbors_sparse as BraunBlanquetSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import BulkAllBitSimilarity as BulkAllBitSimilarity
from rdkit.DataStructs.cDataStructs import BulkAsymmetricSimilarity as BulkAsymmetricSimilarity
from rdkit.DataStructs.cDataStructs import BulkBraunBlanquetSimilarity as BulkBraunBlanquetSimilarity
from rdkit.DataStructs.cDataStructs import BulkCosineSimilarity as BulkCosineSimilarity
from rdkit.DataStructs.cDataStructs import BulkDiceSimilarity as BulkDiceSimilarity
from rdkit.DataStructs.cDataStructs import BulkKulczynskiSimilarity as BulkKulczynskiSimilarity
from rdkit.DataStructs.cDataStructs import BulkMcConnaugheySimilarity as BulkMcConnaugheySimilarity
from rdkit.DataStructs.cDataStructs import BulkOnBitSimilarity as BulkOnBitSimilarity
from rdkit.DataStructs.cDataStructs import BulkRogotGoldbergSimilarity as BulkRogotGoldbergSimilarity
from rdkit.DataStructs.cDataStructs import BulkRusselSimilarity as BulkRusselSimilarity
from rdkit.DataStructs.cDataStructs import BulkSokalSimilarity as BulkSokalSimilarity
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity as BulkTanimotoSimilarity
from rdkit.DataStructs.cDataStructs import BulkTverskySimilarity as BulkTverskySimilarity
from rdkit.DataStructs.cDataStructs import ComputeL1Norm as ComputeL1Norm
from rdkit.DataStructs.cDataStructs import ConvertToExplicit as ConvertToExplicit
from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray as ConvertToNumpyArray
from rdkit.DataStructs.cDataStructs import CosineSimilarity as CosineSimilarity
from rdkit.DataStructs.cDataStructs import CosineSimilarityNeighbors as CosineSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import CosineSimilarityNeighbors_sparse as CosineSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import CreateFromBinaryText as CreateFromBinaryText
from rdkit.DataStructs.cDataStructs import CreateFromBitString as CreateFromBitString
from rdkit.DataStructs.cDataStructs import CreateFromFPSText as CreateFromFPSText
from rdkit.DataStructs.cDataStructs import DiceSimilarity as DiceSimilarity
from rdkit.DataStructs.cDataStructs import DiceSimilarityNeighbors as DiceSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import DiceSimilarityNeighbors_sparse as DiceSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import FoldFingerprint as FoldFingerprint
from rdkit.DataStructs.cDataStructs import InitFromDaylightString as InitFromDaylightString
from rdkit.DataStructs.cDataStructs import KulczynskiSimilarity as KulczynskiSimilarity
from rdkit.DataStructs.cDataStructs import KulczynskiSimilarityNeighbors as KulczynskiSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import KulczynskiSimilarityNeighbors_sparse as KulczynskiSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import McConnaugheySimilarity as McConnaugheySimilarity
from rdkit.DataStructs.cDataStructs import McConnaugheySimilarityNeighbors as McConnaugheySimilarityNeighbors
from rdkit.DataStructs.cDataStructs import McConnaugheySimilarityNeighbors_sparse as McConnaugheySimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import NumBitsInCommon as NumBitsInCommon
from rdkit.DataStructs.cDataStructs import OffBitProjSimilarity as OffBitProjSimilarity
from rdkit.DataStructs.cDataStructs import OffBitsInCommon as OffBitsInCommon
from rdkit.DataStructs.cDataStructs import OnBitProjSimilarity as OnBitProjSimilarity
from rdkit.DataStructs.cDataStructs import OnBitSimilarity as OnBitSimilarity
from rdkit.DataStructs.cDataStructs import OnBitsInCommon as OnBitsInCommon
from rdkit.DataStructs.cDataStructs import RogotGoldbergSimilarity as RogotGoldbergSimilarity
from rdkit.DataStructs.cDataStructs import RogotGoldbergSimilarityNeighbors as RogotGoldbergSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import RogotGoldbergSimilarityNeighbors_sparse as RogotGoldbergSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import RusselSimilarity as RusselSimilarity
from rdkit.DataStructs.cDataStructs import RusselSimilarityNeighbors as RusselSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import RusselSimilarityNeighbors_sparse as RusselSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import SokalSimilarity as SokalSimilarity
from rdkit.DataStructs.cDataStructs import SokalSimilarityNeighbors as SokalSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import SokalSimilarityNeighbors_sparse as SokalSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity as TanimotoSimilarity
from rdkit.DataStructs.cDataStructs import TanimotoSimilarityNeighbors as TanimotoSimilarityNeighbors
from rdkit.DataStructs.cDataStructs import TanimotoSimilarityNeighbors_sparse as TanimotoSimilarityNeighbors_sparse
from rdkit.DataStructs.cDataStructs import TverskySimilarity as TverskySimilarity
from . import cDataStructs as cDataStructs
