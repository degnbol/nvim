# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Import all RDKit chemistry modules

"""
from __future__ import annotations
from collections import namedtuple
import numpy as numpy
from rdkit.Chem import CanonSmiles
from rdkit.Chem.ChemicalFeatures import MCFF_GetFeaturesForMol
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem import MultiConfMolFromSDF
from rdkit.Chem import QuickSmartsMatch
from rdkit.Chem import SupplierFromFilename
from rdkit.Chem import inchi
from rdkit.Chem.inchi import InchiReadWriteError
from rdkit.Chem.inchi import InchiToInchiKey
from rdkit.Chem.inchi import MolBlockToInchi
from rdkit.Chem.inchi import MolBlockToInchiAndAuxInfo
from rdkit.Chem.inchi import MolFromInchi
from rdkit.Chem.inchi import MolFromInchiAndAuxInfo
from rdkit.Chem.inchi import MolToInchi
from rdkit.Chem.inchi import MolToInchiAndAuxInfo
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem import rdCIPLabeler
import rdkit.Chem.rdChemReactions
from rdkit.Chem.rdChemReactions import CartesianProductStrategy
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem.rdChemReactions import EnumerateLibrary
from rdkit.Chem.rdChemReactions import EnumerateLibraryBase
from rdkit.Chem.rdChemReactions import EnumerationParams
from rdkit.Chem.rdChemReactions import EnumerationStrategyBase
from rdkit.Chem.rdChemReactions import EvenSamplePairsStrategy
from rdkit.Chem.rdChemReactions import FingerprintType
from rdkit.Chem.rdChemReactions import MOL_SPTR_VECT
from rdkit.Chem.rdChemReactions import RandomSampleAllBBsStrategy
from rdkit.Chem.rdChemReactions import RandomSampleStrategy
from rdkit.Chem.rdChemReactions import ReactionFingerprintParams
from rdkit.Chem.rdChemReactions import SanitizeFlags
from rdkit.Chem.rdChemReactions import VectMolVect
from rdkit.Chem.rdChemicalFeatures import FreeChemicalFeature
from rdkit.Chem import rdCoordGen
from rdkit.Chem.rdDepictor import ConstrainedDepictionParams
from rdkit.Chem.rdDepictor import UsingCoordGen
import rdkit.Chem.rdDistGeom
from rdkit.Chem.rdDistGeom import EmbedFailureCauses
from rdkit.Chem.rdDistGeom import EmbedParameters
import rdkit.Chem.rdFingerprintGenerator
from rdkit.Chem.rdFingerprintGenerator import AdditionalOutput
from rdkit.Chem.rdFingerprintGenerator import AtomInvariantsGenerator
from rdkit.Chem.rdFingerprintGenerator import AtomPairFingerprintOptions
from rdkit.Chem.rdFingerprintGenerator import BondInvariantsGenerator
from rdkit.Chem.rdFingerprintGenerator import FPType
from rdkit.Chem.rdFingerprintGenerator import FingerprintGenerator32
from rdkit.Chem.rdFingerprintGenerator import FingerprintGenerator64
from rdkit.Chem.rdFingerprintGenerator import FingerprintOptions
from rdkit.Chem.rdFingerprintGenerator import MorganFingerprintOptions
from rdkit.Chem.rdFingerprintGenerator import RDKitFingerprintOptions
from rdkit.Chem.rdFingerprintGenerator import TopologicalTorsionFingerprintOptions
from rdkit.Chem.rdMolAlign import BestAlignmentParams
from rdkit.Chem.rdMolAlign import O3A
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeatureFactory
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Chem.rdMolDescriptors import DoubleCubicLatticeVolume
from rdkit.Chem.rdMolDescriptors import NumRotatableBondsOptions
from rdkit.Chem.rdMolDescriptors import Properties
from rdkit.Chem.rdMolDescriptors import PropertyFunctor
from rdkit.Chem.rdMolDescriptors import PropertyRangeQuery
from rdkit.Chem.rdMolDescriptors import PythonPropertyFunctor
from rdkit.Chem.rdMolEnumerator import EnumeratorType
from rdkit.Chem.rdMolEnumerator import MolEnumeratorParams
from rdkit.Chem import rdMolInterchange
from rdkit.Chem.rdMolInterchange import JSONParseParameters
from rdkit.Chem.rdMolInterchange import JSONWriteParameters
import rdkit.Chem.rdchem
from rdkit.Chem import rdchem
from rdkit.Chem.rdchem import Atom
from rdkit.Chem.rdchem import AtomCoordsMatcher
from rdkit.Chem.rdchem import AtomKekulizeException
from rdkit.Chem.rdchem import AtomMonomerInfo
from rdkit.Chem.rdchem import AtomMonomerType
from rdkit.Chem.rdchem import AtomPDBResidueInfo
from rdkit.Chem.rdchem import AtomSanitizeException
from rdkit.Chem.rdchem import AtomValenceException
from rdkit.Chem.rdchem import Bond
from rdkit.Chem.rdchem import BondDir
from rdkit.Chem.rdchem import BondStereo
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdchem import ChiralType
from rdkit.Chem.rdchem import CompositeQueryType
from rdkit.Chem.rdchem import Conformer
from rdkit.Chem.rdchem import EditableMol
from rdkit.Chem.rdchem import FixedMolSizeMolBundle
from rdkit.Chem.rdchem import HybridizationType
from rdkit.Chem.rdchem import KekulizeException
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdchem import MolBundle
from rdkit.Chem.rdchem import MolSanitizeException
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem.rdchem import PropertyPickleOptions
from rdkit.Chem.rdchem import QueryAtom
from rdkit.Chem.rdchem import QueryBond
from rdkit.Chem.rdchem import RWMol
from rdkit.Chem.rdchem import ResonanceFlags
from rdkit.Chem.rdchem import ResonanceMolSupplier
from rdkit.Chem.rdchem import ResonanceMolSupplierCallback
from rdkit.Chem.rdchem import RingInfo
from rdkit.Chem.rdchem import StereoDescriptor
from rdkit.Chem.rdchem import StereoGroup
from rdkit.Chem.rdchem import StereoGroupType
from rdkit.Chem.rdchem import StereoGroup_vect
from rdkit.Chem.rdchem import StereoInfo
from rdkit.Chem.rdchem import StereoSpecified
from rdkit.Chem.rdchem import StereoType
from rdkit.Chem.rdchem import SubstanceGroup
from rdkit.Chem.rdchem import SubstanceGroupAttach
from rdkit.Chem.rdchem import SubstanceGroupCState
from rdkit.Chem.rdchem import SubstanceGroup_VECT
from rdkit.Chem.rdchem import SubstructMatchParameters
from rdkit.Chem.rdchem import ValenceType
from rdkit.Chem import rdinchi
from rdkit.Chem import rdmolfiles
from rdkit.Chem.rdmolfiles import CDXMLFormat
from rdkit.Chem.rdmolfiles import CDXMLParserParams
from rdkit.Chem.rdmolfiles import CXSmilesFields
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
from rdkit.Chem.rdmolfiles import MaeMolSupplier
from rdkit.Chem.rdmolfiles import MaeWriter
from rdkit.Chem.rdmolfiles import MolFromSCSRParams
from rdkit.Chem.rdmolfiles import MolWriterParams
from rdkit.Chem.rdmolfiles import MultithreadedSDMolSupplier
from rdkit.Chem.rdmolfiles import MultithreadedSmilesMolSupplier
from rdkit.Chem.rdmolfiles import PDBWriter
from rdkit.Chem.rdmolfiles import PNGMetadataParams
from rdkit.Chem.rdmolfiles import RestoreBondDirOption
from rdkit.Chem.rdmolfiles import SCSRBaseHbondOptions
from rdkit.Chem.rdmolfiles import SCSRTemplateNames
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem.rdmolfiles import SmartsParserParams
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from rdkit.Chem.rdmolfiles import SmilesParserParams
from rdkit.Chem.rdmolfiles import SmilesWriteParams
from rdkit.Chem.rdmolfiles import SmilesWriter
from rdkit.Chem.rdmolfiles import TDTMolSupplier
from rdkit.Chem.rdmolfiles import TDTWriter
import rdkit.Chem.rdmolops
from rdkit.Chem import rdmolops
from rdkit.Chem.rdmolops import AddHsParameters
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.Chem.rdmolops import AdjustQueryWhichFlags
from rdkit.Chem.rdmolops import AromaticityModel
from rdkit.Chem.rdmolops import BondWedgingParameters
from rdkit.Chem.rdmolops import BoolVector
from rdkit.Chem.rdmolops import MolzipLabel
from rdkit.Chem.rdmolops import MolzipParams
from rdkit.Chem.rdmolops import RemoveHsParameters
from rdkit.Chem.rdmolops import StereoBondThresholds
from rdkit.Chem.rdmolops import StereoGroupAbsOptions
from rdkit.Chem.rdmolops import SubsetInfo
from rdkit.Chem.rdmolops import SubsetMethod
from rdkit.Chem.rdmolops import SubsetOptions
from rdkit.Chem.rdmolops import UIntUIntMap
from rdkit.Chem.rdmolops import map_indexing_suite_UIntUIntMap_entry
from rdkit import DataStructs
from rdkit import ForceField
from rdkit.Geometry import rdGeometry
from rdkit import RDConfig
import rdkit.RDLogger
from rdkit import rdBase
import sys as sys
import warnings as warnings
__all__: list[str] = ['ADJUST_IGNOREALL', 'ADJUST_IGNORECHAINS', 'ADJUST_IGNOREDUMMIES', 'ADJUST_IGNOREMAPPED', 'ADJUST_IGNORENONDUMMIES', 'ADJUST_IGNORENONE', 'ADJUST_IGNORERINGS', 'ALLOW_CHARGE_SEPARATION', 'ALLOW_INCOMPLETE_OCTETS', 'AROMATICITY_CUSTOM', 'AROMATICITY_DEFAULT', 'AROMATICITY_MDL', 'AROMATICITY_MMFF94', 'AROMATICITY_RDKIT', 'AROMATICITY_SIMPLE', 'AddHsParameters', 'AdditionalOutput', 'AdjustQueryParameters', 'AdjustQueryWhichFlags', 'AllProps', 'AromaticityModel', 'AssignBondOrdersFromTemplate', 'Atom', 'AtomCoordsMatcher', 'AtomInvariantsGenerator', 'AtomKekulizeException', 'AtomMonomerInfo', 'AtomMonomerType', 'AtomPDBResidueInfo', 'AtomPairFP', 'AtomPairFingerprintOptions', 'AtomPairsParameters', 'AtomProps', 'AtomSanitizeException', 'AtomValenceException', 'BAD_DOUBLE_BOND_STEREO', 'BestAlignmentParams', 'Bond', 'BondDir', 'BondInvariantsGenerator', 'BondProps', 'BondStereo', 'BondType', 'BondWedgingParameters', 'BoolVector', 'CDXMLFormat', 'CDXMLParserParams', 'CHECK_CHIRAL_CENTERS', 'CHECK_CHIRAL_CENTERS2', 'CHECK_TETRAHEDRAL_CENTERS', 'CHI_ALLENE', 'CHI_OCTAHEDRAL', 'CHI_OTHER', 'CHI_SQUAREPLANAR', 'CHI_TETRAHEDRAL', 'CHI_TETRAHEDRAL_CCW', 'CHI_TETRAHEDRAL_CW', 'CHI_TRIGONALBIPYRAMIDAL', 'CHI_UNSPECIFIED', 'CLASH', 'COMPOSITE_AND', 'COMPOSITE_OR', 'COMPOSITE_XOR', 'CXSmilesFields', 'CanonSmiles', 'CartesianProductStrategy', 'ChemicalReaction', 'ChiralType', 'CompositeQueryType', 'ComputeMolShape', 'ComputeMolVolume', 'ComputedProps', 'Conformer', 'ConstrainedDepictionParams', 'ConstrainedEmbed', 'CoordsAsDouble', 'DataStructs', 'DoubleCubicLatticeVolume', 'ETK_MINIMIZATION', 'EXCEEDED_TIMEOUT', 'EXPLICIT', 'EditableMol', 'EmbedFailureCauses', 'EmbedParameters', 'EnumerateLibrary', 'EnumerateLibraryBase', 'EnumerateLibraryFromReaction', 'EnumerateStereoisomers', 'EnumerationParams', 'EnumerationStrategyBase', 'EnumeratorType', 'EvenSamplePairsStrategy', 'FINAL_CENTER_IN_VOLUME', 'FINAL_CHIRAL_BOUNDS', 'FIRST_MINIMIZATION', 'FPType', 'FindMolChiralCenters', 'FingerprintGenerator32', 'FingerprintGenerator64', 'FingerprintOptions', 'FingerprintType', 'FixedMolSizeMolBundle', 'ForceField', 'ForwardSDMolSupplier', 'FreeChemicalFeature', 'GetConformerRMS', 'GetConformerRMSMatrix', 'HybridizationType', 'IMPLICIT', 'INCHI_AVAILABLE', 'INITIAL_COORDS', 'InchiReadWriteError', 'InchiToInchiKey', 'JSONParseParameters', 'JSONWriteParameters', 'KEKULE_ALL', 'KTERM_VIOLATION', 'KekulizeException', 'LINEAR_DOUBLE_BOND', 'LayeredFingerprint_substructLayers', 'MCFF_GetFeaturesForMol', 'MINIMIZATION', 'MINIMIZE_FOURTH_DIMENSION', 'MOL_SPTR_VECT', 'MaeMolSupplier', 'MaeWriter', 'Mol', 'MolBlockToInchi', 'MolBlockToInchiAndAuxInfo', 'MolBundle', 'MolChemicalFeature', 'MolChemicalFeatureFactory', 'MolEnumeratorParams', 'MolFromInchi', 'MolFromInchiAndAuxInfo', 'MolFromSCSRParams', 'MolProps', 'MolSanitizeException', 'MolToInchi', 'MolToInchiAndAuxInfo', 'MolToInchiKey', 'MolWriterParams', 'MolzipLabel', 'MolzipParams', 'MorganFP', 'MorganFingerprintOptions', 'MultiConfMolFromSDF', 'MultithreadedSDMolSupplier', 'MultithreadedSmilesMolSupplier', 'NoConformers', 'NoProps', 'NumRotatableBondsOptions', 'O3A', 'PDBWriter', 'PNGMetadataParams', 'PeriodicTable', 'PrivateProps', 'Properties', 'PropertyFunctor', 'PropertyPickleOptions', 'PropertyRangeQuery', 'PythonPropertyFunctor', 'QueryAtom', 'QueryAtomData', 'QueryBond', 'QuickSmartsMatch', 'RDConfig', 'RDKitFP', 'RDKitFingerprintOptions', 'RWMol', 'RandomSampleAllBBsStrategy', 'RandomSampleStrategy', 'ReactionFingerprintParams', 'RemoveHsParameters', 'ResonanceFlags', 'ResonanceMolSupplier', 'ResonanceMolSupplierCallback', 'RestoreBondDirOption', 'RingInfo', 'SANITIZE_ADJUSTHS', 'SANITIZE_ADJUST_REACTANTS', 'SANITIZE_ALL', 'SANITIZE_ATOM_MAPS', 'SANITIZE_CLEANUP', 'SANITIZE_CLEANUPATROPISOMERS', 'SANITIZE_CLEANUPCHIRALITY', 'SANITIZE_CLEANUP_ORGANOMETALLICS', 'SANITIZE_FINDRADICALS', 'SANITIZE_KEKULIZE', 'SANITIZE_MERGEHS', 'SANITIZE_NONE', 'SANITIZE_PROPERTIES', 'SANITIZE_RGROUP_NAMES', 'SANITIZE_SETAROMATICITY', 'SANITIZE_SETCONJUGATION', 'SANITIZE_SETHYBRIDIZATION', 'SANITIZE_SYMMRINGS', 'SCSRBaseHbondOptions', 'SCSRTemplateNames', 'SDMolSupplier', 'SDWriter', 'STEREO_ABSOLUTE', 'STEREO_AND', 'STEREO_OR', 'SanitizeFlags', 'SmartsParserParams', 'SmilesMolSupplier', 'SmilesParserParams', 'SmilesWriteParams', 'SmilesWriter', 'StereoBondThresholds', 'StereoDescriptor', 'StereoEnumerationOptions', 'StereoGroup', 'StereoGroupAbsOptions', 'StereoGroupType', 'StereoGroup_vect', 'StereoInfo', 'StereoSpecified', 'StereoType', 'SubsetInfo', 'SubsetMethod', 'SubsetOptions', 'SubstanceGroup', 'SubstanceGroupAttach', 'SubstanceGroupCState', 'SubstanceGroup_VECT', 'SubstructMatchParameters', 'SupplierFromFilename', 'TDTMolSupplier', 'TDTWriter', 'TopologicalTorsionFP', 'TopologicalTorsionFingerprintOptions', 'TransformMol', 'UIntUIntMap', 'UNCONSTRAINED_ANIONS', 'UNCONSTRAINED_CATIONS', 'UsingCoordGen', 'ValenceType', 'VectMolVect', 'inchi', 'logger', 'map_indexing_suite_UIntUIntMap_entry', 'namedtuple', 'numpy', 'rdBase', 'rdCIPLabeler', 'rdCoordGen', 'rdGeometry', 'rdMolInterchange', 'rdchem', 'rdinchi', 'rdmolfiles', 'rdmolops', 'sys', 'templDir', 'warnings']
def AssignBondOrdersFromTemplate(refmol, mol):
    """
     assigns bond orders to a molecule based on the
        bond orders in a template molecule
    
        Arguments
          - refmol: the template molecule
          - mol: the molecule to assign bond orders to
    
        An example, start by generating a template from a SMILES
        and read in the PDB structure of the molecule
    
        >>> import os
        >>> from rdkit.Chem import AllChem
        >>> template = AllChem.MolFromSmiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")
        >>> mol = AllChem.MolFromPDBFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4DJU_lig.pdb'))
        >>> len([1 for b in template.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
        8
        >>> len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
        22
    
        Now assign the bond orders based on the template molecule
    
        >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        >>> len([1 for b in newMol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
        8
    
        Note that the template molecule should have no explicit hydrogens
        else the algorithm will fail.
    
        It also works if there are different formal charges (this was github issue 235):
    
        >>> template=AllChem.MolFromSmiles('CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]')
        >>> mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4FTR_lig.mol'))
        >>> AllChem.MolToSmiles(mol)
        'COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O'
        >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
        >>> AllChem.MolToSmiles(newMol)
        'COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]'
    
        
    """
def ComputeMolShape(mol, confId = -1, boxDim = (20, 20, 20), spacing = 0.5, **kwargs):
    """
     returns a grid representation of the molecule's shape
    
      Arguments:
        - mol: the molecule
        - confId: (optional) the conformer id to use
        - boxDim: (optional) the dimensions of the box to use
        - spacing: (optional) the spacing to use
        - kwargs: additional arguments to pass to the encoding function
    
      Returns:
        a UniformGrid3D object
      
    """
def ComputeMolVolume(mol, confId = -1, gridSpacing = 0.2, boxMargin = 2.0):
    """
     Calculates the volume of a particular conformer of a molecule
        based on a grid-encoding of the molecular shape.
    
    
      Arguments:
        - mol: the molecule
        - confId: (optional) the conformer id to use
        - gridSpacing: (optional) the spacing to use 
        - boxMargin: (optional) the margin to use around the molecule
    
        A bit of demo as well as a test of github #1883:
    
        >>> from rdkit import Chem
        >>> from rdkit.Chem import AllChem
        >>> mol = Chem.AddHs(Chem.MolFromSmiles('C'))
        >>> AllChem.EmbedMolecule(mol)
        0
        >>> ComputeMolVolume(mol)
        28...
        >>> mol = Chem.AddHs(Chem.MolFromSmiles('O'))
        >>> AllChem.EmbedMolecule(mol)
        0
        >>> ComputeMolVolume(mol)
        20...
    
        
    """
def ConstrainedEmbed(mol, core, useTethers = True, coreConfId = -1, randomseed = 2342, getForceField = ..., **kwargs):
    """
     generates an embedding of a molecule where part of the molecule
        is constrained to have particular coordinates
    
        Arguments
          - mol: the molecule to embed
          - core: the molecule to use as a source of constraints
          - useTethers: (optional) if True, the final conformation will be
              optimized subject to a series of extra forces that pull the
              matching atoms to the positions of the core atoms. Otherwise
              simple distance constraints based on the core atoms will be
              used in the optimization.
          - coreConfId: (optional) id of the core conformation to use
          - randomSeed: (optional) seed for the random number generator
          - getForceField: (optional) a function to use to get a force field
              for the final cleanup
          - kwargs: additional arguments to pass to the embedding function
    
        An example, start by generating a template with a 3D structure:
    
        >>> from rdkit.Chem import AllChem
        >>> template = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1")
        >>> AllChem.EmbedMolecule(template)
        0
        >>> AllChem.UFFOptimizeMolecule(template)
        0
    
        Here's a molecule:
    
        >>> mol = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")
    
        Now do the constrained embedding
      
        >>> mol = AllChem.ConstrainedEmbed(mol, template)
    
        Demonstrate that the positions are nearly the same with template:
    
        >>> import math
        >>> molp = mol.GetConformer().GetAtomPosition(0)
        >>> templatep = template.GetConformer().GetAtomPosition(0)
        >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
        True
        >>> molp = mol.GetConformer().GetAtomPosition(1)
        >>> templatep = template.GetConformer().GetAtomPosition(1)
        >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
        True
    
        
    """
def EnumerateLibraryFromReaction(reaction, sidechainSets, returnReactants = False):
    """
     Returns a generator for the virtual library defined by
        a reaction and a sequence of sidechain sets
    
        Arguments:
          - reaction: the reaction
          - sidechainSets: a sequence of sequences of sidechains
          - returnReactants: (optional) if True, the generator will return information about the reactants
                             as well as the products
    
        >>> from rdkit import Chem
        >>> from rdkit.Chem import AllChem
        >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
        >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
        >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
        >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
        >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
        ['CNC=O', 'CCNC=O', 'CNC(C)=O', 'CCNC(C)=O']
    
        Note that this is all done in a lazy manner, so "infinitely" large libraries can
        be done without worrying about running out of memory. Your patience will run out first:
    
        Define a set of 10000 amines:
    
        >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))
    
        ... a set of 10000 acids
    
        >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))
    
        ... now the virtual library (1e8 compounds in principle):
    
        >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])
    
        ... look at the first 4 compounds:
    
        >>> [Chem.MolToSmiles(next(r)[0]) for x in range(4)]
        ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']
    
        Here's what returnReactants does:
    
        >>> l = list(AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1],returnReactants=True))
        >>> type(l[0])
        <class 'rdkit.Chem.AllChem.ProductReactants'>
        >>> [Chem.MolToSmiles(x) for x in l[0].reactants]
        ['O=CO', 'CN']
        >>> [Chem.MolToSmiles(x) for x in l[0].products]
        ['CNC=O']
    
        
    """
def GetConformerRMS(mol, confId1, confId2, atomIds = None, prealigned = False):
    """
     Returns the RMS between two conformations.
        By default, the conformers will be aligned to the first conformer
        before the RMS calculation and, as a side-effect, the second will be left
        in the aligned state.
    
        Arguments:
          - mol:        the molecule
          - confId1:    the id of the first conformer
          - confId2:    the id of the second conformer
          - atomIds:    (optional) list of atom ids to use a points for
                        alingment - defaults to all atoms
          - prealigned: (optional) by default the conformers are assumed
                        be unaligned and the second conformer be aligned
                        to the first
    
        
    """
def GetConformerRMSMatrix(mol, atomIds = None, prealigned = False):
    """
     Returns the RMS matrix of the conformers of a molecule.
        As a side-effect, the conformers will be aligned to the first
        conformer (i.e. the reference) and will left in the aligned state.
    
        Arguments:
          - mol:     the molecule
          - atomIds: (optional) list of atom ids to use a points for
                     alingment - defaults to all atoms
          - prealigned: (optional) by default the conformers are assumed
                        be unaligned and will therefore be aligned to the
                        first conformer
    
        Note that the returned RMS matrix is symmetrical, i.e. it is the
        lower half of the matrix, e.g. for 5 conformers::
    
          rmsmatrix = [ a,
                        b, c,
                        d, e, f,
                        g, h, i, j]
    
        where a is the RMS between conformers 0 and 1, b is the RMS between
        conformers 0 and 2, etc.
        This way it can be directly used as distance matrix in e.g. Butina
        clustering.
    
        
    """
def TransformMol(mol, tform, confId = -1, keepConfs = False):
    """
      Applies the transformation (usually a 4x4 double matrix) to a molecule
        
      Arguments:
        - mol: the molecule to be transformed
        - tform: the transformation to apply
        - confId: (optional) the conformer id to transform
        - keepConfs: (optional) if keepConfs is False then all but that conformer are removed
      
      
    """
ADJUST_IGNOREALL: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL
ADJUST_IGNORECHAINS: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
ADJUST_IGNOREDUMMIES: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
ADJUST_IGNOREMAPPED: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED
ADJUST_IGNORENONDUMMIES: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES
ADJUST_IGNORENONE: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE
ADJUST_IGNORERINGS: rdkit.Chem.rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS
ALLOW_CHARGE_SEPARATION: rdkit.Chem.rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
ALLOW_INCOMPLETE_OCTETS: rdkit.Chem.rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
AROMATICITY_CUSTOM: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM
AROMATICITY_DEFAULT: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT
AROMATICITY_MDL: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL
AROMATICITY_MMFF94: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94
AROMATICITY_RDKIT: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT
AROMATICITY_SIMPLE: rdkit.Chem.rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE
AllProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
AtomPairFP: rdkit.Chem.rdFingerprintGenerator.FPType  # value = rdkit.Chem.rdFingerprintGenerator.FPType.AtomPairFP
AtomProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
BAD_DOUBLE_BOND_STEREO: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
BondProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
CHECK_CHIRAL_CENTERS: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
CHECK_CHIRAL_CENTERS2: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS2
CHECK_TETRAHEDRAL_CENTERS: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
CHI_ALLENE: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
CHI_OCTAHEDRAL: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
CHI_OTHER: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
CHI_SQUAREPLANAR: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
CHI_TETRAHEDRAL: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
CHI_TETRAHEDRAL_CCW: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
CHI_TETRAHEDRAL_CW: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
CHI_TRIGONALBIPYRAMIDAL: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
CHI_UNSPECIFIED: rdkit.Chem.rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
CLASH: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CLASH
COMPOSITE_AND: rdkit.Chem.rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
COMPOSITE_OR: rdkit.Chem.rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
COMPOSITE_XOR: rdkit.Chem.rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
ComputedProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
CoordsAsDouble: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
ETK_MINIMIZATION: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
EXCEEDED_TIMEOUT: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.EXCEEDED_TIMEOUT
EXPLICIT: rdkit.Chem.rdchem.ValenceType  # value = rdkit.Chem.rdchem.ValenceType.EXPLICIT
FINAL_CENTER_IN_VOLUME: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
FINAL_CHIRAL_BOUNDS: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
FIRST_MINIMIZATION: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
IMPLICIT: rdkit.Chem.rdchem.ValenceType  # value = rdkit.Chem.rdchem.ValenceType.IMPLICIT
INCHI_AVAILABLE: bool = True
INITIAL_COORDS: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
KEKULE_ALL: rdkit.Chem.rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
KTERM_VIOLATION: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.KTERM_VIOLATION
LINEAR_DOUBLE_BOND: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
LayeredFingerprint_substructLayers: int = 7
MINIMIZATION: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZATION
MINIMIZE_FOURTH_DIMENSION: rdkit.Chem.rdDistGeom.EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
MolProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
MorganFP: rdkit.Chem.rdFingerprintGenerator.FPType  # value = rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP
NoConformers: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData: rdkit.Chem.rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
RDKitFP: rdkit.Chem.rdFingerprintGenerator.FPType  # value = rdkit.Chem.rdFingerprintGenerator.FPType.RDKitFP
SANITIZE_ADJUSTHS: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
SANITIZE_ADJUST_REACTANTS: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ADJUST_REACTANTS
SANITIZE_ALL: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ALL
SANITIZE_ATOM_MAPS: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_ATOM_MAPS
SANITIZE_CLEANUP: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
SANITIZE_CLEANUPATROPISOMERS: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS
SANITIZE_CLEANUPCHIRALITY: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
SANITIZE_CLEANUP_ORGANOMETALLICS: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
SANITIZE_FINDRADICALS: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
SANITIZE_KEKULIZE: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
SANITIZE_MERGEHS: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_MERGEHS
SANITIZE_NONE: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_NONE
SANITIZE_PROPERTIES: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
SANITIZE_RGROUP_NAMES: rdkit.Chem.rdChemReactions.SanitizeFlags  # value = rdkit.Chem.rdChemReactions.SanitizeFlags.SANITIZE_RGROUP_NAMES
SANITIZE_SETAROMATICITY: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
SANITIZE_SETCONJUGATION: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
SANITIZE_SETHYBRIDIZATION: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
SANITIZE_SYMMRINGS: rdkit.Chem.rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
STEREO_ABSOLUTE: rdkit.Chem.rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND: rdkit.Chem.rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR: rdkit.Chem.rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
TopologicalTorsionFP: rdkit.Chem.rdFingerprintGenerator.FPType  # value = rdkit.Chem.rdFingerprintGenerator.FPType.TopologicalTorsionFP
UNCONSTRAINED_ANIONS: rdkit.Chem.rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS: rdkit.Chem.rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
logger: rdkit.RDLogger.logger  # value = <rdkit.RDLogger.logger object>
templDir: str = '/Users/runner/work/rdkit-pypi/rdkit-pypi/build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/share/RDKit/Data/'

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdqueries import AAtomQueryAtom as AAtomQueryAtom
from rdkit.Chem.rdqueries import AHAtomQueryAtom as AHAtomQueryAtom
from rdkit.Chem.rdmolops import AddHs as AddHs
from rdkit.Chem.rdmolfiles import AddMetadataToPNGFile as AddMetadataToPNGFile
from rdkit.Chem.rdmolfiles import AddMetadataToPNGString as AddMetadataToPNGString
from rdkit.Chem.rdchem import AddMolSubstanceGroup as AddMolSubstanceGroup
from rdkit.Chem.rdmolops import AddRecursiveQuery as AddRecursiveQuery
from rdkit.Chem.rdDepictor import AddRingSystemTemplates as AddRingSystemTemplates
from rdkit.Chem.rdmolops import AddStereoAnnotations as AddStereoAnnotations
from rdkit.Chem.rdmolops import AddWavyBondsForStereoAny as AddWavyBondsForStereoAny
from rdkit.Chem.rdmolops import AdjustQueryProperties as AdjustQueryProperties
from rdkit.Chem.rdmolops import AdjustQueryPropertiesWithGenericGroups as AdjustQueryPropertiesWithGenericGroups
from rdkit.Chem.rdMolAlign import AlignMol as AlignMol
from rdkit.Chem.rdMolAlign import AlignMolConformers as AlignMolConformers
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromMolParity as AssignAtomChiralTagsFromMolParity
from rdkit.Chem.rdmolops import AssignAtomChiralTagsFromStructure as AssignAtomChiralTagsFromStructure
from rdkit.Chem.rdCIPLabeler import AssignCIPLabels as AssignCIPLabels
from rdkit.Chem.rdmolops import AssignChiralTypesFromBondDirs as AssignChiralTypesFromBondDirs
from rdkit.Chem.rdmolops import AssignRadicals as AssignRadicals
from rdkit.Chem.rdmolops import AssignStereochemistry as AssignStereochemistry
from rdkit.Chem.rdmolops import AssignStereochemistryFrom3D as AssignStereochemistryFrom3D
from rdkit.Chem.rdmolfiles import AtomFromSmarts as AtomFromSmarts
from rdkit.Chem.rdmolfiles import AtomFromSmiles as AtomFromSmiles
from rdkit.Chem.rdmolops import AtomHasConjugatedBond as AtomHasConjugatedBond
from rdkit.Chem.rdqueries import AtomNumEqualsQueryAtom as AtomNumEqualsQueryAtom
from rdkit.Chem.rdqueries import AtomNumGreaterQueryAtom as AtomNumGreaterQueryAtom
from rdkit.Chem.rdqueries import AtomNumLessQueryAtom as AtomNumLessQueryAtom
from rdkit.Chem.rdMolDescriptors import BCUT2D as BCUT2D
from rdkit.Chem.rdmolfiles import BondFromSmarts as BondFromSmarts
from rdkit.Chem.rdmolfiles import BondFromSmiles as BondFromSmiles
from rdkit.Chem.rdMolChemicalFeatures import BuildFeatureFactory as BuildFeatureFactory
from rdkit.Chem.rdMolChemicalFeatures import BuildFeatureFactoryFromString as BuildFeatureFactoryFromString
from rdkit.Chem.rdMolDescriptors import CalcAUTOCORR2D as CalcAUTOCORR2D
from rdkit.Chem.rdMolDescriptors import CalcAUTOCORR3D as CalcAUTOCORR3D
from rdkit.Chem.rdMolDescriptors import CalcAsphericity as CalcAsphericity
from rdkit.Chem.rdMolDescriptors import CalcChi0n as CalcChi0n
from rdkit.Chem.rdMolDescriptors import CalcChi0v as CalcChi0v
from rdkit.Chem.rdMolDescriptors import CalcChi1n as CalcChi1n
from rdkit.Chem.rdMolDescriptors import CalcChi1v as CalcChi1v
from rdkit.Chem.rdMolDescriptors import CalcChi2n as CalcChi2n
from rdkit.Chem.rdMolDescriptors import CalcChi2v as CalcChi2v
from rdkit.Chem.rdMolDescriptors import CalcChi3n as CalcChi3n
from rdkit.Chem.rdMolDescriptors import CalcChi3v as CalcChi3v
from rdkit.Chem.rdMolDescriptors import CalcChi4n as CalcChi4n
from rdkit.Chem.rdMolDescriptors import CalcChi4v as CalcChi4v
from rdkit.Chem.rdMolDescriptors import CalcChiNn as CalcChiNn
from rdkit.Chem.rdMolDescriptors import CalcChiNv as CalcChiNv
from rdkit.Chem.rdMolDescriptors import CalcCoulombMat as CalcCoulombMat
from rdkit.Chem.rdMolDescriptors import CalcCrippenDescriptors as CalcCrippenDescriptors
from rdkit.Chem.rdMolDescriptors import CalcEEMcharges as CalcEEMcharges
from rdkit.Chem.rdMolDescriptors import CalcEccentricity as CalcEccentricity
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt as CalcExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcFractionCSP3 as CalcFractionCSP3
from rdkit.Chem.rdMolDescriptors import CalcGETAWAY as CalcGETAWAY
from rdkit.Chem.rdMolDescriptors import CalcHallKierAlpha as CalcHallKierAlpha
from rdkit.Chem.rdMolDescriptors import CalcInertialShapeFactor as CalcInertialShapeFactor
from rdkit.Chem.rdMolDescriptors import CalcKappa1 as CalcKappa1
from rdkit.Chem.rdMolDescriptors import CalcKappa2 as CalcKappa2
from rdkit.Chem.rdMolDescriptors import CalcKappa3 as CalcKappa3
from rdkit.Chem.rdMolDescriptors import CalcLabuteASA as CalcLabuteASA
from rdkit.Chem.rdMolDescriptors import CalcMORSE as CalcMORSE
from rdkit.Chem.rdMolDescriptors import CalcMolFormula as CalcMolFormula
from rdkit.Chem.rdMolDescriptors import CalcNPR1 as CalcNPR1
from rdkit.Chem.rdMolDescriptors import CalcNPR2 as CalcNPR2
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticCarbocycles as CalcNumAliphaticCarbocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticHeterocycles as CalcNumAliphaticHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAliphaticRings as CalcNumAliphaticRings
from rdkit.Chem.rdMolDescriptors import CalcNumAmideBonds as CalcNumAmideBonds
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticCarbocycles as CalcNumAromaticCarbocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticHeterocycles as CalcNumAromaticHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumAromaticRings as CalcNumAromaticRings
from rdkit.Chem.rdMolDescriptors import CalcNumAtomStereoCenters as CalcNumAtomStereoCenters
from rdkit.Chem.rdMolDescriptors import CalcNumAtoms as CalcNumAtoms
from rdkit.Chem.rdMolDescriptors import CalcNumBridgeheadAtoms as CalcNumBridgeheadAtoms
from rdkit.Chem.rdMolDescriptors import CalcNumHBA as CalcNumHBA
from rdkit.Chem.rdMolDescriptors import CalcNumHBD as CalcNumHBD
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms as CalcNumHeavyAtoms
from rdkit.Chem.rdMolDescriptors import CalcNumHeteroatoms as CalcNumHeteroatoms
from rdkit.Chem.rdMolDescriptors import CalcNumHeterocycles as CalcNumHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumLipinskiHBA as CalcNumLipinskiHBA
from rdkit.Chem.rdMolDescriptors import CalcNumLipinskiHBD as CalcNumLipinskiHBD
from rdkit.Chem.rdMolDescriptors import CalcNumRings as CalcNumRings
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds as CalcNumRotatableBonds
from rdkit.Chem.rdMolDescriptors import CalcNumSaturatedCarbocycles as CalcNumSaturatedCarbocycles
from rdkit.Chem.rdMolDescriptors import CalcNumSaturatedHeterocycles as CalcNumSaturatedHeterocycles
from rdkit.Chem.rdMolDescriptors import CalcNumSaturatedRings as CalcNumSaturatedRings
from rdkit.Chem.rdMolDescriptors import CalcNumSpiroAtoms as CalcNumSpiroAtoms
from rdkit.Chem.rdMolDescriptors import CalcNumUnspecifiedAtomStereoCenters as CalcNumUnspecifiedAtomStereoCenters
from rdkit.Chem.rdMolDescriptors import CalcOxidationNumbers as CalcOxidationNumbers
from rdkit.Chem.rdMolDescriptors import CalcPBF as CalcPBF
from rdkit.Chem.rdMolDescriptors import CalcPMI1 as CalcPMI1
from rdkit.Chem.rdMolDescriptors import CalcPMI2 as CalcPMI2
from rdkit.Chem.rdMolDescriptors import CalcPMI3 as CalcPMI3
from rdkit.Chem.rdMolDescriptors import CalcPhi as CalcPhi
from rdkit.Chem.rdMolDescriptors import CalcRDF as CalcRDF
from rdkit.Chem.rdMolAlign import CalcRMS as CalcRMS
from rdkit.Chem.rdMolDescriptors import CalcRadiusOfGyration as CalcRadiusOfGyration
from rdkit.Chem.rdMolDescriptors import CalcSpherocityIndex as CalcSpherocityIndex
from rdkit.Chem.rdMolDescriptors import CalcTPSA as CalcTPSA
from rdkit.Chem.rdMolDescriptors import CalcWHIM as CalcWHIM
from rdkit.Chem.rdmolfiles import CanonicalRankAtoms as CanonicalRankAtoms
from rdkit.Chem.rdmolfiles import CanonicalRankAtomsInFragment as CanonicalRankAtomsInFragment
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer as CanonicalizeConformer
from rdkit.Chem.rdmolfiles import CanonicalizeEnhancedStereo as CanonicalizeEnhancedStereo
from rdkit.Chem.rdMolTransforms import CanonicalizeMol as CanonicalizeMol
from rdkit.Chem.rdmolops import CanonicalizeStereoGroups as CanonicalizeStereoGroups
from rdkit.Chem.rdmolops import Cleanup as Cleanup
from rdkit.Chem.rdmolops import CleanupAtropisomers as CleanupAtropisomers
from rdkit.Chem.rdmolops import CleanupChirality as CleanupChirality
from rdkit.Chem.rdmolops import CleanupOrganometallics as CleanupOrganometallics
from rdkit.Chem.rdmolops import CleanupStereoGroups as CleanupStereoGroups
from rdkit.Chem.rdchem import ClearMolSubstanceGroups as ClearMolSubstanceGroups
from rdkit.Chem.rdmolops import CollapseAttachmentPoints as CollapseAttachmentPoints
from rdkit.Chem.rdmolops import CombineMols as CombineMols
from rdkit.Chem.rdDepictor import Compute2DCoords as Compute2DCoords
from rdkit.Chem.rdChemReactions import Compute2DCoordsForReaction as Compute2DCoordsForReaction
from rdkit.Chem.rdDepictor import Compute2DCoordsMimicDistmat as Compute2DCoordsMimicDistmat
from rdkit.Chem.rdmolops import ComputeAtomCIPRanks as ComputeAtomCIPRanks
from rdkit.Chem.rdMolTransforms import ComputeCanonicalTransform as ComputeCanonicalTransform
from rdkit.Chem.rdMolTransforms import ComputeCentroid as ComputeCentroid
from rdkit.Chem.rdShapeHelpers import ComputeConfBox as ComputeConfBox
from rdkit.Chem.rdShapeHelpers import ComputeConfDimsAndOffset as ComputeConfDimsAndOffset
from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges as ComputeGasteigerCharges
from rdkit.Chem.rdMolTransforms import ComputePrincipalAxesAndMoments as ComputePrincipalAxesAndMoments
from rdkit.Chem.rdMolTransforms import ComputePrincipalAxesAndMomentsFromGyrationMatrix as ComputePrincipalAxesAndMomentsFromGyrationMatrix
from rdkit.Chem.rdShapeHelpers import ComputeUnionBox as ComputeUnionBox
from rdkit.Chem.rdmolops import ConvertGenericQueriesToSubstanceGroups as ConvertGenericQueriesToSubstanceGroups
from rdkit.Chem.rdmolops import CopyMolSubset as CopyMolSubset
from rdkit.Chem.rdmolops import CountAtomElec as CountAtomElec
from rdkit.Chem.rdmolfiles import CreateAtomBoolPropertyList as CreateAtomBoolPropertyList
from rdkit.Chem.rdmolfiles import CreateAtomDoublePropertyList as CreateAtomDoublePropertyList
from rdkit.Chem.rdmolfiles import CreateAtomIntPropertyList as CreateAtomIntPropertyList
from rdkit.Chem.rdmolfiles import CreateAtomStringPropertyList as CreateAtomStringPropertyList
from rdkit.Chem.rdmolfiles import CreateBondBoolPropertyList as CreateBondBoolPropertyList
from rdkit.Chem.rdmolfiles import CreateBondDoublePropertyList as CreateBondDoublePropertyList
from rdkit.Chem.rdmolfiles import CreateBondIntPropertyList as CreateBondIntPropertyList
from rdkit.Chem.rdmolfiles import CreateBondStringPropertyList as CreateBondStringPropertyList
from rdkit.Chem.rdChemReactions import CreateDifferenceFingerprintForReaction as CreateDifferenceFingerprintForReaction
from rdkit.Chem.rdForceFieldHelpers import CreateEmptyForceFieldForMol as CreateEmptyForceFieldForMol
from rdkit.Chem.rdchem import CreateMolDataSubstanceGroup as CreateMolDataSubstanceGroup
from rdkit.Chem.rdchem import CreateMolSubstanceGroup as CreateMolSubstanceGroup
from rdkit.Chem.rdchem import CreateStereoGroup as CreateStereoGroup
from rdkit.Chem.rdChemReactions import CreateStructuralFingerprintForReaction as CreateStructuralFingerprintForReaction
from rdkit.Chem.rdMolDescriptors import CustomProp_VSA_ as CustomProp_VSA_
from rdkit.Chem.rdDistGeom import DG as DG
from rdkit.Chem.rdmolops import DativeBondsToHaptic as DativeBondsToHaptic
from rdkit.Chem.rdmolops import DeleteSubstructs as DeleteSubstructs
from rdkit.Chem.rdmolops import DetectBondStereoChemistry as DetectBondStereoChemistry
from rdkit.Chem.rdmolops import DetectBondStereochemistry as DetectBondStereochemistry
from rdkit.Chem.rdmolops import DetectChemistryProblems as DetectChemistryProblems
from rdkit.Chem.rdDistGeom import ETDG as ETDG
from rdkit.Chem.rdDistGeom import ETDGv2 as ETDGv2
from rdkit.Chem.rdDistGeom import ETKDG as ETKDG
from rdkit.Chem.rdDistGeom import ETKDGv2 as ETKDGv2
from rdkit.Chem.rdDistGeom import ETKDGv3 as ETKDGv3
from rdkit.Chem.rdDistGeom import EmbedMolecule as EmbedMolecule
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs as EmbedMultipleConfs
from rdkit.Chem.rdDistGeom import EmbedParametersToJSON as EmbedParametersToJSON
from rdkit.Chem.rdShapeHelpers import EncodeShape as EncodeShape
from rdkit.Chem.rdMolEnumerator import Enumerate as Enumerate
from rdkit.Chem.rdChemReactions import EnumerateLibraryCanSerialize as EnumerateLibraryCanSerialize
from rdkit.Chem.rdmolops import ExpandAttachmentPoints as ExpandAttachmentPoints
from rdkit.Chem.rdqueries import ExplicitDegreeEqualsQueryAtom as ExplicitDegreeEqualsQueryAtom
from rdkit.Chem.rdqueries import ExplicitDegreeGreaterQueryAtom as ExplicitDegreeGreaterQueryAtom
from rdkit.Chem.rdqueries import ExplicitDegreeLessQueryAtom as ExplicitDegreeLessQueryAtom
from rdkit.Chem.rdqueries import ExplicitValenceEqualsQueryAtom as ExplicitValenceEqualsQueryAtom
from rdkit.Chem.rdqueries import ExplicitValenceGreaterQueryAtom as ExplicitValenceGreaterQueryAtom
from rdkit.Chem.rdqueries import ExplicitValenceLessQueryAtom as ExplicitValenceLessQueryAtom
from rdkit.Chem.rdmolops import FastFindRings as FastFindRings
from rdkit.Chem.rdmolops import FindAllPathsOfLengthN as FindAllPathsOfLengthN
from rdkit.Chem.rdmolops import FindAllSubgraphsOfLengthMToN as FindAllSubgraphsOfLengthMToN
from rdkit.Chem.rdmolops import FindAllSubgraphsOfLengthN as FindAllSubgraphsOfLengthN
from rdkit.Chem.rdmolops import FindAtomEnvironmentOfRadiusN as FindAtomEnvironmentOfRadiusN
from rdkit.Chem.rdmolops import FindMesoCenters as FindMesoCenters
from rdkit.Chem.rdmolops import FindPotentialStereo as FindPotentialStereo
from rdkit.Chem.rdmolops import FindPotentialStereoBonds as FindPotentialStereoBonds
from rdkit.Chem.rdmolops import FindRingFamilies as FindRingFamilies
from rdkit.Chem.rdmolops import FindUniqueSubgraphsOfLengthN as FindUniqueSubgraphsOfLengthN
from rdkit.Chem.rdFingerprintGenerator import FingerprintGeneratorFromJSON as FingerprintGeneratorFromJSON
from rdkit.Chem.rdqueries import FormalChargeEqualsQueryAtom as FormalChargeEqualsQueryAtom
from rdkit.Chem.rdqueries import FormalChargeGreaterQueryAtom as FormalChargeGreaterQueryAtom
from rdkit.Chem.rdqueries import FormalChargeLessQueryAtom as FormalChargeLessQueryAtom
from rdkit.Chem.rdchem import ForwardStereoGroupIds as ForwardStereoGroupIds
from rdkit.Chem.rdmolops import FragmentOnBRICSBonds as FragmentOnBRICSBonds
from rdkit.Chem.rdmolops import FragmentOnBonds as FragmentOnBonds
from rdkit.Chem.rdmolops import FragmentOnSomeBonds as FragmentOnSomeBonds
from rdkit.Chem.rdDepictor import GenerateDepictionMatching2DStructure as GenerateDepictionMatching2DStructure
from rdkit.Chem.rdDepictor import GenerateDepictionMatching3DStructure as GenerateDepictionMatching3DStructure
from rdkit.Chem.rdReducedGraphs import GenerateErGFingerprintForReducedGraph as GenerateErGFingerprintForReducedGraph
from rdkit.Chem.rdReducedGraphs import GenerateMolExtendedReducedGraph as GenerateMolExtendedReducedGraph
from rdkit.Chem.rdmolops import Get3DDistanceMatrix as Get3DDistanceMatrix
from rdkit.Chem.rdmolops import GetAdjacencyMatrix as GetAdjacencyMatrix
from rdkit.Chem.rdMolAlign import GetAlignmentTransform as GetAlignmentTransform
from rdkit.Chem.rdMolAlign import GetAllConformerBestRMS as GetAllConformerBestRMS
from rdkit.Chem.rdmolops import GetAllowNontetrahedralChirality as GetAllowNontetrahedralChirality
from rdkit.Chem.rdMolTransforms import GetAngleDeg as GetAngleDeg
from rdkit.Chem.rdMolTransforms import GetAngleRad as GetAngleRad
from rdkit.Chem.rdchem import GetAtomAlias as GetAtomAlias
from rdkit.Chem.rdMolDescriptors import GetAtomFeatures as GetAtomFeatures
from rdkit.Chem.rdMolChemicalFeatures import GetAtomMatch as GetAtomMatch
from rdkit.Chem.rdMolDescriptors import GetAtomPairAtomCode as GetAtomPairAtomCode
from rdkit.Chem.rdFingerprintGenerator import GetAtomPairAtomInvGen as GetAtomPairAtomInvGen
from rdkit.Chem.rdMolDescriptors import GetAtomPairCode as GetAtomPairCode
from rdkit.Chem.rdMolDescriptors import GetAtomPairFingerprint as GetAtomPairFingerprint
from rdkit.Chem.rdFingerprintGenerator import GetAtomPairGenerator as GetAtomPairGenerator
from rdkit.Chem.rdchem import GetAtomRLabel as GetAtomRLabel
from rdkit.Chem.rdchem import GetAtomValue as GetAtomValue
from rdkit.Chem.rdMolAlign import GetBestAlignmentTransform as GetBestAlignmentTransform
from rdkit.Chem.rdMolAlign import GetBestRMS as GetBestRMS
from rdkit.Chem.rdMolTransforms import GetBondLength as GetBondLength
from rdkit.Chem.rdChemReactions import GetChemDrawRxnAdjustParams as GetChemDrawRxnAdjustParams
from rdkit.Chem.rdMolDescriptors import GetConnectivityInvariants as GetConnectivityInvariants
from rdkit.Chem.rdFingerprintGenerator import GetCountFPs as GetCountFPs
from rdkit.Chem.rdMolAlign import GetCrippenO3A as GetCrippenO3A
from rdkit.Chem.rdMolAlign import GetCrippenO3AForProbeConfs as GetCrippenO3AForProbeConfs
from rdkit.Chem.rdChemReactions import GetDefaultAdjustParams as GetDefaultAdjustParams
from rdkit.Chem.rdchem import GetDefaultPickleProperties as GetDefaultPickleProperties
from rdkit.Chem.rdMolTransforms import GetDihedralDeg as GetDihedralDeg
from rdkit.Chem.rdMolTransforms import GetDihedralRad as GetDihedralRad
from rdkit.Chem.rdmolops import GetDistanceMatrix as GetDistanceMatrix
from rdkit.Chem.rdReducedGraphs import GetErGFingerprint as GetErGFingerprint
from rdkit.Chem.rdDistGeom import GetExperimentalTorsions as GetExperimentalTorsions
from rdkit.Chem.rdFingerprintGenerator import GetFPs as GetFPs
from rdkit.Chem.rdMolDescriptors import GetFeatureInvariants as GetFeatureInvariants
from rdkit.Chem.rdmolops import GetFormalCharge as GetFormalCharge
from rdkit.Chem.rdMolDescriptors import GetHashedAtomPairFingerprint as GetHashedAtomPairFingerprint
from rdkit.Chem.rdMolDescriptors import GetHashedAtomPairFingerprintAsBitVect as GetHashedAtomPairFingerprintAsBitVect
from rdkit.Chem.rdMolDescriptors import GetHashedMorganFingerprint as GetHashedMorganFingerprint
from rdkit.Chem.rdMolDescriptors import GetHashedTopologicalTorsionFingerprint as GetHashedTopologicalTorsionFingerprint
from rdkit.Chem.rdMolDescriptors import GetHashedTopologicalTorsionFingerprintAsBitVect as GetHashedTopologicalTorsionFingerprintAsBitVect
from rdkit.Chem.rdinchi import GetInchiVersion as GetInchiVersion
from rdkit.Chem.rdMolDescriptors import GetMACCSKeysFingerprint as GetMACCSKeysFingerprint
from rdkit.Chem.rdmolops import GetMolFrags as GetMolFrags
from rdkit.Chem.rdchem import GetMolSubstanceGroupWithIdx as GetMolSubstanceGroupWithIdx
from rdkit.Chem.rdchem import GetMolSubstanceGroups as GetMolSubstanceGroups
from rdkit.Chem.rdDistGeom import GetMoleculeBoundsMatrix as GetMoleculeBoundsMatrix
from rdkit.Chem.rdFingerprintGenerator import GetMorganAtomInvGen as GetMorganAtomInvGen
from rdkit.Chem.rdFingerprintGenerator import GetMorganBondInvGen as GetMorganBondInvGen
from rdkit.Chem.rdFingerprintGenerator import GetMorganFeatureAtomInvGen as GetMorganFeatureAtomInvGen
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint as GetMorganFingerprint
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect as GetMorganFingerprintAsBitVect
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator as GetMorganGenerator
from rdkit.Chem.rdmolops import GetMostSubstitutedCoreMatch as GetMostSubstitutedCoreMatch
from rdkit.Chem.rdchem import GetNumPiElectrons as GetNumPiElectrons
from rdkit.Chem.rdMolAlign import GetO3A as GetO3A
from rdkit.Chem.rdMolAlign import GetO3AForProbeConfs as GetO3AForProbeConfs
from rdkit.Chem.rdchem import GetPeriodicTable as GetPeriodicTable
from rdkit.Chem.rdDepictor import GetPreferCoordGen as GetPreferCoordGen
from rdkit.Chem.rdFingerprintGenerator import GetRDKitAtomInvGen as GetRDKitAtomInvGen
from rdkit.Chem.rdFingerprintGenerator import GetRDKitFPGenerator as GetRDKitFPGenerator
from rdkit.Chem.rdmolops import GetSSSR as GetSSSR
from rdkit.Chem.rdmolops import GetShortestPath as GetShortestPath
from rdkit.Chem.rdFingerprintGenerator import GetSparseCountFPs as GetSparseCountFPs
from rdkit.Chem.rdFingerprintGenerator import GetSparseFPs as GetSparseFPs
from rdkit.Chem.rdchem import GetSupplementalSmilesLabel as GetSupplementalSmilesLabel
from rdkit.Chem.rdmolops import GetSymmSSSR as GetSymmSSSR
from rdkit.Chem.rdMolDescriptors import GetTopologicalTorsionFingerprint as GetTopologicalTorsionFingerprint
from rdkit.Chem.rdFingerprintGenerator import GetTopologicalTorsionGenerator as GetTopologicalTorsionGenerator
from rdkit.Chem.rdForceFieldHelpers import GetUFFAngleBendParams as GetUFFAngleBendParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFBondStretchParams as GetUFFBondStretchParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFInversionParams as GetUFFInversionParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFTorsionParams as GetUFFTorsionParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFVdWParams as GetUFFVdWParams
from rdkit.Chem.rdMolDescriptors import GetUSR as GetUSR
from rdkit.Chem.rdMolDescriptors import GetUSRCAT as GetUSRCAT
from rdkit.Chem.rdMolDescriptors import GetUSRDistributions as GetUSRDistributions
from rdkit.Chem.rdMolDescriptors import GetUSRDistributionsFromPoints as GetUSRDistributionsFromPoints
from rdkit.Chem.rdMolDescriptors import GetUSRFromDistributions as GetUSRFromDistributions
from rdkit.Chem.rdMolDescriptors import GetUSRScore as GetUSRScore
from rdkit.Chem.rdmolops import GetUseLegacyStereoPerception as GetUseLegacyStereoPerception
from rdkit.Chem.rdqueries import HCountEqualsQueryAtom as HCountEqualsQueryAtom
from rdkit.Chem.rdqueries import HCountGreaterQueryAtom as HCountGreaterQueryAtom
from rdkit.Chem.rdqueries import HCountLessQueryAtom as HCountLessQueryAtom
from rdkit.Chem.rdmolops import HapticBondsToDative as HapticBondsToDative
from rdkit.Chem.rdChemReactions import HasAgentTemplateSubstructMatch as HasAgentTemplateSubstructMatch
from rdkit.Chem.rdqueries import HasBitVectPropWithValueQueryAtom as HasBitVectPropWithValueQueryAtom
from rdkit.Chem.rdqueries import HasBoolPropWithValueQueryAtom as HasBoolPropWithValueQueryAtom
from rdkit.Chem.rdqueries import HasBoolPropWithValueQueryBond as HasBoolPropWithValueQueryBond
from rdkit.Chem.rdmolfiles import HasChemDrawCDXSupport as HasChemDrawCDXSupport
from rdkit.Chem.rdqueries import HasChiralTagQueryAtom as HasChiralTagQueryAtom
from rdkit.Chem.rdqueries import HasDoublePropWithValueQueryAtom as HasDoublePropWithValueQueryAtom
from rdkit.Chem.rdqueries import HasDoublePropWithValueQueryBond as HasDoublePropWithValueQueryBond
from rdkit.Chem.rdqueries import HasIntPropWithValueQueryAtom as HasIntPropWithValueQueryAtom
from rdkit.Chem.rdqueries import HasIntPropWithValueQueryBond as HasIntPropWithValueQueryBond
from rdkit.Chem.rdChemReactions import HasProductTemplateSubstructMatch as HasProductTemplateSubstructMatch
from rdkit.Chem.rdqueries import HasPropQueryAtom as HasPropQueryAtom
from rdkit.Chem.rdqueries import HasPropQueryBond as HasPropQueryBond
from rdkit.Chem.rdmolops import HasQueryHs as HasQueryHs
from rdkit.Chem.rdChemReactions import HasReactantTemplateSubstructMatch as HasReactantTemplateSubstructMatch
from rdkit.Chem.rdChemReactions import HasReactionAtomMapping as HasReactionAtomMapping
from rdkit.Chem.rdChemReactions import HasReactionSubstructMatch as HasReactionSubstructMatch
from rdkit.Chem.rdqueries import HasStringPropWithValueQueryAtom as HasStringPropWithValueQueryAtom
from rdkit.Chem.rdqueries import HasStringPropWithValueQueryBond as HasStringPropWithValueQueryBond
from rdkit.Chem.rdqueries import HybridizationEqualsQueryAtom as HybridizationEqualsQueryAtom
from rdkit.Chem.rdqueries import HybridizationGreaterQueryAtom as HybridizationGreaterQueryAtom
from rdkit.Chem.rdqueries import HybridizationLessQueryAtom as HybridizationLessQueryAtom
from rdkit.Chem.rdqueries import InNRingsEqualsQueryAtom as InNRingsEqualsQueryAtom
from rdkit.Chem.rdqueries import InNRingsGreaterQueryAtom as InNRingsGreaterQueryAtom
from rdkit.Chem.rdqueries import InNRingsLessQueryAtom as InNRingsLessQueryAtom
from rdkit.Chem.rdqueries import IsAliphaticQueryAtom as IsAliphaticQueryAtom
from rdkit.Chem.rdqueries import IsAromaticQueryAtom as IsAromaticQueryAtom
from rdkit.Chem.rdqueries import IsBridgeheadQueryAtom as IsBridgeheadQueryAtom
from rdkit.Chem.rdDepictor import IsCoordGenSupportAvailable as IsCoordGenSupportAvailable
from rdkit.Chem.rdqueries import IsInRingQueryAtom as IsInRingQueryAtom
from rdkit.Chem.rdChemReactions import IsReactionTemplateMoleculeAgent as IsReactionTemplateMoleculeAgent
from rdkit.Chem.rdqueries import IsUnsaturatedQueryAtom as IsUnsaturatedQueryAtom
from rdkit.Chem.rdqueries import IsotopeEqualsQueryAtom as IsotopeEqualsQueryAtom
from rdkit.Chem.rdqueries import IsotopeGreaterQueryAtom as IsotopeGreaterQueryAtom
from rdkit.Chem.rdqueries import IsotopeLessQueryAtom as IsotopeLessQueryAtom
from rdkit.Chem.rdMolInterchange import JSONToMols as JSONToMols
from rdkit.Chem.rdDistGeom import KDG as KDG
from rdkit.Chem.rdmolops import Kekulize as Kekulize
from rdkit.Chem.rdmolops import KekulizeIfPossible as KekulizeIfPossible
from rdkit.Chem.rdmolops import LayeredFingerprint as LayeredFingerprint
from rdkit.Chem.rdDepictor import LoadDefaultRingSystemTemplates as LoadDefaultRingSystemTemplates
from rdkit.Chem.rdqueries import MAtomQueryAtom as MAtomQueryAtom
from rdkit.Chem.rdqueries import MHAtomQueryAtom as MHAtomQueryAtom
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField as MMFFGetMoleculeForceField
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties as MMFFGetMoleculeProperties
from rdkit.Chem.rdForceFieldHelpers import MMFFHasAllMoleculeParams as MMFFHasAllMoleculeParams
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule as MMFFOptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs as MMFFOptimizeMoleculeConfs
from rdkit.Chem.rdForceFieldHelpers import MMFFSanitizeMolecule as MMFFSanitizeMolecule
from rdkit.Chem.rdMolDescriptors import MQNs_ as MQNs_
from rdkit.Chem.rdMolDescriptors import MakePropertyRangeQuery as MakePropertyRangeQuery
from rdkit.Chem.rdqueries import MassEqualsQueryAtom as MassEqualsQueryAtom
from rdkit.Chem.rdqueries import MassGreaterQueryAtom as MassGreaterQueryAtom
from rdkit.Chem.rdqueries import MassLessQueryAtom as MassLessQueryAtom
from rdkit.Chem.rdChemReactions import MatchOnlyAtRgroupsAdjustParams as MatchOnlyAtRgroupsAdjustParams
from rdkit.Chem.rdmolops import MergeQueryHs as MergeQueryHs
from rdkit.Chem.rdmolfiles import MetadataFromPNGFile as MetadataFromPNGFile
from rdkit.Chem.rdmolfiles import MetadataFromPNGString as MetadataFromPNGString
from rdkit.Chem.rdqueries import MinRingSizeEqualsQueryAtom as MinRingSizeEqualsQueryAtom
from rdkit.Chem.rdqueries import MinRingSizeGreaterQueryAtom as MinRingSizeGreaterQueryAtom
from rdkit.Chem.rdqueries import MinRingSizeLessQueryAtom as MinRingSizeLessQueryAtom
from rdkit.Chem.rdqueries import MissingChiralTagQueryAtom as MissingChiralTagQueryAtom
from rdkit.Chem.rdmolops import MolAddRecursiveQueries as MolAddRecursiveQueries
from rdkit.Chem.rdchem import MolBundleCanSerialize as MolBundleCanSerialize
from rdkit.Chem.rdmolfiles import MolFragmentToCXSmarts as MolFragmentToCXSmarts
from rdkit.Chem.rdmolfiles import MolFragmentToCXSmiles as MolFragmentToCXSmiles
from rdkit.Chem.rdmolfiles import MolFragmentToSmarts as MolFragmentToSmarts
from rdkit.Chem.rdmolfiles import MolFragmentToSmiles as MolFragmentToSmiles
from rdkit.Chem.rdmolfiles import MolFromFASTA as MolFromFASTA
from rdkit.Chem.rdmolfiles import MolFromHELM as MolFromHELM
from rdkit.Chem.rdmolfiles import MolFromMol2Block as MolFromMol2Block
from rdkit.Chem.rdmolfiles import MolFromMol2File as MolFromMol2File
from rdkit.Chem.rdmolfiles import MolFromMolBlock as MolFromMolBlock
from rdkit.Chem.rdmolfiles import MolFromMolFile as MolFromMolFile
from rdkit.Chem.rdmolfiles import MolFromMrvBlock as MolFromMrvBlock
from rdkit.Chem.rdmolfiles import MolFromMrvFile as MolFromMrvFile
from rdkit.Chem.rdmolfiles import MolFromPDBBlock as MolFromPDBBlock
from rdkit.Chem.rdmolfiles import MolFromPDBFile as MolFromPDBFile
from rdkit.Chem.rdmolfiles import MolFromPNGFile as MolFromPNGFile
from rdkit.Chem.rdmolfiles import MolFromPNGString as MolFromPNGString
from rdkit.Chem.rdSLNParse import MolFromQuerySLN as MolFromQuerySLN
from rdkit.Chem.rdmolfiles import MolFromRDKitSVG as MolFromRDKitSVG
from rdkit.Chem.rdmolfiles import MolFromSCSRBlock as MolFromSCSRBlock
from rdkit.Chem.rdmolfiles import MolFromSCSRFile as MolFromSCSRFile
from rdkit.Chem.rdSLNParse import MolFromSLN as MolFromSLN
from rdkit.Chem.rdmolfiles import MolFromSequence as MolFromSequence
from rdkit.Chem.rdmolfiles import MolFromSmarts as MolFromSmarts
from rdkit.Chem.rdmolfiles import MolFromSmiles as MolFromSmiles
from rdkit.Chem.rdmolfiles import MolFromTPLBlock as MolFromTPLBlock
from rdkit.Chem.rdmolfiles import MolFromTPLFile as MolFromTPLFile
from rdkit.Chem.rdmolfiles import MolFromXYZBlock as MolFromXYZBlock
from rdkit.Chem.rdmolfiles import MolFromXYZFile as MolFromXYZFile
from rdkit.Chem.rdmolfiles import MolMetadataToPNGFile as MolMetadataToPNGFile
from rdkit.Chem.rdmolfiles import MolMetadataToPNGString as MolMetadataToPNGString
from rdkit.Chem.rdmolfiles import MolToCDXMLBlock as MolToCDXMLBlock
from rdkit.Chem.rdmolfiles import MolToCMLBlock as MolToCMLBlock
from rdkit.Chem.rdmolfiles import MolToCMLFile as MolToCMLFile
from rdkit.Chem.rdmolfiles import MolToCXSmarts as MolToCXSmarts
from rdkit.Chem.rdmolfiles import MolToCXSmiles as MolToCXSmiles
from rdkit.Chem.rdmolfiles import MolToFASTA as MolToFASTA
from rdkit.Chem.rdmolfiles import MolToHELM as MolToHELM
from rdkit.Chem.rdMolInterchange import MolToJSON as MolToJSON
from rdkit.Chem.rdmolfiles import MolToMolBlock as MolToMolBlock
from rdkit.Chem.rdmolfiles import MolToMolFile as MolToMolFile
from rdkit.Chem.rdmolfiles import MolToMrvBlock as MolToMrvBlock
from rdkit.Chem.rdmolfiles import MolToMrvFile as MolToMrvFile
from rdkit.Chem.rdmolfiles import MolToPDBBlock as MolToPDBBlock
from rdkit.Chem.rdmolfiles import MolToPDBFile as MolToPDBFile
from rdkit.Chem.rdmolfiles import MolToRandomSmilesVect as MolToRandomSmilesVect
from rdkit.Chem.rdmolfiles import MolToSequence as MolToSequence
from rdkit.Chem.rdmolfiles import MolToSmarts as MolToSmarts
from rdkit.Chem.rdmolfiles import MolToSmiles as MolToSmiles
from rdkit.Chem.rdmolfiles import MolToTPLBlock as MolToTPLBlock
from rdkit.Chem.rdmolfiles import MolToTPLFile as MolToTPLFile
from rdkit.Chem.rdmolfiles import MolToV2KMolBlock as MolToV2KMolBlock
from rdkit.Chem.rdmolfiles import MolToV3KMolBlock as MolToV3KMolBlock
from rdkit.Chem.rdmolfiles import MolToV3KMolFile as MolToV3KMolFile
from rdkit.Chem.rdmolfiles import MolToXYZBlock as MolToXYZBlock
from rdkit.Chem.rdmolfiles import MolToXYZFile as MolToXYZFile
from rdkit.Chem.rdmolfiles import MolsFromCDXML as MolsFromCDXML
from rdkit.Chem.rdmolfiles import MolsFromCDXMLFile as MolsFromCDXMLFile
from rdkit.Chem.rdmolfiles import MolsFromPNGFile as MolsFromPNGFile
from rdkit.Chem.rdmolfiles import MolsFromPNGString as MolsFromPNGString
from rdkit.Chem.rdMolInterchange import MolsToJSON as MolsToJSON
from rdkit.Chem.rdChemReactions import MrvBlockIsReaction as MrvBlockIsReaction
from rdkit.Chem.rdChemReactions import MrvFileIsReaction as MrvFileIsReaction
from rdkit.Chem.rdmolops import MurckoDecompose as MurckoDecompose
from rdkit.Chem.rdmolops import NeedsHs as NeedsHs
from rdkit.Chem.rdqueries import NonHydrogenDegreeEqualsQueryAtom as NonHydrogenDegreeEqualsQueryAtom
from rdkit.Chem.rdqueries import NonHydrogenDegreeGreaterQueryAtom as NonHydrogenDegreeGreaterQueryAtom
from rdkit.Chem.rdqueries import NonHydrogenDegreeLessQueryAtom as NonHydrogenDegreeLessQueryAtom
from rdkit.Chem.rdDepictor import NormalizeDepiction as NormalizeDepiction
from rdkit.Chem.rdqueries import NumAliphaticHeteroatomNeighborsEqualsQueryAtom as NumAliphaticHeteroatomNeighborsEqualsQueryAtom
from rdkit.Chem.rdqueries import NumAliphaticHeteroatomNeighborsGreaterQueryAtom as NumAliphaticHeteroatomNeighborsGreaterQueryAtom
from rdkit.Chem.rdqueries import NumAliphaticHeteroatomNeighborsLessQueryAtom as NumAliphaticHeteroatomNeighborsLessQueryAtom
from rdkit.Chem.rdqueries import NumHeteroatomNeighborsEqualsQueryAtom as NumHeteroatomNeighborsEqualsQueryAtom
from rdkit.Chem.rdqueries import NumHeteroatomNeighborsGreaterQueryAtom as NumHeteroatomNeighborsGreaterQueryAtom
from rdkit.Chem.rdqueries import NumHeteroatomNeighborsLessQueryAtom as NumHeteroatomNeighborsLessQueryAtom
from rdkit.Chem.rdqueries import NumRadicalElectronsEqualsQueryAtom as NumRadicalElectronsEqualsQueryAtom
from rdkit.Chem.rdqueries import NumRadicalElectronsGreaterQueryAtom as NumRadicalElectronsGreaterQueryAtom
from rdkit.Chem.rdqueries import NumRadicalElectronsLessQueryAtom as NumRadicalElectronsLessQueryAtom
from rdkit.Chem.rdForceFieldHelpers import OptimizeMolecule as OptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import OptimizeMoleculeConfs as OptimizeMoleculeConfs
from rdkit.Chem.rdDistGeom import OrderedEmbedFailureCauses as OrderedEmbedFailureCauses
from rdkit.Chem.rdMolDescriptors import PEOE_VSA_ as PEOE_VSA_
from rdkit.Chem.rdmolops import ParseMolQueryDefFile as ParseMolQueryDefFile
from rdkit.Chem.rdmolops import PathToSubmol as PathToSubmol
from rdkit.Chem.rdmolops import PatternFingerprint as PatternFingerprint
from rdkit.Chem.rdChemReactions import PreprocessReaction as PreprocessReaction
from rdkit.Chem.rdqueries import QAtomQueryAtom as QAtomQueryAtom
from rdkit.Chem.rdqueries import QHAtomQueryAtom as QHAtomQueryAtom
from rdkit.Chem.rdmolops import RDKFingerprint as RDKFingerprint
from rdkit.Chem.rdMolAlign import RandomTransform as RandomTransform
from rdkit.Chem.rdChemReactions import ReactionFromMolecule as ReactionFromMolecule
from rdkit.Chem.rdChemReactions import ReactionFromMrvBlock as ReactionFromMrvBlock
from rdkit.Chem.rdChemReactions import ReactionFromMrvFile as ReactionFromMrvFile
from rdkit.Chem.rdChemReactions import ReactionFromPNGFile as ReactionFromPNGFile
from rdkit.Chem.rdChemReactions import ReactionFromPNGString as ReactionFromPNGString
from rdkit.Chem.rdChemReactions import ReactionFromRxnBlock as ReactionFromRxnBlock
from rdkit.Chem.rdChemReactions import ReactionFromRxnFile as ReactionFromRxnFile
from rdkit.Chem.rdChemReactions import ReactionFromSmarts as ReactionFromSmarts
from rdkit.Chem.rdChemReactions import ReactionFromSmiles as ReactionFromSmiles
from rdkit.Chem.rdChemReactions import ReactionMetadataToPNGFile as ReactionMetadataToPNGFile
from rdkit.Chem.rdChemReactions import ReactionMetadataToPNGString as ReactionMetadataToPNGString
from rdkit.Chem.rdChemReactions import ReactionToCXSmarts as ReactionToCXSmarts
from rdkit.Chem.rdChemReactions import ReactionToCXSmiles as ReactionToCXSmiles
from rdkit.Chem.rdChemReactions import ReactionToMolecule as ReactionToMolecule
from rdkit.Chem.rdChemReactions import ReactionToMrvBlock as ReactionToMrvBlock
from rdkit.Chem.rdChemReactions import ReactionToMrvFile as ReactionToMrvFile
from rdkit.Chem.rdChemReactions import ReactionToRxnBlock as ReactionToRxnBlock
from rdkit.Chem.rdChemReactions import ReactionToSmarts as ReactionToSmarts
from rdkit.Chem.rdChemReactions import ReactionToSmiles as ReactionToSmiles
from rdkit.Chem.rdChemReactions import ReactionToV3KRxnBlock as ReactionToV3KRxnBlock
from rdkit.Chem.rdChemReactions import ReactionsFromCDXMLBlock as ReactionsFromCDXMLBlock
from rdkit.Chem.rdChemReactions import ReactionsFromCDXMLFile as ReactionsFromCDXMLFile
from rdkit.Chem.rdmolops import ReapplyMolBlockWedging as ReapplyMolBlockWedging
from rdkit.Chem.rdChemReactions import ReduceProductToSideChains as ReduceProductToSideChains
from rdkit.Chem.rdmolops import RemoveAllHs as RemoveAllHs
from rdkit.Chem.rdmolops import RemoveHs as RemoveHs
from rdkit.Chem.rdChemReactions import RemoveMappingNumbersFromReactions as RemoveMappingNumbersFromReactions
from rdkit.Chem.rdmolops import RemoveNonExplicit3DChirality as RemoveNonExplicit3DChirality
from rdkit.Chem.rdmolops import RemoveStereochemistry as RemoveStereochemistry
from rdkit.Chem.rdmolops import RenumberAtoms as RenumberAtoms
from rdkit.Chem.rdqueries import ReplaceAtomWithQueryAtom as ReplaceAtomWithQueryAtom
from rdkit.Chem.rdmolops import ReplaceCore as ReplaceCore
from rdkit.Chem.rdmolops import ReplaceSidechains as ReplaceSidechains
from rdkit.Chem.rdmolops import ReplaceSubstructs as ReplaceSubstructs
from rdkit.Chem.rdqueries import RingBondCountEqualsQueryAtom as RingBondCountEqualsQueryAtom
from rdkit.Chem.rdqueries import RingBondCountGreaterQueryAtom as RingBondCountGreaterQueryAtom
from rdkit.Chem.rdqueries import RingBondCountLessQueryAtom as RingBondCountLessQueryAtom
from rdkit.Chem.rdMolDescriptors import SMR_VSA_ as SMR_VSA_
from rdkit.Chem.rdmolops import SanitizeMol as SanitizeMol
from rdkit.Chem.rdChemReactions import SanitizeRxn as SanitizeRxn
from rdkit.Chem.rdChemReactions import SanitizeRxnAsMols as SanitizeRxnAsMols
from rdkit.Chem.rdmolops import SetAllowNontetrahedralChirality as SetAllowNontetrahedralChirality
from rdkit.Chem.rdMolTransforms import SetAngleDeg as SetAngleDeg
from rdkit.Chem.rdMolTransforms import SetAngleRad as SetAngleRad
from rdkit.Chem.rdmolops import SetAromaticity as SetAromaticity
from rdkit.Chem.rdchem import SetAtomAlias as SetAtomAlias
from rdkit.Chem.rdchem import SetAtomRLabel as SetAtomRLabel
from rdkit.Chem.rdchem import SetAtomValue as SetAtomValue
from rdkit.Chem.rdMolTransforms import SetBondLength as SetBondLength
from rdkit.Chem.rdmolops import SetBondStereoFromDirections as SetBondStereoFromDirections
from rdkit.Chem.rdmolops import SetConjugation as SetConjugation
from rdkit.Chem.rdchem import SetDefaultPickleProperties as SetDefaultPickleProperties
from rdkit.Chem.rdMolTransforms import SetDihedralDeg as SetDihedralDeg
from rdkit.Chem.rdMolTransforms import SetDihedralRad as SetDihedralRad
from rdkit.Chem.rdmolops import SetDoubleBondNeighborDirections as SetDoubleBondNeighborDirections
from rdkit.Chem.rdmolops import SetGenericQueriesFromProperties as SetGenericQueriesFromProperties
from rdkit.Chem.rdmolops import SetHybridization as SetHybridization
from rdkit.Chem.rdDepictor import SetPreferCoordGen as SetPreferCoordGen
from rdkit.Chem.rdDepictor import SetRingSystemTemplates as SetRingSystemTemplates
from rdkit.Chem.rdchem import SetSupplementalSmilesLabel as SetSupplementalSmilesLabel
from rdkit.Chem.rdmolops import SetTerminalAtomCoords as SetTerminalAtomCoords
from rdkit.Chem.rdmolops import SetUseLegacyStereoPerception as SetUseLegacyStereoPerception
from rdkit.Chem.rdShapeHelpers import ShapeProtrudeDist as ShapeProtrudeDist
from rdkit.Chem.rdShapeHelpers import ShapeTanimotoDist as ShapeTanimotoDist
from rdkit.Chem.rdShapeHelpers import ShapeTverskyIndex as ShapeTverskyIndex
from rdkit.Chem.rdmolops import SimplifyEnhancedStereo as SimplifyEnhancedStereo
from rdkit.Chem.rdMolDescriptors import SlogP_VSA_ as SlogP_VSA_
from rdkit.Chem.rdmolfiles import SmilesMolSupplierFromText as SmilesMolSupplierFromText
from rdkit.Chem.rdmolops import SortMatchesByDegreeOfCoreSubstitution as SortMatchesByDegreeOfCoreSubstitution
from rdkit.Chem.rdmolops import SplitMolByPDBChainId as SplitMolByPDBChainId
from rdkit.Chem.rdmolops import SplitMolByPDBResidues as SplitMolByPDBResidues
from rdkit.Chem.rdDepictor import StraightenDepiction as StraightenDepiction
from rdkit.Chem.rdqueries import TotalDegreeEqualsQueryAtom as TotalDegreeEqualsQueryAtom
from rdkit.Chem.rdqueries import TotalDegreeGreaterQueryAtom as TotalDegreeGreaterQueryAtom
from rdkit.Chem.rdqueries import TotalDegreeLessQueryAtom as TotalDegreeLessQueryAtom
from rdkit.Chem.rdqueries import TotalValenceEqualsQueryAtom as TotalValenceEqualsQueryAtom
from rdkit.Chem.rdqueries import TotalValenceGreaterQueryAtom as TotalValenceGreaterQueryAtom
from rdkit.Chem.rdqueries import TotalValenceLessQueryAtom as TotalValenceLessQueryAtom
from rdkit.Chem.rdMolTransforms import TransformConformer as TransformConformer
from rdkit.Chem.rdmolops import TranslateChiralFlagToStereoGroups as TranslateChiralFlagToStereoGroups
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField as UFFGetMoleculeForceField
from rdkit.Chem.rdForceFieldHelpers import UFFHasAllMoleculeParams as UFFHasAllMoleculeParams
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule as UFFOptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMoleculeConfs as UFFOptimizeMoleculeConfs
from rdkit.Chem.rdmolops import UnfoldedRDKFingerprintCountBased as UnfoldedRDKFingerprintCountBased
from rdkit.Chem.rdChemReactions import UpdateProductsStereochemistry as UpdateProductsStereochemistry
from rdkit.Chem.rdmolops import WedgeBond as WedgeBond
from rdkit.Chem.rdmolops import WedgeMolBonds as WedgeMolBonds
from rdkit.Chem.rdqueries import XAtomQueryAtom as XAtomQueryAtom
from rdkit.Chem.rdqueries import XHAtomQueryAtom as XHAtomQueryAtom
from rdkit.Chem.rdmolops import molzip as molzip
from rdkit.Chem.rdmolops import molzipFragments as molzipFragments
from rdkit.Chem.rdDistGeom import srETKDGv3 as srETKDGv3
from rdkit.Chem.rdchem import tossit as tossit
