# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation

"""
from __future__ import annotations
import typing
_RDKitObj = typing.TypeVar('_RDKitObj')
from rdkit.Chem.inchi import InchiReadWriteError
from rdkit.Chem.inchi import InchiToInchiKey
from rdkit.Chem.inchi import MolBlockToInchi
from rdkit.Chem.inchi import MolBlockToInchiAndAuxInfo
from rdkit.Chem.inchi import MolFromInchi
from rdkit.Chem.inchi import MolFromInchiAndAuxInfo
from rdkit.Chem.inchi import MolToInchi
from rdkit.Chem.inchi import MolToInchiAndAuxInfo
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdMolInterchange import JSONParseParameters
from rdkit.Chem.rdMolInterchange import JSONWriteParameters
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
from rdkit.Chem.rdmolops import AddHsParameters
from rdkit.Chem.rdmolops import AdjustQueryParameters
from rdkit.Chem.rdmolops import AdjustQueryWhichFlags
from rdkit.Chem.rdmolops import AromaticityModel
from rdkit.Chem.rdmolops import BondWedgingParameters
from rdkit.Chem.rdmolops import BoolVector
from rdkit.Chem.rdmolops import MolzipLabel
from rdkit.Chem.rdmolops import MolzipParams
from rdkit.Chem.rdmolops import RemoveHsParameters
from rdkit.Chem.rdmolops import SanitizeFlags
from rdkit.Chem.rdmolops import StereoBondThresholds
from rdkit.Chem.rdmolops import StereoGroupAbsOptions
from rdkit.Chem.rdmolops import SubsetInfo
from rdkit.Chem.rdmolops import SubsetMethod
from rdkit.Chem.rdmolops import SubsetOptions
from rdkit.Chem.rdmolops import UIntUIntMap
from rdkit.Chem.rdmolops import map_indexing_suite_UIntUIntMap_entry
from rdkit import DataStructs
from rdkit.Geometry import rdGeometry
from rdkit import RDConfig
from rdkit import rdBase
from .inchi import *
from .rdCIPLabeler import *
from .rdCoordGen import *
from .rdMolInterchange import *
from .rdchem import *
from .rdinchi import *
from .rdmolfiles import *
from .rdmolops import *
__all__: list[str] = ['ADJUST_IGNOREALL', 'ADJUST_IGNORECHAINS', 'ADJUST_IGNOREDUMMIES', 'ADJUST_IGNOREMAPPED', 'ADJUST_IGNORENONDUMMIES', 'ADJUST_IGNORENONE', 'ADJUST_IGNORERINGS', 'ALLOW_CHARGE_SEPARATION', 'ALLOW_INCOMPLETE_OCTETS', 'AROMATICITY_CUSTOM', 'AROMATICITY_DEFAULT', 'AROMATICITY_MDL', 'AROMATICITY_MMFF94', 'AROMATICITY_RDKIT', 'AROMATICITY_SIMPLE', 'AddHsParameters', 'AdjustQueryParameters', 'AdjustQueryWhichFlags', 'AllProps', 'AromaticityModel', 'Atom', 'AtomCoordsMatcher', 'AtomKekulizeException', 'AtomMonomerInfo', 'AtomMonomerType', 'AtomPDBResidueInfo', 'AtomProps', 'AtomSanitizeException', 'AtomValenceException', 'Bond', 'BondDir', 'BondProps', 'BondStereo', 'BondType', 'BondWedgingParameters', 'BoolVector', 'CDXMLFormat', 'CDXMLParserParams', 'CHI_ALLENE', 'CHI_OCTAHEDRAL', 'CHI_OTHER', 'CHI_SQUAREPLANAR', 'CHI_TETRAHEDRAL', 'CHI_TETRAHEDRAL_CCW', 'CHI_TETRAHEDRAL_CW', 'CHI_TRIGONALBIPYRAMIDAL', 'CHI_UNSPECIFIED', 'COMPOSITE_AND', 'COMPOSITE_OR', 'COMPOSITE_XOR', 'CXSmilesFields', 'CanonSmiles', 'ChiralType', 'CompositeQueryType', 'ComputedProps', 'Conformer', 'CoordsAsDouble', 'DataStructs', 'EXPLICIT', 'EditableMol', 'FindMolChiralCenters', 'FixedMolSizeMolBundle', 'ForwardSDMolSupplier', 'HybridizationType', 'IMPLICIT', 'INCHI_AVAILABLE', 'InchiReadWriteError', 'InchiToInchiKey', 'JSONParseParameters', 'JSONWriteParameters', 'KEKULE_ALL', 'KekulizeException', 'LayeredFingerprint_substructLayers', 'MaeMolSupplier', 'MaeWriter', 'Mol', 'MolBlockToInchi', 'MolBlockToInchiAndAuxInfo', 'MolBundle', 'MolFromInchi', 'MolFromInchiAndAuxInfo', 'MolFromSCSRParams', 'MolProps', 'MolSanitizeException', 'MolToInchi', 'MolToInchiAndAuxInfo', 'MolToInchiKey', 'MolWriterParams', 'MolzipLabel', 'MolzipParams', 'MultiConfMolFromSDF', 'MultithreadedSDMolSupplier', 'MultithreadedSmilesMolSupplier', 'NoConformers', 'NoProps', 'PDBWriter', 'PNGMetadataParams', 'PeriodicTable', 'PrivateProps', 'PropertyPickleOptions', 'QueryAtom', 'QueryAtomData', 'QueryBond', 'QuickSmartsMatch', 'RDConfig', 'RWMol', 'RemoveHsParameters', 'ResonanceFlags', 'ResonanceMolSupplier', 'ResonanceMolSupplierCallback', 'RestoreBondDirOption', 'RingInfo', 'SANITIZE_ADJUSTHS', 'SANITIZE_ALL', 'SANITIZE_CLEANUP', 'SANITIZE_CLEANUPATROPISOMERS', 'SANITIZE_CLEANUPCHIRALITY', 'SANITIZE_CLEANUP_ORGANOMETALLICS', 'SANITIZE_FINDRADICALS', 'SANITIZE_KEKULIZE', 'SANITIZE_NONE', 'SANITIZE_PROPERTIES', 'SANITIZE_SETAROMATICITY', 'SANITIZE_SETCONJUGATION', 'SANITIZE_SETHYBRIDIZATION', 'SANITIZE_SYMMRINGS', 'SCSRBaseHbondOptions', 'SCSRTemplateNames', 'SDMolSupplier', 'SDWriter', 'STEREO_ABSOLUTE', 'STEREO_AND', 'STEREO_OR', 'SanitizeFlags', 'SmartsParserParams', 'SmilesMolSupplier', 'SmilesParserParams', 'SmilesWriteParams', 'SmilesWriter', 'StereoBondThresholds', 'StereoDescriptor', 'StereoGroup', 'StereoGroupAbsOptions', 'StereoGroupType', 'StereoGroup_vect', 'StereoInfo', 'StereoSpecified', 'StereoType', 'SubsetInfo', 'SubsetMethod', 'SubsetOptions', 'SubstanceGroup', 'SubstanceGroupAttach', 'SubstanceGroupCState', 'SubstanceGroup_VECT', 'SubstructMatchParameters', 'SupplierFromFilename', 'TDTMolSupplier', 'TDTWriter', 'UIntUIntMap', 'UNCONSTRAINED_ANIONS', 'UNCONSTRAINED_CATIONS', 'ValenceType', 'inchi', 'map_indexing_suite_UIntUIntMap_entry', 'rdBase', 'rdCIPLabeler', 'rdCoordGen', 'rdGeometry', 'rdMolInterchange', 'rdchem', 'rdinchi', 'rdmolfiles', 'rdmolops', 'templDir']
class _GetAtomsIterator(_GetRDKitObjIterator[Atom]): ...
class _GetBondsIterator(_GetRDKitObjIterator[Bond]): ...
class _GetRDKitObjIterator(typing.Generic[_RDKitObj]):
    def __init__(self, mol: Mol) -> None: ...
    def __iter__(self) -> _GetRDKitObjIterator[_RDKitObj]: ...
    def __next__(self) -> _RDKitObj: ...
    def __getitem__(self, i: int) -> _RDKitObj: ...
    def __len__(self) -> int: ...
def CanonSmiles(smi, useChiral = 1):
    """
     A convenience function for canonicalizing SMILES
    
      Arguments:
        - smi: the SMILES to canonicalize
        - useChiral: (optional) determines whether or not chiral information is included in the canonicalization and SMILES
    
      Returns:
        the canonical SMILES
    
      
    """
def FindMolChiralCenters(mol, force = True, includeUnassigned = False, includeCIP = True, useLegacyImplementation = None):
    """
     returns information about the chiral centers in a molecule
    
      Arguments:
        - mol: the molecule to work with
        - force: (optional) if True, stereochemistry will be assigned even if it has been already
        - includeUnassigned: (optional) if True, unassigned stereo centers will be included in the output
        - includeCIP: (optional) if True, the CIP code for each chiral center will be included in the output
        - useLegacyImplementation: (optional) if True, the legacy stereochemistry perception code will be used
    
      Returns:
        a list of tuples of the form (atomId, CIPCode)
    
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles('[C@H](Cl)(F)Br')
        >>> Chem.FindMolChiralCenters(mol)
        [(0, 'R')]
        >>> mol = Chem.MolFromSmiles('[C@@H](Cl)(F)Br')
        >>> Chem.FindMolChiralCenters(mol)
        [(0, 'S')]
    
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('CCC'))
        []
    
        By default unassigned stereo centers are not reported:
    
        >>> mol = Chem.MolFromSmiles('C[C@H](F)C(F)(Cl)Br')
        >>> Chem.FindMolChiralCenters(mol,force=True)
        [(1, 'S')]
    
        but this can be changed:
    
        >>> Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True)
        [(1, 'S'), (3, '?')]
    
        The handling of unassigned stereocenters for dependent stereochemistry is not correct 
        using the legacy implementation:
    
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1CC(C)C(C)C(C)C1'),includeUnassigned=True, useLegacyImplementation=True)
        [(2, '?'), (6, '?')]
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)C(C)[C@H](C)C1'),includeUnassigned=True, useLegacyImplementation=True)
        [(2, 'S'), (4, '?'), (6, 'R')]
    
        But works with the new implementation:
    
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1CC(C)C(C)C(C)C1'),includeUnassigned=True, useLegacyImplementation=False)
        [(2, '?'), (4, '?'), (6, '?')]
    
        Note that the new implementation also gets the correct descriptors for para-stereochemistry:
    
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)[C@H](C)[C@H](C)C1'),useLegacyImplementation=False)
        [(2, 'S'), (4, 's'), (6, 'R')]
    
        With the new implementation, if you don't care about the CIP labels of stereocenters, you can save
        some time by disabling those:
    
        >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)[C@H](C)[C@H](C)C1'), includeCIP=False, useLegacyImplementation=False)
        [(2, 'Tet_CCW'), (4, 'Tet_CCW'), (6, 'Tet_CCW')]
    
      
    """
def MultiConfMolFromSDF(sdf, sanitize = True, removeHs = True, strictParsing = True):
    """
     A convenience function for creating a multi conformer molecule from an SD file.
    
      This assumes that all structures in the same file are actually the same molecule.
      
      Arguments:
        - sdf: the name of the file to read from
        - sanitize: (optional) Sanitize the molecule after constructing it
        - removeHs: (optional) Remove Hs after constructing the molecule
        - strictParsing: (optional) If set to false, the parser is more lax about the correctness of the contents
    
      Returns:
        A molecule with all conformers found in the SD file.
    
      
    """
def QuickSmartsMatch(smi, sma, unique = True, display = False):
    """
     A convenience function for quickly matching a SMARTS against a SMILES
    
      Arguments:
        - smi: the SMILES to match
        - sma: the SMARTS to match
        - unique: (optional) determines whether or not only unique matches are returned
        - display: (optional) IGNORED
    
      Returns:
        a list of list of the indices of the atoms in the molecule that match the SMARTS  
      
      
    """
def SupplierFromFilename(fileN, delim = '', **kwargs):
    """
     A convenience function for creating a molecule supplier from a filename 
      
      Arguments:
        - fileN: the name of the file to read from
        - delim: (optional) the delimiter to use for reading the file (only for csv and txt files)
        - kwargs: additional keyword arguments to be passed to the supplier constructor
    
      Returns:
        a molecule supplier
    
      
    """
ADJUST_IGNOREALL: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL
ADJUST_IGNORECHAINS: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
ADJUST_IGNOREDUMMIES: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
ADJUST_IGNOREMAPPED: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED
ADJUST_IGNORENONDUMMIES: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES
ADJUST_IGNORENONE: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE
ADJUST_IGNORERINGS: rdmolops.AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS
ALLOW_CHARGE_SEPARATION: rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
ALLOW_INCOMPLETE_OCTETS: rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
AROMATICITY_CUSTOM: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM
AROMATICITY_DEFAULT: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT
AROMATICITY_MDL: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL
AROMATICITY_MMFF94: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94
AROMATICITY_RDKIT: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT
AROMATICITY_SIMPLE: rdmolops.AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE
AllProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
AtomProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
BondProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
CHI_ALLENE: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
CHI_OCTAHEDRAL: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
CHI_OTHER: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
CHI_SQUAREPLANAR: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
CHI_TETRAHEDRAL: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
CHI_TETRAHEDRAL_CCW: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
CHI_TETRAHEDRAL_CW: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
CHI_TRIGONALBIPYRAMIDAL: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
CHI_UNSPECIFIED: rdchem.ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
COMPOSITE_AND: rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
COMPOSITE_OR: rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
COMPOSITE_XOR: rdchem.CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
ComputedProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
CoordsAsDouble: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
EXPLICIT: rdchem.ValenceType  # value = rdkit.Chem.rdchem.ValenceType.EXPLICIT
IMPLICIT: rdchem.ValenceType  # value = rdkit.Chem.rdchem.ValenceType.IMPLICIT
INCHI_AVAILABLE: bool = True
KEKULE_ALL: rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
LayeredFingerprint_substructLayers: int = 7
MolProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
NoConformers: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData: rdchem.PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
SANITIZE_ADJUSTHS: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
SANITIZE_ALL: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL
SANITIZE_CLEANUP: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
SANITIZE_CLEANUPATROPISOMERS: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS
SANITIZE_CLEANUPCHIRALITY: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
SANITIZE_CLEANUP_ORGANOMETALLICS: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
SANITIZE_FINDRADICALS: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
SANITIZE_KEKULIZE: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
SANITIZE_NONE: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
SANITIZE_PROPERTIES: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
SANITIZE_SETAROMATICITY: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
SANITIZE_SETCONJUGATION: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
SANITIZE_SETHYBRIDIZATION: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
SANITIZE_SYMMRINGS: rdmolops.SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
STEREO_ABSOLUTE: rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND: rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR: rdchem.StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
UNCONSTRAINED_ANIONS: rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS: rdchem.ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
templDir: str = '/Users/runner/work/rdkit-pypi/rdkit-pypi/build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/share/RDKit/Data/'

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdmolops import AddHs as AddHs
from rdkit.Chem.rdmolfiles import AddMetadataToPNGFile as AddMetadataToPNGFile
from rdkit.Chem.rdmolfiles import AddMetadataToPNGString as AddMetadataToPNGString
from rdkit.Chem.rdchem import AddMolSubstanceGroup as AddMolSubstanceGroup
from rdkit.Chem.rdmolops import AddRecursiveQuery as AddRecursiveQuery
from rdkit.Chem.rdmolops import AddStereoAnnotations as AddStereoAnnotations
from rdkit.Chem.rdmolops import AddWavyBondsForStereoAny as AddWavyBondsForStereoAny
from rdkit.Chem.rdmolops import AdjustQueryProperties as AdjustQueryProperties
from rdkit.Chem.rdmolops import AdjustQueryPropertiesWithGenericGroups as AdjustQueryPropertiesWithGenericGroups
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
from rdkit.Chem.rdmolfiles import BondFromSmarts as BondFromSmarts
from rdkit.Chem.rdmolfiles import BondFromSmiles as BondFromSmiles
from rdkit.Chem.rdmolfiles import CanonicalRankAtoms as CanonicalRankAtoms
from rdkit.Chem.rdmolfiles import CanonicalRankAtomsInFragment as CanonicalRankAtomsInFragment
from rdkit.Chem.rdmolfiles import CanonicalizeEnhancedStereo as CanonicalizeEnhancedStereo
from rdkit.Chem.rdmolops import CanonicalizeStereoGroups as CanonicalizeStereoGroups
from rdkit.Chem.rdmolops import Cleanup as Cleanup
from rdkit.Chem.rdmolops import CleanupAtropisomers as CleanupAtropisomers
from rdkit.Chem.rdmolops import CleanupChirality as CleanupChirality
from rdkit.Chem.rdmolops import CleanupOrganometallics as CleanupOrganometallics
from rdkit.Chem.rdmolops import CleanupStereoGroups as CleanupStereoGroups
from rdkit.Chem.rdchem import ClearMolSubstanceGroups as ClearMolSubstanceGroups
from rdkit.Chem.rdmolops import CollapseAttachmentPoints as CollapseAttachmentPoints
from rdkit.Chem.rdmolops import CombineMols as CombineMols
from rdkit.Chem.rdmolops import ComputeAtomCIPRanks as ComputeAtomCIPRanks
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
from rdkit.Chem.rdchem import CreateMolDataSubstanceGroup as CreateMolDataSubstanceGroup
from rdkit.Chem.rdchem import CreateMolSubstanceGroup as CreateMolSubstanceGroup
from rdkit.Chem.rdchem import CreateStereoGroup as CreateStereoGroup
from rdkit.Chem.rdmolops import DativeBondsToHaptic as DativeBondsToHaptic
from rdkit.Chem.rdmolops import DeleteSubstructs as DeleteSubstructs
from rdkit.Chem.rdmolops import DetectBondStereoChemistry as DetectBondStereoChemistry
from rdkit.Chem.rdmolops import DetectBondStereochemistry as DetectBondStereochemistry
from rdkit.Chem.rdmolops import DetectChemistryProblems as DetectChemistryProblems
from rdkit.Chem.rdmolops import ExpandAttachmentPoints as ExpandAttachmentPoints
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
from rdkit.Chem.rdchem import ForwardStereoGroupIds as ForwardStereoGroupIds
from rdkit.Chem.rdmolops import FragmentOnBRICSBonds as FragmentOnBRICSBonds
from rdkit.Chem.rdmolops import FragmentOnBonds as FragmentOnBonds
from rdkit.Chem.rdmolops import FragmentOnSomeBonds as FragmentOnSomeBonds
from rdkit.Chem.rdmolops import Get3DDistanceMatrix as Get3DDistanceMatrix
from rdkit.Chem.rdmolops import GetAdjacencyMatrix as GetAdjacencyMatrix
from rdkit.Chem.rdmolops import GetAllowNontetrahedralChirality as GetAllowNontetrahedralChirality
from rdkit.Chem.rdchem import GetAtomAlias as GetAtomAlias
from rdkit.Chem.rdchem import GetAtomRLabel as GetAtomRLabel
from rdkit.Chem.rdchem import GetAtomValue as GetAtomValue
from rdkit.Chem.rdchem import GetDefaultPickleProperties as GetDefaultPickleProperties
from rdkit.Chem.rdmolops import GetDistanceMatrix as GetDistanceMatrix
from rdkit.Chem.rdmolops import GetFormalCharge as GetFormalCharge
from rdkit.Chem.rdinchi import GetInchiVersion as GetInchiVersion
from rdkit.Chem.rdmolops import GetMolFrags as GetMolFrags
from rdkit.Chem.rdchem import GetMolSubstanceGroupWithIdx as GetMolSubstanceGroupWithIdx
from rdkit.Chem.rdchem import GetMolSubstanceGroups as GetMolSubstanceGroups
from rdkit.Chem.rdmolops import GetMostSubstitutedCoreMatch as GetMostSubstitutedCoreMatch
from rdkit.Chem.rdchem import GetNumPiElectrons as GetNumPiElectrons
from rdkit.Chem.rdchem import GetPeriodicTable as GetPeriodicTable
from rdkit.Chem.rdmolops import GetSSSR as GetSSSR
from rdkit.Chem.rdmolops import GetShortestPath as GetShortestPath
from rdkit.Chem.rdchem import GetSupplementalSmilesLabel as GetSupplementalSmilesLabel
from rdkit.Chem.rdmolops import GetSymmSSSR as GetSymmSSSR
from rdkit.Chem.rdmolops import GetUseLegacyStereoPerception as GetUseLegacyStereoPerception
from rdkit.Chem.rdmolops import HapticBondsToDative as HapticBondsToDative
from rdkit.Chem.rdmolfiles import HasChemDrawCDXSupport as HasChemDrawCDXSupport
from rdkit.Chem.rdmolops import HasQueryHs as HasQueryHs
from rdkit.Chem.rdMolInterchange import JSONToMols as JSONToMols
from rdkit.Chem.rdmolops import Kekulize as Kekulize
from rdkit.Chem.rdmolops import KekulizeIfPossible as KekulizeIfPossible
from rdkit.Chem.rdmolops import LayeredFingerprint as LayeredFingerprint
from rdkit.Chem.rdmolops import MergeQueryHs as MergeQueryHs
from rdkit.Chem.rdmolfiles import MetadataFromPNGFile as MetadataFromPNGFile
from rdkit.Chem.rdmolfiles import MetadataFromPNGString as MetadataFromPNGString
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
from rdkit.Chem.rdmolfiles import MolFromRDKitSVG as MolFromRDKitSVG
from rdkit.Chem.rdmolfiles import MolFromSCSRBlock as MolFromSCSRBlock
from rdkit.Chem.rdmolfiles import MolFromSCSRFile as MolFromSCSRFile
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
from rdkit.Chem.rdmolops import MurckoDecompose as MurckoDecompose
from rdkit.Chem.rdmolops import NeedsHs as NeedsHs
from rdkit.Chem.rdmolops import ParseMolQueryDefFile as ParseMolQueryDefFile
from rdkit.Chem.rdmolops import PathToSubmol as PathToSubmol
from rdkit.Chem.rdmolops import PatternFingerprint as PatternFingerprint
from rdkit.Chem.rdmolops import RDKFingerprint as RDKFingerprint
from rdkit.Chem.rdmolops import ReapplyMolBlockWedging as ReapplyMolBlockWedging
from rdkit.Chem.rdmolops import RemoveAllHs as RemoveAllHs
from rdkit.Chem.rdmolops import RemoveHs as RemoveHs
from rdkit.Chem.rdmolops import RemoveNonExplicit3DChirality as RemoveNonExplicit3DChirality
from rdkit.Chem.rdmolops import RemoveStereochemistry as RemoveStereochemistry
from rdkit.Chem.rdmolops import RenumberAtoms as RenumberAtoms
from rdkit.Chem.rdmolops import ReplaceCore as ReplaceCore
from rdkit.Chem.rdmolops import ReplaceSidechains as ReplaceSidechains
from rdkit.Chem.rdmolops import ReplaceSubstructs as ReplaceSubstructs
from rdkit.Chem.rdmolops import SanitizeMol as SanitizeMol
from rdkit.Chem.rdmolops import SetAllowNontetrahedralChirality as SetAllowNontetrahedralChirality
from rdkit.Chem.rdmolops import SetAromaticity as SetAromaticity
from rdkit.Chem.rdchem import SetAtomAlias as SetAtomAlias
from rdkit.Chem.rdchem import SetAtomRLabel as SetAtomRLabel
from rdkit.Chem.rdchem import SetAtomValue as SetAtomValue
from rdkit.Chem.rdmolops import SetBondStereoFromDirections as SetBondStereoFromDirections
from rdkit.Chem.rdmolops import SetConjugation as SetConjugation
from rdkit.Chem.rdchem import SetDefaultPickleProperties as SetDefaultPickleProperties
from rdkit.Chem.rdmolops import SetDoubleBondNeighborDirections as SetDoubleBondNeighborDirections
from rdkit.Chem.rdmolops import SetGenericQueriesFromProperties as SetGenericQueriesFromProperties
from rdkit.Chem.rdmolops import SetHybridization as SetHybridization
from rdkit.Chem.rdchem import SetSupplementalSmilesLabel as SetSupplementalSmilesLabel
from rdkit.Chem.rdmolops import SetTerminalAtomCoords as SetTerminalAtomCoords
from rdkit.Chem.rdmolops import SetUseLegacyStereoPerception as SetUseLegacyStereoPerception
from rdkit.Chem.rdmolops import SimplifyEnhancedStereo as SimplifyEnhancedStereo
from rdkit.Chem.rdmolfiles import SmilesMolSupplierFromText as SmilesMolSupplierFromText
from rdkit.Chem.rdmolops import SortMatchesByDegreeOfCoreSubstitution as SortMatchesByDegreeOfCoreSubstitution
from rdkit.Chem.rdmolops import SplitMolByPDBChainId as SplitMolByPDBChainId
from rdkit.Chem.rdmolops import SplitMolByPDBResidues as SplitMolByPDBResidues
from rdkit.Chem.rdmolops import TranslateChiralFlagToStereoGroups as TranslateChiralFlagToStereoGroups
from rdkit.Chem.rdmolops import UnfoldedRDKFingerprintCountBased as UnfoldedRDKFingerprintCountBased
from rdkit.Chem.rdmolops import WedgeBond as WedgeBond
from rdkit.Chem.rdmolops import WedgeMolBonds as WedgeMolBonds
from . import inchi as inchi
from rdkit.Chem.rdmolops import molzip as molzip
from rdkit.Chem.rdmolops import molzipFragments as molzipFragments
from . import rdCIPLabeler as rdCIPLabeler
from . import rdCoordGen as rdCoordGen
from . import rdMolInterchange as rdMolInterchange
from . import rdchem as rdchem
from . import rdinchi as rdinchi
from . import rdmolfiles as rdmolfiles
from . import rdmolops as rdmolops
from rdkit.Chem.rdchem import tossit as tossit
