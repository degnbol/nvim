# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing the core chemistry functionality of the RDKit
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['ALLOW_CHARGE_SEPARATION', 'ALLOW_INCOMPLETE_OCTETS', 'AddMolSubstanceGroup', 'AllProps', 'Atom', 'AtomCoordsMatcher', 'AtomKekulizeException', 'AtomMonomerInfo', 'AtomMonomerType', 'AtomPDBResidueInfo', 'AtomProps', 'AtomSanitizeException', 'AtomValenceException', 'Bond', 'BondDir', 'BondProps', 'BondStereo', 'BondType', 'CHI_ALLENE', 'CHI_OCTAHEDRAL', 'CHI_OTHER', 'CHI_SQUAREPLANAR', 'CHI_TETRAHEDRAL', 'CHI_TETRAHEDRAL_CCW', 'CHI_TETRAHEDRAL_CW', 'CHI_TRIGONALBIPYRAMIDAL', 'CHI_UNSPECIFIED', 'COMPOSITE_AND', 'COMPOSITE_OR', 'COMPOSITE_XOR', 'ChiralType', 'ClearMolSubstanceGroups', 'CompositeQueryType', 'ComputedProps', 'Conformer', 'CoordsAsDouble', 'CreateMolDataSubstanceGroup', 'CreateMolSubstanceGroup', 'CreateStereoGroup', 'EXPLICIT', 'EditableMol', 'FixedMolSizeMolBundle', 'ForwardStereoGroupIds', 'GetAtomAlias', 'GetAtomRLabel', 'GetAtomValue', 'GetDefaultPickleProperties', 'GetMolSubstanceGroupWithIdx', 'GetMolSubstanceGroups', 'GetNumPiElectrons', 'GetPeriodicTable', 'GetSupplementalSmilesLabel', 'HybridizationType', 'IMPLICIT', 'KEKULE_ALL', 'KekulizeException', 'Mol', 'MolBundle', 'MolBundleCanSerialize', 'MolProps', 'MolSanitizeException', 'NoConformers', 'NoProps', 'PeriodicTable', 'PrivateProps', 'PropertyPickleOptions', 'QueryAtom', 'QueryAtomData', 'QueryBond', 'RWMol', 'ResonanceFlags', 'ResonanceMolSupplier', 'ResonanceMolSupplierCallback', 'RingInfo', 'STEREO_ABSOLUTE', 'STEREO_AND', 'STEREO_OR', 'SetAtomAlias', 'SetAtomRLabel', 'SetAtomValue', 'SetDefaultPickleProperties', 'SetSupplementalSmilesLabel', 'StereoDescriptor', 'StereoGroup', 'StereoGroupType', 'StereoGroup_vect', 'StereoInfo', 'StereoSpecified', 'StereoType', 'SubstanceGroup', 'SubstanceGroupAttach', 'SubstanceGroupCState', 'SubstanceGroup_VECT', 'SubstructMatchParameters', 'UNCONSTRAINED_ANIONS', 'UNCONSTRAINED_CATIONS', 'ValenceType', 'tossit']
class Atom(Boost.Python.instance):
    """
    The class to store Atoms.
    Note that, though it is possible to create one, having an Atom on its own
    (i.e not associated with a molecule) is not particularly useful.
    """
    NOATOM: typing.ClassVar[int] = 4294967295
    __instance_size__: typing.ClassVar[int] = 104
    @staticmethod
    def HasValenceViolation(arg1: Atom) -> bool:
        """
            Returns whether the atom has a valence violation or not.
            
            
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def ClearProp(self, key: str) -> None:
        """
            Removes a particular property from an Atom (does nothing if not already set).
            
              ARGUMENTS:
                - key: the name of the property to be removed.
            
        
        """
    def ClearPropertyCache(self) -> None:
        """
            Clears implicit and explicit valence information.
            
            
        
        """
    def DescribeQuery(self) -> str:
        """
            returns a text description of the query. Primarily intended for debugging purposes.
            
            
        
        """
    def GetAtomMapNum(self) -> int:
        """
            Gets the atoms map number, returns 0 if not set
        
        """
    def GetAtomicNum(self) -> int:
        """
            Returns the atomic number.
        
        """
    def GetBonds(self) -> tuple:
        """
            Returns a read-only sequence of the atom's bonds
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a bool).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str, default: bool) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a bool).
            
                - default: value to return if the property is not present.
            
              RETURNS: a bool, or default if the property is not present.
            
        
        """
    def GetChiralTag(self) -> ChiralType:
        """
        """
    def GetDegree(self) -> int:
        """
            Returns the degree of the atom in the molecule.
            
              The degree of an atom is defined to be its number of
              directly-bonded neighbors.
              The degree is independent of bond orders, but is dependent
                on whether or not Hs are explicit in the graph.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str, default: float) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
                - default: value to return if the property is not present.
            
              RETURNS: a double, or default if the property is not present.
            
        
        """
    def GetExplicitBitVectProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a ExplicitBitVect).
            
              RETURNS: an ExplicitBitVect 
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    def GetExplicitValence(self) -> int:
        """
            DEPRECATED, please use GetValence(Chem.ValenceType.EXPLICIT) instead.
            Returns the explicit valence of the atom.
            
        
        """
    def GetFormalCharge(self) -> int:
        """
        """
    def GetHybridization(self) -> HybridizationType:
        """
            Returns the atom's hybridization.
            
        
        """
    def GetIdx(self) -> int:
        """
            Returns the atom's index (ordering in the molecule)
            
        
        """
    def GetImplicitValence(self) -> int:
        """
            DEPRECATED, please use GetValence(Chem.ValenceType.IMPLICIT) instead.
            Returns the number of implicit Hs on the atom.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
              RETURNS: an int
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
                - default: value to return if the property is not present.
            
              RETURNS: an int, or default if the property is not present.
            
        
        """
    def GetIsAromatic(self) -> bool:
        """
        """
    def GetIsotope(self) -> int:
        """
        """
    def GetMass(self) -> float:
        """
        """
    def GetMonomerInfo(self) -> AtomMonomerInfo:
        """
            Returns the atom's MonomerInfo object, if there is one.
            
            
        
        """
    def GetNeighbors(self) -> tuple:
        """
            Returns a read-only sequence of the atom's neighbors
            
        
        """
    def GetNoImplicit(self) -> bool:
        """
            Returns whether or not the atom is *allowed* to have implicit Hs.
            
        
        """
    def GetNumExplicitHs(self) -> int:
        """
        """
    def GetNumImplicitHs(self) -> int:
        """
            Returns the total number of implicit Hs on the atom.
            
        
        """
    def GetNumRadicalElectrons(self) -> int:
        """
        """
    def GetOwningMol(self) -> Mol:
        """
            Returns the Mol that owns this atom.
            
        
        """
    def GetPDBResidueInfo(self) -> AtomPDBResidueInfo:
        """
            Returns the atom's MonomerInfo object, if there is one.
            
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False, default: typing.Any = ...) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
                - default: value to return if the property is not present.
            
              RETURNS: the property value, or default if the property is not present.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a list of the properties set on the Atom.
            
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary populated with properties.
            When possible, string values will be converted to integers or doubles (trimming if necessary)
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
                - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                                  Defaults to True.
            
              RETURNS: a dictionary
            
        
        """
    def GetQueryType(self) -> str:
        """
        """
    def GetSmarts(self, doKekule: bool = False, allHsExplicit: bool = False, isomericSmiles: bool = True) -> str:
        """
            returns the SMARTS (or SMILES) string for an Atom
            
            
        
        """
    def GetSymbol(self) -> str:
        """
            Returns the atomic symbol (a string)
            
        
        """
    def GetTotalDegree(self) -> int:
        """
            Returns the degree of the atom in the molecule including Hs.
            
              The degree of an atom is defined to be its number of
              directly-bonded neighbors.
              The degree is independent of bond orders.
            
        
        """
    def GetTotalNumHs(self, includeNeighbors: bool = False) -> int:
        """
            Returns the total number of Hs (explicit and implicit) on the atom.
            
              ARGUMENTS:
            
                - includeNeighbors: (optional) toggles inclusion of neighboring H atoms in the sum.
                  Defaults to 0.
            
        
        """
    def GetTotalValence(self) -> int:
        """
            Returns the total valence (explicit + implicit) of the atom.
            
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
              RETURNS: an integer (Python has no unsigned type)
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
                - default: value to return if the property is not present.
            
              RETURNS: an integer, or default if the property is not present.
            
        
        """
    def GetValence(self, which: ValenceType) -> int:
        """
            Returns the valence (Chem.ValenceType.EXPLICIT or Chem.ValenceType.IMPLICIT) of the atom.
            
        
        """
    def HasOwningMol(self) -> bool:
        """
            Returns whether or not this instance belongs to a molecule.
            
        
        """
    def HasProp(self, key: str) -> int:
        """
            Queries a Atom to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            
        
        """
    def HasQuery(self) -> bool:
        """
            Returns whether or not the atom has an associated query
            
            
        
        """
    def InvertChirality(self) -> bool:
        """
        """
    def IsInRing(self) -> bool:
        """
            Returns whether or not the atom is in a ring
            
            
        
        """
    def IsInRingSize(self, size: int) -> bool:
        """
            Returns whether or not the atom is in a ring of a particular size.
            
              ARGUMENTS:
                - size: the ring size to look for
            
        
        """
    def Match(self, what: Atom) -> bool:
        """
            Returns whether or not this atom matches another Atom.
            
              Each Atom (or query Atom) has a query function which is
              used for this type of matching.
            
              ARGUMENTS:
                - other: the other Atom to which to compare
            
        
        """
    def NeedsUpdatePropertyCache(self) -> bool:
        """
            Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
            
            
        
        """
    def SetAtomMapNum(self, mapno: int, strict: bool = False) -> None:
        """
            Sets the atoms map number, a value of 0 clears the atom map
        
        """
    def SetAtomicNum(self, newNum: int) -> None:
        """
            Sets the atomic number, takes an integer value as an argument
        
        """
    def SetBoolProp(self, key: str, val: bool) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a bool).
                - value: the property value (a bool).
            
            
        
        """
    def SetChiralTag(self, what: ChiralType) -> None:
        """
        """
    def SetDoubleProp(self, key: str, val: float) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a double).
                - value: the property value (a double).
            
            
        
        """
    def SetExplicitBitVectProp(self, key: str, val: ExplicitBitVect) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (an ExplicitBitVect).
                - value: the property value (an ExplicitBitVect).
            
            
        
        """
    def SetFormalCharge(self, what: int) -> None:
        """
        """
    def SetHybridization(self, what: HybridizationType) -> None:
        """
            Sets the hybridization of the atom.
              The argument should be a HybridizationType
            
        
        """
    def SetIntProp(self, key: str, val: int) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a int).
                - value: the property value (a int).
            
            
        
        """
    def SetIsAromatic(self, what: bool) -> None:
        """
        """
    def SetIsotope(self, what: int) -> None:
        """
        """
    def SetMonomerInfo(self, info: AtomMonomerInfo) -> None:
        """
            Sets the atom's MonomerInfo object.
            
            
        
        """
    def SetNoImplicit(self, what: bool) -> None:
        """
            Sets a marker on the atom that *disallows* implicit Hs.
              This holds even if the atom would otherwise have implicit Hs added.
            
        
        """
    def SetNumExplicitHs(self, what: int) -> None:
        """
        """
    def SetNumRadicalElectrons(self, num: int) -> None:
        """
        """
    def SetPDBResidueInfo(self, info: AtomMonomerInfo) -> None:
        """
            Sets the atom's MonomerInfo object.
            
            
        
        """
    def SetProp(self, key: str, val: str) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
            
            
        
        """
    def SetUnsignedProp(self, key: str, val: int) -> None:
        """
            Sets an atomic property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned integer).
                - value: the property value (a int >= 0).
            
            
        
        """
    def UpdatePropertyCache(self, strict: bool = True) -> None:
        """
            Regenerates computed properties like implicit valence and ring information.
            
            
        
        """
    def __copy__(self) -> Atom:
        """
            Create a copy of the atom
        
        """
    @typing.overload
    def __init__(self, what: str) -> None:
        ...
    @typing.overload
    def __init__(self, other: Atom) -> None:
        ...
    @typing.overload
    def __init__(self, num: int) -> None:
        ...
class AtomCoordsMatcher(Boost.Python.instance):
    """
    Allows using atom coordinates as part of substructure matching
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, arg1: Atom, arg2: Atom) -> bool:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, refConfId: int = -1, queryConfId: int = -1, tol: float = 0.0001) -> None:
        ...
    @property
    def queryConfId(self) -> int:
        """query conformer ID (default: -1)"""
    @queryConfId.setter
    def queryConfId(self, value: int) -> None: ...
    @property
    def refConfId(self) -> int:
        """reference conformer ID (default: -1)"""
    @refConfId.setter
    def refConfId(self, value: int) -> None: ...
    @property
    def tol2(self) -> float:
        """squared distance tolerance (default: 1e-08)"""
    @tol2.setter
    def tol2(self, value: float) -> None: ...
class AtomKekulizeException(AtomSanitizeException):
    pass
class AtomMonomerInfo(Boost.Python.instance):
    """
    The class to store monomer information attached to Atoms
    """
    __instance_size__: typing.ClassVar[int] = 144
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetChainId(self) -> str:
        """
        """
    def GetMonomerClass(self) -> str:
        """
        """
    def GetMonomerType(self) -> AtomMonomerType:
        """
        """
    def GetName(self) -> str:
        """
        """
    def GetResidueName(self) -> str:
        """
        """
    def GetResidueNumber(self) -> int:
        """
        """
    def SetChainId(self, val: str) -> None:
        """
        """
    def SetMonomerClass(self, val: str) -> None:
        """
        """
    def SetMonomerType(self, typ: AtomMonomerType) -> None:
        """
        """
    def SetName(self, nm: str) -> None:
        """
        """
    def SetResidueName(self, val: str) -> None:
        """
        """
    def SetResidueNumber(self, val: int) -> None:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, type: AtomMonomerType, name: str = '', residueName: str = '', resNum: int = 0, chainId: str = '', monomerClass: str = '') -> None:
        ...
class AtomMonomerType(Boost.Python.enum):
    OTHER: typing.ClassVar[AtomMonomerType]  # value = rdkit.Chem.rdchem.AtomMonomerType.OTHER
    PDBRESIDUE: typing.ClassVar[AtomMonomerType]  # value = rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE
    UNKNOWN: typing.ClassVar[AtomMonomerType]  # value = rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'UNKNOWN': rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN, 'PDBRESIDUE': rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE, 'OTHER': rdkit.Chem.rdchem.AtomMonomerType.OTHER}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.AtomMonomerType.UNKNOWN, 1: rdkit.Chem.rdchem.AtomMonomerType.PDBRESIDUE, 2: rdkit.Chem.rdchem.AtomMonomerType.OTHER}
class AtomPDBResidueInfo(AtomMonomerInfo):
    """
    The class to store PDB residue information attached to Atoms
    """
    __instance_size__: typing.ClassVar[int] = 232
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetAltLoc(self) -> str:
        """
        """
    def GetChainId(self) -> str:
        """
        """
    def GetInsertionCode(self) -> str:
        """
        """
    def GetIsHeteroAtom(self) -> bool:
        """
        """
    def GetMonomerClass(self) -> str:
        """
        """
    def GetOccupancy(self) -> float:
        """
        """
    def GetResidueName(self) -> str:
        """
        """
    def GetResidueNumber(self) -> int:
        """
        """
    def GetSecondaryStructure(self) -> int:
        """
        """
    def GetSegmentNumber(self) -> int:
        """
        """
    def GetSerialNumber(self) -> int:
        """
        """
    def GetTempFactor(self) -> float:
        """
        """
    def SetAltLoc(self, val: str) -> None:
        """
        """
    def SetChainId(self, val: str) -> None:
        """
        """
    def SetInsertionCode(self, val: str) -> None:
        """
        """
    def SetIsHeteroAtom(self, val: bool) -> None:
        """
        """
    def SetMonomerClass(self, val: str) -> None:
        """
        """
    def SetOccupancy(self, val: float) -> None:
        """
        """
    def SetResidueName(self, val: str) -> None:
        """
        """
    def SetResidueNumber(self, val: int) -> None:
        """
        """
    def SetSecondaryStructure(self, val: int) -> None:
        """
        """
    def SetSegmentNumber(self, val: int) -> None:
        """
        """
    def SetSerialNumber(self, val: int) -> None:
        """
        """
    def SetTempFactor(self, val: float) -> None:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atomName: str, serialNumber: int = 1, altLoc: str = '', residueName: str = '', residueNumber: int = 0, chainId: str = '', insertionCode: str = '', occupancy: float = 1.0, tempFactor: float = 0.0, isHeteroAtom: bool = False, secondaryStructure: int = 0, segmentNumber: int = 0, monomerClass: str = '') -> None:
        ...
class AtomSanitizeException(MolSanitizeException):
    pass
class AtomValenceException(AtomSanitizeException):
    pass
class Bond(Boost.Python.instance):
    """
    The class to store Bonds.
    Note: unlike Atoms, is it currently impossible to construct Bonds from
    Python.
    """
    @staticmethod
    def GetSmarts(bond: Bond, allBondsExplicit: bool = False) -> str:
        """
            returns the SMARTS (or SMILES) string for a Bond
        
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
    def ClearProp(self, key: str) -> None:
        """
            Removes a particular property from an Bond (does nothing if not already set).
            
              ARGUMENTS:
                - key: the name of the property to be removed.
            
        
        """
    def DescribeQuery(self) -> str:
        """
            returns a text description of the query. Primarily intended for debugging purposes.
            
            
        
        """
    def GetBeginAtom(self) -> Atom:
        """
            Returns the bond's first atom.
            
        
        """
    def GetBeginAtomIdx(self) -> int:
        """
            Returns the index of the bond's first atom.
            
        
        """
    def GetBondDir(self) -> BondDir:
        """
            Returns the type of the bond as a BondDir
            
        
        """
    def GetBondType(self) -> BondType:
        """
            Returns the type of the bond as a BondType
            
        
        """
    def GetBondTypeAsDouble(self) -> float:
        """
            Returns the type of the bond as a double (i.e. 1.0 for SINGLE, 1.5 for AROMATIC, 2.0 for DOUBLE)
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a boolean).
            
              RETURNS: a boolean
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str, default: bool) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a boolean).
            
                - default: value to return if the property is not present.
            
              RETURNS: a bool, or default if the property is not present.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str, default: float) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a double).
            
                - default: value to return if the property is not present.
            
              RETURNS: a double, or default if the property is not present.
            
        
        """
    def GetEndAtom(self) -> Atom:
        """
            Returns the bond's second atom.
            
        
        """
    def GetEndAtomIdx(self) -> int:
        """
            Returns the index of the bond's first atom.
            
        
        """
    def GetIdx(self) -> int:
        """
            Returns the bond's index (ordering in the molecule)
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
              RETURNS: an int
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an int).
            
                - default: value to return if the property is not present.
            
              RETURNS: an int, or default if the property is not present.
            
        
        """
    def GetIsAromatic(self) -> bool:
        """
        """
    def GetIsConjugated(self) -> bool:
        """
            Returns whether or not the bond is considered to be conjugated.
        
        """
    def GetOtherAtom(self, what: Atom) -> Atom:
        """
            Given one of the bond's atoms, returns the other one.
            
        
        """
    def GetOtherAtomIdx(self, thisIdx: int) -> int:
        """
            Given the index of one of the bond's atoms, returns the
            index of the other.
            
        
        """
    def GetOwningMol(self) -> Mol:
        """
            Returns the Mol that owns this bond.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False, default: typing.Any = ...) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
                - default: value to return if the property is not present.
            
              RETURNS: the property value, or default if the property is not present.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a list of the properties set on the Bond.
            
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary populated with properties.
            When possible, string values will be converted to integers or doubles (trimming if necessary)
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
                - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                                  Defaults to True.
            
              RETURNS: a dictionary
            
        
        """
    def GetStereo(self) -> BondStereo:
        """
            Returns the stereo configuration of the bond as a BondStereo
            
        
        """
    def GetStereoAtoms(self) -> typing.Sequence[int]:
        """
            Returns the indices of the atoms setting this bond's stereochemistry.
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
              RETURNS: an int (Python has no unsigned type)
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (an unsigned integer).
            
                - default: value to return if the property is not present.
            
              RETURNS: an integer, or default if the property is not present.
            
        
        """
    def GetValenceContrib(self, at: Atom) -> float:
        """
            Returns the contribution of the bond to the valence of an Atom.
            
              ARGUMENTS:
            
                - atom: the Atom to consider.
            
        
        """
    def HasOwningMol(self) -> bool:
        """
            Returns whether or not this instance belongs to a molecule.
            
        
        """
    def HasProp(self, key: str) -> int:
        """
            Queries a Bond to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            
        
        """
    def HasQuery(self) -> bool:
        """
            Returns whether or not the bond has an associated query
            
            
        
        """
    def InvertChirality(self) -> bool:
        """
        """
    def IsInRing(self) -> bool:
        """
            Returns whether or not the bond is in a ring of any size.
            
            
        
        """
    def IsInRingSize(self, size: int) -> bool:
        """
            Returns whether or not the bond is in a ring of a particular size.
            
              ARGUMENTS:
                - size: the ring size to look for
            
        
        """
    def Match(self, what: Bond) -> bool:
        """
            Returns whether or not this bond matches another Bond.
            
              Each Bond (or query Bond) has a query function which is
              used for this type of matching.
            
              ARGUMENTS:
                - other: the other Bond to which to compare
            
        
        """
    def SetBondDir(self, what: BondDir) -> None:
        """
            Set the type of the bond as a BondDir
            
        
        """
    def SetBondType(self, bT: BondType) -> None:
        """
            Set the type of the bond as a BondType
            
        
        """
    def SetBoolProp(self, key: str, val: bool) -> None:
        """
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a boolean).
            
            
        
        """
    def SetDoubleProp(self, key: str, val: float) -> None:
        """
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a double).
            
            
        
        """
    def SetIntProp(self, key: str, val: int) -> None:
        """
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (an int).
            
            
        
        """
    def SetIsAromatic(self, what: bool) -> None:
        """
        """
    def SetIsConjugated(self, what: bool) -> None:
        """
        """
    def SetProp(self, key: str, val: str) -> None:
        """
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
            
            
        
        """
    def SetStereo(self, what: BondStereo) -> None:
        """
            Set the stereo configuration of the bond as a BondStereo
            
        
        """
    def SetStereoAtoms(self, bgnIdx: int, endIdx: int) -> None:
        """
            Set the indices of the atoms setting this bond's stereochemistry.
            
        
        """
    def SetUnsignedProp(self, key: str, val: int) -> None:
        """
            Sets a bond property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (an int >= 0).
            
            
        
        """
class BondDir(Boost.Python.enum):
    BEGINDASH: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.BEGINDASH
    BEGINWEDGE: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.BEGINWEDGE
    EITHERDOUBLE: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.EITHERDOUBLE
    ENDDOWNRIGHT: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT
    ENDUPRIGHT: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.ENDUPRIGHT
    NONE: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.NONE
    UNKNOWN: typing.ClassVar[BondDir]  # value = rdkit.Chem.rdchem.BondDir.UNKNOWN
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NONE': rdkit.Chem.rdchem.BondDir.NONE, 'BEGINWEDGE': rdkit.Chem.rdchem.BondDir.BEGINWEDGE, 'BEGINDASH': rdkit.Chem.rdchem.BondDir.BEGINDASH, 'ENDDOWNRIGHT': rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT, 'ENDUPRIGHT': rdkit.Chem.rdchem.BondDir.ENDUPRIGHT, 'EITHERDOUBLE': rdkit.Chem.rdchem.BondDir.EITHERDOUBLE, 'UNKNOWN': rdkit.Chem.rdchem.BondDir.UNKNOWN}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.BondDir.NONE, 1: rdkit.Chem.rdchem.BondDir.BEGINWEDGE, 2: rdkit.Chem.rdchem.BondDir.BEGINDASH, 3: rdkit.Chem.rdchem.BondDir.ENDDOWNRIGHT, 4: rdkit.Chem.rdchem.BondDir.ENDUPRIGHT, 5: rdkit.Chem.rdchem.BondDir.EITHERDOUBLE, 6: rdkit.Chem.rdchem.BondDir.UNKNOWN}
class BondStereo(Boost.Python.enum):
    STEREOANY: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOANY
    STEREOATROPCCW: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOATROPCCW
    STEREOATROPCW: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOATROPCW
    STEREOCIS: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOCIS
    STEREOE: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOE
    STEREONONE: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREONONE
    STEREOTRANS: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOTRANS
    STEREOZ: typing.ClassVar[BondStereo]  # value = rdkit.Chem.rdchem.BondStereo.STEREOZ
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'STEREONONE': rdkit.Chem.rdchem.BondStereo.STEREONONE, 'STEREOANY': rdkit.Chem.rdchem.BondStereo.STEREOANY, 'STEREOZ': rdkit.Chem.rdchem.BondStereo.STEREOZ, 'STEREOE': rdkit.Chem.rdchem.BondStereo.STEREOE, 'STEREOCIS': rdkit.Chem.rdchem.BondStereo.STEREOCIS, 'STEREOTRANS': rdkit.Chem.rdchem.BondStereo.STEREOTRANS, 'STEREOATROPCW': rdkit.Chem.rdchem.BondStereo.STEREOATROPCW, 'STEREOATROPCCW': rdkit.Chem.rdchem.BondStereo.STEREOATROPCCW}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.BondStereo.STEREONONE, 1: rdkit.Chem.rdchem.BondStereo.STEREOANY, 2: rdkit.Chem.rdchem.BondStereo.STEREOZ, 3: rdkit.Chem.rdchem.BondStereo.STEREOE, 4: rdkit.Chem.rdchem.BondStereo.STEREOCIS, 5: rdkit.Chem.rdchem.BondStereo.STEREOTRANS, 6: rdkit.Chem.rdchem.BondStereo.STEREOATROPCW, 7: rdkit.Chem.rdchem.BondStereo.STEREOATROPCCW}
class BondType(Boost.Python.enum):
    AROMATIC: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.AROMATIC
    DATIVE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.DATIVE
    DATIVEL: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.DATIVEL
    DATIVEONE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.DATIVEONE
    DATIVER: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.DATIVER
    DOUBLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.DOUBLE
    FIVEANDAHALF: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.FIVEANDAHALF
    FOURANDAHALF: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.FOURANDAHALF
    HEXTUPLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.HEXTUPLE
    HYDROGEN: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.HYDROGEN
    IONIC: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.IONIC
    ONEANDAHALF: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.ONEANDAHALF
    OTHER: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.OTHER
    QUADRUPLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.QUADRUPLE
    QUINTUPLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.QUINTUPLE
    SINGLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.SINGLE
    THREEANDAHALF: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.THREEANDAHALF
    THREECENTER: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.THREECENTER
    TRIPLE: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.TRIPLE
    TWOANDAHALF: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.TWOANDAHALF
    UNSPECIFIED: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.UNSPECIFIED
    ZERO: typing.ClassVar[BondType]  # value = rdkit.Chem.rdchem.BondType.ZERO
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'UNSPECIFIED': rdkit.Chem.rdchem.BondType.UNSPECIFIED, 'SINGLE': rdkit.Chem.rdchem.BondType.SINGLE, 'DOUBLE': rdkit.Chem.rdchem.BondType.DOUBLE, 'TRIPLE': rdkit.Chem.rdchem.BondType.TRIPLE, 'QUADRUPLE': rdkit.Chem.rdchem.BondType.QUADRUPLE, 'QUINTUPLE': rdkit.Chem.rdchem.BondType.QUINTUPLE, 'HEXTUPLE': rdkit.Chem.rdchem.BondType.HEXTUPLE, 'ONEANDAHALF': rdkit.Chem.rdchem.BondType.ONEANDAHALF, 'TWOANDAHALF': rdkit.Chem.rdchem.BondType.TWOANDAHALF, 'THREEANDAHALF': rdkit.Chem.rdchem.BondType.THREEANDAHALF, 'FOURANDAHALF': rdkit.Chem.rdchem.BondType.FOURANDAHALF, 'FIVEANDAHALF': rdkit.Chem.rdchem.BondType.FIVEANDAHALF, 'AROMATIC': rdkit.Chem.rdchem.BondType.AROMATIC, 'IONIC': rdkit.Chem.rdchem.BondType.IONIC, 'HYDROGEN': rdkit.Chem.rdchem.BondType.HYDROGEN, 'THREECENTER': rdkit.Chem.rdchem.BondType.THREECENTER, 'DATIVEONE': rdkit.Chem.rdchem.BondType.DATIVEONE, 'DATIVE': rdkit.Chem.rdchem.BondType.DATIVE, 'DATIVEL': rdkit.Chem.rdchem.BondType.DATIVEL, 'DATIVER': rdkit.Chem.rdchem.BondType.DATIVER, 'OTHER': rdkit.Chem.rdchem.BondType.OTHER, 'ZERO': rdkit.Chem.rdchem.BondType.ZERO}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.BondType.UNSPECIFIED, 1: rdkit.Chem.rdchem.BondType.SINGLE, 2: rdkit.Chem.rdchem.BondType.DOUBLE, 3: rdkit.Chem.rdchem.BondType.TRIPLE, 4: rdkit.Chem.rdchem.BondType.QUADRUPLE, 5: rdkit.Chem.rdchem.BondType.QUINTUPLE, 6: rdkit.Chem.rdchem.BondType.HEXTUPLE, 7: rdkit.Chem.rdchem.BondType.ONEANDAHALF, 8: rdkit.Chem.rdchem.BondType.TWOANDAHALF, 9: rdkit.Chem.rdchem.BondType.THREEANDAHALF, 10: rdkit.Chem.rdchem.BondType.FOURANDAHALF, 11: rdkit.Chem.rdchem.BondType.FIVEANDAHALF, 12: rdkit.Chem.rdchem.BondType.AROMATIC, 13: rdkit.Chem.rdchem.BondType.IONIC, 14: rdkit.Chem.rdchem.BondType.HYDROGEN, 15: rdkit.Chem.rdchem.BondType.THREECENTER, 16: rdkit.Chem.rdchem.BondType.DATIVEONE, 17: rdkit.Chem.rdchem.BondType.DATIVE, 18: rdkit.Chem.rdchem.BondType.DATIVEL, 19: rdkit.Chem.rdchem.BondType.DATIVER, 20: rdkit.Chem.rdchem.BondType.OTHER, 21: rdkit.Chem.rdchem.BondType.ZERO}
class ChiralType(Boost.Python.enum):
    CHI_ALLENE: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
    CHI_OCTAHEDRAL: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
    CHI_OTHER: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
    CHI_SQUAREPLANAR: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
    CHI_TETRAHEDRAL: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
    CHI_TETRAHEDRAL_CCW: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
    CHI_TETRAHEDRAL_CW: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
    CHI_TRIGONALBIPYRAMIDAL: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
    CHI_UNSPECIFIED: typing.ClassVar[ChiralType]  # value = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'CHI_UNSPECIFIED': rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED, 'CHI_TETRAHEDRAL_CW': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, 'CHI_TETRAHEDRAL_CCW': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 'CHI_OTHER': rdkit.Chem.rdchem.ChiralType.CHI_OTHER, 'CHI_TETRAHEDRAL': rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL, 'CHI_ALLENE': rdkit.Chem.rdchem.ChiralType.CHI_ALLENE, 'CHI_SQUAREPLANAR': rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR, 'CHI_TRIGONALBIPYRAMIDAL': rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL, 'CHI_OCTAHEDRAL': rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED, 1: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW, 2: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW, 3: rdkit.Chem.rdchem.ChiralType.CHI_OTHER, 4: rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL, 5: rdkit.Chem.rdchem.ChiralType.CHI_ALLENE, 6: rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR, 7: rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL, 8: rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL}
class CompositeQueryType(Boost.Python.enum):
    COMPOSITE_AND: typing.ClassVar[CompositeQueryType]  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
    COMPOSITE_OR: typing.ClassVar[CompositeQueryType]  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
    COMPOSITE_XOR: typing.ClassVar[CompositeQueryType]  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'COMPOSITE_AND': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, 'COMPOSITE_OR': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR, 'COMPOSITE_XOR': rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND, 1: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR, 2: rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR}
class Conformer(Boost.Python.instance):
    """
    The class to store 2D or 3D conformation of a molecule
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def ClearComputedProps(self) -> None:
        """
            Removes all computed properties from the conformer.
            
            
        
        """
    def ClearProp(self, key: str) -> None:
        """
            Removes a property from the conformer.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            
        
        """
    def GetAtomPosition(self, aid: int) -> Point3D:
        """
            Get the posistion of an atom
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str) -> typing.Any:
        """
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str, default: bool) -> typing.Any:
        """
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: a bool, or default if the property is not present.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str) -> typing.Any:
        """
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str, default: float) -> typing.Any:
        """
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: a double, or default if the property is not present.
            
        
        """
    def GetId(self) -> int:
        """
            Get the ID of the conformer
        
        """
    @typing.overload
    def GetIntProp(self, key: str) -> typing.Any:
        """
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: an integer, or default if the property is not present.
            
        
        """
    def GetNumAtoms(self) -> int:
        """
            Get the number of atoms in the conformer
            
        
        """
    def GetOwningMol(self) -> Mol:
        """
            Get the owning molecule
            
        
        """
    def GetPositions(self) -> typing.Any:
        """
            Get positions of all the atoms
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False, default: typing.Any = ...) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
                - default: value to return if the property is not present.
            
              RETURNS: the property value, or default if the property is not present.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a tuple with all property names for this conformer.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary populated with properties.
            When possible, string values will be converted to integers or doubles (trimming if necessary)
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
                - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                                  Defaults to True.
            
              RETURNS: a dictionary
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str) -> typing.Any:
        """
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: an unsigned integer, or default if the property is not present.
            
        
        """
    def HasOwningMol(self) -> bool:
        """
            Returns whether or not this instance belongs to a molecule.
            
        
        """
    def HasProp(self, key: str) -> int:
        """
            Queries a conformer to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            
        
        """
    def Is3D(self) -> bool:
        """
            returns the 3D flag of the conformer
            
        
        """
    def Set3D(self, v: bool) -> None:
        """
            Set the 3D flag of the conformer
            
        
        """
    @typing.overload
    def SetAtomPosition(self, aid: int, loc: typing.Any) -> None:
        """
            Set the position of the specified atom
            
        
        """
    @typing.overload
    def SetAtomPosition(self, atomId: int, position: Point3D) -> None:
        """
            Set the position of the specified atom
            
        
        """
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            
        
        """
    def SetId(self, id: int) -> None:
        """
            Set the ID of the conformer
            
        
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetPositions(self, positions: typing.Any) -> None:
        """
            Set positions of all the atoms given a 2D or 3D numpy array of type double
            
        
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, numAtoms: int) -> None:
        ...
    @typing.overload
    def __init__(self, other: Conformer) -> None:
        ...
class EditableMol(Boost.Python.instance):
    """
    an editable molecule class
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddAtom(self, atom: Atom) -> int:
        """
            add an atom, returns the index of the newly added atom
        
        """
    def AddBond(self, beginAtomIdx: int, endAtomIdx: int, order: BondType = ...) -> int:
        """
            add a bond, returns the total number of bonds
        
        """
    def BeginBatchEdit(self) -> None:
        """
            starts batch editing
        
        """
    def CommitBatchEdit(self) -> None:
        """
            finishes batch editing and makes the actual edits
        
        """
    def GetMol(self) -> Mol:
        """
            Returns a Mol (a normal molecule)
        
        """
    def RemoveAtom(self, idx: int) -> None:
        """
            Remove the specified atom from the molecule
        
        """
    def RemoveBond(self, idx1: int, idx2: int) -> None:
        """
            Remove the specified bond from the molecule
        
        """
    def ReplaceAtom(self, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None:
        """
            replaces the specified atom with the provided one
            If updateLabel is True, the new atom becomes the active atom
            If preserveProps is True preserve keep the existing props unless explicit set on the new atom
        
        """
    def ReplaceBond(self, index: int, newBond: Bond, preserveProps: bool = False) -> None:
        """
            replaces the specified bond with the provided one.
            If preserveProps is True preserve keep the existing props unless explicit set on the new bond
        
        """
    def RollbackBatchEdit(self) -> None:
        """
            cancels batch editing
        
        """
    def __init__(self, m: Mol) -> None:
        ...
class FixedMolSizeMolBundle(MolBundle):
    """
    A class for storing groups of related molecules.
        Here related means that the molecules have to have the same number of atoms.
    
    """
    __instance_size__: typing.ClassVar[int] = 88
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class HybridizationType(Boost.Python.enum):
    OTHER: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.OTHER
    S: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.S
    SP: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP
    SP2: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP2
    SP2D: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP2D
    SP3: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP3
    SP3D: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP3D
    SP3D2: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.SP3D2
    UNSPECIFIED: typing.ClassVar[HybridizationType]  # value = rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'UNSPECIFIED': rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED, 'S': rdkit.Chem.rdchem.HybridizationType.S, 'SP': rdkit.Chem.rdchem.HybridizationType.SP, 'SP2': rdkit.Chem.rdchem.HybridizationType.SP2, 'SP3': rdkit.Chem.rdchem.HybridizationType.SP3, 'SP2D': rdkit.Chem.rdchem.HybridizationType.SP2D, 'SP3D': rdkit.Chem.rdchem.HybridizationType.SP3D, 'SP3D2': rdkit.Chem.rdchem.HybridizationType.SP3D2, 'OTHER': rdkit.Chem.rdchem.HybridizationType.OTHER}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.HybridizationType.UNSPECIFIED, 1: rdkit.Chem.rdchem.HybridizationType.S, 2: rdkit.Chem.rdchem.HybridizationType.SP, 3: rdkit.Chem.rdchem.HybridizationType.SP2, 4: rdkit.Chem.rdchem.HybridizationType.SP3, 5: rdkit.Chem.rdchem.HybridizationType.SP2D, 6: rdkit.Chem.rdchem.HybridizationType.SP3D, 7: rdkit.Chem.rdchem.HybridizationType.SP3D2, 8: rdkit.Chem.rdchem.HybridizationType.OTHER}
class KekulizeException(MolSanitizeException):
    pass
class Mol(Boost.Python.instance):
    """
    The Molecule class.
    
      In addition to the expected Atoms and Bonds, molecules contain:
        - a collection of Atom and Bond bookmarks indexed with integers
            that can be used to flag and retrieve particular Atoms or Bonds
            using the {get|set}{Atom|Bond}Bookmark() methods.
    
        - a set of string-valued properties. These can have arbitrary string
            labels and can be set and retrieved using the {set|get}Prop() methods
            Molecular properties can be tagged as being *computed*, in which case
              they will be automatically cleared under certain circumstances (when the
              molecule itself is modified, for example).
            Molecules also have the concept of *private* properties, which are tagged
              by beginning the property name with an underscore (_).
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 40
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddConformer(self, conf: Conformer, assignId: bool = False) -> int:
        """
            Add a conformer to the molecule and return the conformer ID
        
        """
    def ClearComputedProps(self, includeRings: bool = True) -> None:
        """
            Removes all computed properties from the molecule.
            
            
        
        """
    def ClearProp(self, key: str) -> None:
        """
            Removes a property from the molecule.
            
              ARGUMENTS:
                - key: the name of the property to clear (a string).
            
        
        """
    def ClearPropertyCache(self) -> None:
        """
            Clears implicit and explicit valence information from all atoms.
            
            
        
        """
    def Debug(self, useStdout: bool = True) -> None:
        """
            Prints debugging information about the molecule.
            
        
        """
    def GetAromaticAtoms(self) -> typing.Sequence[rdkit.Chem.QueryAtom]:
        """
            Returns a read-only sequence containing all of the molecule's aromatic Atoms.
            
        
        """
    def GetAtomWithIdx(self, idx: int) -> Atom:
        """
            Returns a particular Atom.
            
              ARGUMENTS:
                - idx: which Atom to return
            
              NOTE: atom indices start at 0
            
        
        """
    def GetAtoms(self) -> rdkit.Chem._GetAtomsIterator:
        """
        returns an iterator over the atoms in the molecule
        """
    def GetAtomsMatchingQuery(self, qa: QueryAtom) -> typing.Sequence[rdkit.Chem.QueryAtom]:
        """
            Returns a read-only sequence containing all of the atoms in a molecule that match the query atom. Atom query options are defined in the rdkit.Chem.rdqueries module.
            
        
        """
    def GetBondBetweenAtoms(self, idx1: int, idx2: int) -> Bond:
        """
            Returns the bond between two atoms, if there is one.
            
              ARGUMENTS:
                - idx1,idx2: the Atom indices
            
              Returns:
                The Bond between the two atoms, if such a bond exists.
                If there is no Bond between the atoms, None is returned instead.
            
              NOTE: bond indices start at 0
            
        
        """
    def GetBondWithIdx(self, idx: int) -> Bond:
        """
            Returns a particular Bond.
            
              ARGUMENTS:
                - idx: which Bond to return
            
              NOTE: bond indices start at 0
            
        
        """
    def GetBonds(self) -> rdkit.Chem._GetBondsIterator:
        """
        returns an iterator over the bonds in the molecule
        """
    @typing.overload
    def GetBoolProp(self, key: str) -> typing.Any:
        """
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a bool
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetBoolProp(self, key: str, default: bool) -> typing.Any:
        """
            Returns the Bool value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: a bool, or default if the property is not present.
            
        
        """
    def GetConformer(self, id: int = -1) -> Conformer:
        """
            Get the conformer with a specified ID
        
        """
    def GetConformers(self) -> typing.Sequence[rdkit.Chem.Conformer]:
        """
            Returns a read-only sequence containing all of the molecule's Conformers.
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str) -> typing.Any:
        """
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: a double
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str, default: float) -> typing.Any:
        """
            Returns the double value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: a double, or default if the property is not present.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str) -> typing.Any:
        """
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetIntProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the integer value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: an integer, or default if the property is not present.
            
        
        """
    def GetName(self) -> str:
        """
            Returns the molecule name stored as the _Name property.
            
              NOTE:
                - If the _Name property has not been set, an empty string will be returned.
            
        
        """
    def GetNumAtoms(self, onlyHeavy: int = -1, onlyExplicit: bool = True) -> int:
        """
            Returns the number of atoms in the molecule.
            
              ARGUMENTS:
                - onlyExplicit: (optional) include only explicit atoms (atoms in the molecular graph)
                                defaults to 1.
              NOTE: the onlyHeavy argument is deprecated
            
        
        """
    def GetNumBonds(self, onlyHeavy: bool = True) -> int:
        """
            Returns the number of Bonds in the molecule.
            
              ARGUMENTS:
                - onlyHeavy: (optional) include only bonds to heavy atoms (not Hs)
                              defaults to 1.
            
        
        """
    def GetNumConformers(self) -> int:
        """
            Return the number of conformations on the molecule
        
        """
    def GetNumHeavyAtoms(self) -> int:
        """
            Returns the number of heavy atoms (atomic number >1) in the molecule.
            
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False, default: typing.Any = ...) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
                - default: value to return if the property is not present.
            
              RETURNS: the property value, or default if the property is not present.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a tuple with all property names for this molecule.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to 0.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to 0.
            
              RETURNS: a tuple of strings
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = False, includeComputed: bool = False, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary populated with properties.
            When possible, string values will be converted to integers or doubles (trimming if necessary)
             n.b. Some properties are not able to be converted to python types.
            
              ARGUMENTS:
                - includePrivate: (optional) toggles inclusion of private properties in the result set.
                                  Defaults to False.
                - includeComputed: (optional) toggles inclusion of computed properties in the result set.
                                  Defaults to False.
            
                - autoConvertStrings: (optional) toggles automatic conversion of string properties to integers or doubles.
                                  Defaults to True.
            
              RETURNS: a dictionary
            
        
        """
    def GetRingInfo(self) -> RingInfo:
        """
            Returns the number of molecule's RingInfo object.
            
            
        
        """
    def GetStereoGroups(self) -> StereoGroup_vect:
        """
            Returns a list of StereoGroups defining the relative stereochemistry of the atoms.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
            Returns the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
            Returns tuples of the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
        """
    @typing.overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns tuples of the indices of the molecule's atoms that match a substructure query.
            
              ARGUMENTS:
                - query: a Molecule.
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters) -> typing.Any:
        """
        """
    @typing.overload
    def GetUnsignedProp(self, key: str) -> typing.Any:
        """
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
              RETURNS: an unsigned integer
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str, default: int) -> typing.Any:
        """
            Returns the unsigned int value of the property if possible.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - default: value to return if the property is not present.
            
              RETURNS: an unsigned integer, or default if the property is not present.
            
        
        """
    def HasProp(self, key: str) -> int:
        """
            Queries a molecule to see if a particular property has been assigned.
            
              ARGUMENTS:
                - key: the name of the property to check for (a string).
            
        
        """
    def HasQuery(self) -> bool:
        """
            Returns if any atom or bond in molecule has a query
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
            Queries whether or not the molecule contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - recursionPossible: (optional)
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
        """
    @typing.overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> bool:
        """
            Queries whether or not the molecule contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - params: parameters controlling the substructure match
            
              RETURNS: True or False
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters = True) -> bool:
        """
        """
    def NeedsUpdatePropertyCache(self) -> bool:
        """
            Returns true or false depending on whether implicit and explicit valence of the molecule have already been calculated.
            
            
        
        """
    def RemoveAllConformers(self) -> None:
        """
            Remove all the conformations on the molecule
        
        """
    def RemoveConformer(self, id: int) -> None:
        """
            Remove the conformer with the specified ID
        
        """
    def SetBoolProp(self, key: str, val: bool, computed: bool = False) -> None:
        """
            Sets a boolean valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a bool.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetDoubleProp(self, key: str, val: float, computed: bool = False) -> None:
        """
            Sets a double valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as a double.
                - computed: (optional) marks the property as being computed.
                            Defaults to 0.
            
            
        
        """
    def SetIntProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (an unsigned number).
                - value: the property value as an integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetName(self, name: str) -> None:
        """
            Sets the molecule name; this is stored as the _Name property.
            
              ARGUMENTS:
                - name: the molecule name (a string).
            
        
        """
    def SetProp(self, key: str, val: str, computed: bool = False) -> None:
        """
            Sets a molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value (a string).
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    def SetUnsignedProp(self, key: str, val: int, computed: bool = False) -> None:
        """
            Sets an unsigned integer valued molecular property
            
              ARGUMENTS:
                - key: the name of the property to be set (a string).
                - value: the property value as an unsigned integer.
                - computed: (optional) marks the property as being computed.
                            Defaults to False.
            
            
        
        """
    @typing.overload
    def ToBinary(self) -> typing.Any:
        """
            Returns a binary string representation of the molecule.
            
        
        """
    @typing.overload
    def ToBinary(self, propertyFlags: int) -> typing.Any:
        """
            Returns a binary string representation of the molecule pickling the specified properties.
            
        
        """
    def UpdatePropertyCache(self, strict: bool = True) -> None:
        """
            Regenerates computed properties like implicit valence and ring information.
            
            
        
        """
    def __copy__(self) -> typing.Any:
        """
        """
    def __deepcopy__(self, memo: dict) -> typing.Any:
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
    def __init__(self, pklString: str) -> None:
        ...
    @typing.overload
    def __init__(self, pklString: str, propertyFlags: int) -> None:
        ...
    @typing.overload
    def __init__(self, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
class MolBundle(Boost.Python.instance):
    """
    A class for storing groups of related molecules.
    
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 88
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddMol(self, nmol: Mol) -> int:
        """
        """
    def GetMol(self, idx: int) -> Mol:
        """
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
            Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
            Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.
            
              ARGUMENTS:
                - query: a MolBundle
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a Molecule
            
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns the indices of the atoms from the first molecule in a bundle that matches a substructure query from a bundle.
            
              ARGUMENTS:
                - query: a MolBundle
            
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: Mol, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
            Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, uniquify: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000) -> typing.Any:
        """
            Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.
            
              ARGUMENTS:
                - query: a MolBundle.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: Mol, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query.
            
              ARGUMENTS:
                - query: a molecule.
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def GetSubstructMatches(self, query: MolBundle, params: SubstructMatchParameters) -> typing.Any:
        """
            Returns tuple of all indices of the atoms from the first molecule in a bundle that matches a substructure query from the second bundle.
            
              ARGUMENTS:
                - query: a MolBundle.
                - params: parameters controlling the substructure match
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
            Queries whether or not any molecule in the bundle contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - recursionPossible: (optional)
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
        """
            Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.
            
              ARGUMENTS:
                - query: a MolBundle
            
                - recursionPossible: (optional)
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: Mol, params: SubstructMatchParameters) -> bool:
        """
            Queries whether or not any molecule in the bundle contains a particular substructure.
            
              ARGUMENTS:
                - query: a Molecule
            
                - params: parameters controlling the substructure match
            
            matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: True or False
            
        
        """
    @typing.overload
    def HasSubstructMatch(self, query: MolBundle, params: SubstructMatchParameters) -> bool:
        """
            Queries whether or not any molecule in the first bundle matches any molecule in the second bundle.
            
              ARGUMENTS:
                - query: a MolBundle
            
                - params: parameters controlling the substructure match
            
              RETURNS: True or False
            
        
        """
    def Size(self) -> int:
        """
        """
    def ToBinary(self) -> typing.Any:
        """
            Returns a binary string representation of the MolBundle.
            
        
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, idx: int) -> Mol:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, pkl: str) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
class MolSanitizeException(ValueError):
    pass
class PeriodicTable(Boost.Python.instance):
    """
    A class which stores information from the Periodic Table.
    
      It is not possible to create a PeriodicTable object directly from Python,
      use GetPeriodicTable() to get the global table.
    
      The PeriodicTable object can be queried for a variety of properties:
    
        - GetAtomicWeight
    
        - GetAtomicNumber
    
        - GetElementSymbol
    
        - GetElementName
    
        - GetRow
    
        - GetRvdw (van der Waals radius)
    
        - GetRCovalent (covalent radius)
    
        - GetDefaultValence
    
        - GetValenceList
    
        - GetNOuterElecs (number of valence electrons)
    
        - GetMostCommonIsotope
    
        - GetMostCommonIsotopeMass
    
        - GetRb0
    
        - GetAbundanceForIsotope
    
        - GetMassForIsotope
    
      When it makes sense, these can be queried using either an atomic number (integer)
      or an atomic symbol (string)
    
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
    @typing.overload
    def GetAbundanceForIsotope(self, atomicNumber: int, isotope: int) -> float:
        """
        """
    @typing.overload
    def GetAbundanceForIsotope(self, elementSymbol: str, isotope: int) -> float:
        """
        """
    def GetAtomicNumber(self, elementSymbol: str) -> int:
        """
        """
    @typing.overload
    def GetAtomicWeight(self, atomicNumber: int) -> float:
        """
        """
    @typing.overload
    def GetAtomicWeight(self, elementSymbol: str) -> float:
        """
        """
    @typing.overload
    def GetDefaultValence(self, atomicNumber: int) -> int:
        """
        """
    @typing.overload
    def GetDefaultValence(self, elementSymbol: str) -> int:
        """
        """
    def GetElementName(self, atomicNumber: int) -> str:
        """
        """
    def GetElementSymbol(self, atomicNumber: int) -> str:
        """
        """
    @typing.overload
    def GetMassForIsotope(self, atomicNumber: int, isotope: int) -> float:
        """
        """
    @typing.overload
    def GetMassForIsotope(self, elementSymbol: str, isotope: int) -> float:
        """
        """
    def GetMaxAtomicNumber(self) -> int:
        """
        """
    @typing.overload
    def GetMostCommonIsotope(self, atomicNumber: int) -> int:
        """
        """
    @typing.overload
    def GetMostCommonIsotope(self, elementSymbol: str) -> int:
        """
        """
    @typing.overload
    def GetMostCommonIsotopeMass(self, atomicNumber: int) -> float:
        """
        """
    @typing.overload
    def GetMostCommonIsotopeMass(self, elementSymbol: str) -> float:
        """
        """
    @typing.overload
    def GetNOuterElecs(self, atomicNumber: int) -> int:
        """
        """
    @typing.overload
    def GetNOuterElecs(self, elementSymbol: str) -> int:
        """
        """
    @typing.overload
    def GetRb0(self, atomicNumber: int) -> float:
        """
        """
    @typing.overload
    def GetRb0(self, elementSymbol: str) -> float:
        """
        """
    @typing.overload
    def GetRcovalent(self, atomicNumber: int) -> float:
        """
        """
    @typing.overload
    def GetRcovalent(self, elementSymbol: str) -> float:
        """
        """
    @typing.overload
    def GetRow(self, atomicNumber: int) -> int:
        """
        """
    @typing.overload
    def GetRow(self, elementSymbol: str) -> int:
        """
        """
    @typing.overload
    def GetRvdw(self, atomicNumber: int) -> float:
        """
        """
    @typing.overload
    def GetRvdw(self, elementSymbol: str) -> float:
        """
        """
    @typing.overload
    def GetValenceList(self, atomicNumber: int) -> typing.Sequence[int]:
        """
        """
    @typing.overload
    def GetValenceList(self, elementSymbol: str) -> typing.Sequence[int]:
        """
        """
class PropertyPickleOptions(Boost.Python.enum):
    AllProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
    AtomProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
    BondProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
    ComputedProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
    CoordsAsDouble: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
    MolProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
    NoConformers: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
    NoProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
    PrivateProps: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
    QueryAtomData: typing.ClassVar[PropertyPickleOptions]  # value = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NoProps': rdkit.Chem.rdchem.PropertyPickleOptions.NoProps, 'MolProps': rdkit.Chem.rdchem.PropertyPickleOptions.MolProps, 'AtomProps': rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps, 'BondProps': rdkit.Chem.rdchem.PropertyPickleOptions.BondProps, 'QueryAtomData': rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData, 'PrivateProps': rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps, 'ComputedProps': rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps, 'AllProps': rdkit.Chem.rdchem.PropertyPickleOptions.AllProps, 'CoordsAsDouble': rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble, 'NoConformers': rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.PropertyPickleOptions.NoProps, 1: rdkit.Chem.rdchem.PropertyPickleOptions.MolProps, 2: rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData, 4: rdkit.Chem.rdchem.PropertyPickleOptions.BondProps, 16: rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps, 32: rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps, 65535: rdkit.Chem.rdchem.PropertyPickleOptions.AllProps, 65536: rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble, 131072: rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers}
class QueryAtom(Atom):
    """
    The class to store QueryAtoms.
    These cannot currently be constructed directly from Python
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
    def ExpandQuery(self, other: QueryAtom, how: CompositeQueryType = ..., maintainOrder: bool = True) -> None:
        """
            combines the query from other with ours
        
        """
    def SetQuery(self, other: QueryAtom) -> None:
        """
            Replace our query with a copy of the other query
        
        """
class QueryBond(Bond):
    """
    The class to store QueryBonds.
    These cannot currently be constructed directly from Python
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
    def ExpandQuery(self, other: QueryBond, how: CompositeQueryType = ..., maintainOrder: bool = True) -> None:
        """
            combines the query from other with ours
        
        """
    def SetQuery(self, other: QueryBond) -> None:
        """
            Replace our query with a copy of the other query
        
        """
class RWMol(Mol):
    """
    The RW molecule class (read/write)
    
      This class is a more-performant version of the EditableMolecule class in that
      it is a 'live' molecule and shares the interface from the Mol class.
      All changes are performed without the need to create a copy of the
      molecule using GetMol() (this is still available, however).
      
      n.b. Eventually this class may become a direct replacement for EditableMol
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 288
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddAtom(self, atom: Atom) -> int:
        """
            add an atom, returns the index of the newly added atom
        
        """
    def AddBond(self, beginAtomIdx: int, endAtomIdx: int, order: BondType = ...) -> int:
        """
            add a bond, returns the new number of bonds
        
        """
    def BeginBatchEdit(self) -> None:
        """
            starts batch editing
        
        """
    def CommitBatchEdit(self) -> None:
        """
            finishes batch editing and makes the actual changes
        
        """
    def GetMol(self) -> Mol:
        """
            Returns a Mol (a normal molecule)
        
        """
    def InsertMol(self, mol: Mol) -> None:
        """
            Insert (add) the given molecule into this one
        
        """
    def RemoveAtom(self, idx: int) -> None:
        """
            Remove the specified atom from the molecule
        
        """
    def RemoveBond(self, idx1: int, idx2: int) -> None:
        """
            Remove the specified bond from the molecule
        
        """
    def ReplaceAtom(self, index: int, newAtom: Atom, updateLabel: bool = False, preserveProps: bool = False) -> None:
        """
            replaces the specified atom with the provided one
            If updateLabel is True, the new atom becomes the active atom
            If preserveProps is True preserve keep the existing props unless explicit set on the new atom
        
        """
    def ReplaceBond(self, index: int, newBond: Bond, preserveProps: bool = False, keepSGroups: bool = True) -> None:
        """
            replaces the specified bond with the provided one.
            If preserveProps is True preserve keep the existing props unless explicit set on the new bond. If keepSGroups is False, allSubstance Groups referencing the bond will be dropped.
        
        """
    def RollbackBatchEdit(self) -> None:
        """
            cancels batch editing
        
        """
    def SetStereoGroups(self, stereo_groups: list) -> None:
        """
            Set the stereo groups
        
        """
    def __copy__(self) -> typing.Any:
        """
        """
    def __deepcopy__(self, memo: dict) -> typing.Any:
        """
        """
    def __enter__(self) -> RWMol:
        """
        """
    def __exit__(self, exc_type: typing.Any, exc_value: typing.Any, traceback: typing.Any) -> bool:
        """
        """
    def __getinitargs__(self: Mol) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __init__(self, m: Mol) -> None:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, pklString: str) -> None:
        ...
    @typing.overload
    def __init__(self, pklString: str, propertyFlags: int) -> None:
        ...
    @typing.overload
    def __init__(self, mol: Mol, quickCopy: bool = False, confId: int = -1) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
class ResonanceFlags(Boost.Python.enum):
    ALLOW_CHARGE_SEPARATION: typing.ClassVar[ResonanceFlags]  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
    ALLOW_INCOMPLETE_OCTETS: typing.ClassVar[ResonanceFlags]  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
    KEKULE_ALL: typing.ClassVar[ResonanceFlags]  # value = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
    UNCONSTRAINED_ANIONS: typing.ClassVar[ResonanceFlags]  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
    UNCONSTRAINED_CATIONS: typing.ClassVar[ResonanceFlags]  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'ALLOW_INCOMPLETE_OCTETS': rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS, 'ALLOW_CHARGE_SEPARATION': rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION, 'KEKULE_ALL': rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL, 'UNCONSTRAINED_CATIONS': rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS, 'UNCONSTRAINED_ANIONS': rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS, 2: rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION, 4: rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL, 8: rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS, 16: rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS}
class ResonanceMolSupplier(Boost.Python.instance):
    """
    A class which supplies resonance structures (as mols) from a mol.
    
      Usage examples:
    
        1) Lazy evaluation: the resonance structures are not constructed
           until we ask for them:
    
           >>> suppl = ResonanceMolSupplier(mol)
           >>> for resMol in suppl:
           ...    resMol.GetNumAtoms()
    
        2) Lazy evaluation 2:
    
           >>> suppl = ResonanceMolSupplier(mol)
           >>> resMol1 = next(suppl)
           >>> resMol2 = next(suppl)
           >>> suppl.reset()
           >>> resMol3 = next(suppl)
           # resMol3 and resMol1 are the same: 
           >>> MolToSmiles(resMol3)==MolToSmiles(resMol1)
    
        3) Random Access:
    
           >>> suppl = ResonanceMolSupplier(mol)
           >>> resMol1 = suppl[0] 
           >>> resMol2 = suppl[1] 
    
           NOTE: this will generate an IndexError if the supplier doesn't have that many
           molecules.
    
        4) Random Access 2: looping over all resonance structures
           >>> suppl = ResonanceMolSupplier(mol)
           >>> nResMols = len(suppl)
           >>> for i in range(nResMols):
           ...   suppl[i].GetNumAtoms()
    
    """
    __instance_size__: typing.ClassVar[int] = 168
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def Enumerate(self) -> None:
        """
            Ask ResonanceMolSupplier to enumerate resonance structures(automatically done as soon as any attempt to access them is made).
            
        
        """
    def GetAtomConjGrpIdx(self, ai: int) -> int:
        """
            Given an atom index, it returns the index of the conjugated groupthe atom belongs to, or -1 if it is not conjugated.
            
        
        """
    def GetBondConjGrpIdx(self, bi: int) -> int:
        """
            Given a bond index, it returns the index of the conjugated groupthe bond belongs to, or -1 if it is not conjugated.
            
        
        """
    def GetIsEnumerated(self) -> bool:
        """
            Returns true if resonance structure enumeration has already happened.
            
        
        """
    def GetNumConjGrps(self) -> int:
        """
            Returns the number of individual conjugated groups in the molecule.
            
        
        """
    def GetProgressCallback(self) -> typing.Any:
        """
            Get the ResonanceMolSupplierCallback subclass instance,
            or None if none was set.
            
        
        """
    def GetSubstructMatch(self, query: Mol, useChirality: bool = False, useQueryQueryMatches: bool = False) -> typing.Any:
        """
            Returns the indices of the molecule's atoms that match a substructure query,
            taking into account all resonance structures in ResonanceMolSupplier.
            
              ARGUMENTS:
                - query: a Molecule
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
              RETURNS: a tuple of integers
            
              NOTES:
                 - only a single match is returned
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    def GetSubstructMatches(self, query: Mol, uniquify: bool = False, useChirality: bool = False, useQueryQueryMatches: bool = False, maxMatches: int = 1000, numThreads: int = 1) -> typing.Any:
        """
            Returns tuples of the indices of the molecule's atoms that match a substructure query,
            taking into account all resonance structures in ResonanceMolSupplier.
            
              ARGUMENTS:
                - query: a Molecule.
                - uniquify: (optional) determines whether or not the matches are uniquified.
                            Defaults to 1.
            
                - useChirality: enables the use of stereochemistry in the matching
            
                - useQueryQueryMatches: use query-query matching logic
            
                - maxMatches: The maximum number of matches that will be returned.
                              In high-symmetry cases with medium-sized molecules, it is
                              very easy to end up with a combinatorial explosion in the
                              number of possible matches. This argument prevents that from
                              having unintended consequences
            
                - numThreads: The number of threads to be used (defaults to 1; 0 selects the
                              number of concurrent threads supported by the hardware; negative
                              values are added to the number of concurrent threads supported
                              by the hardware).
            
              RETURNS: a tuple of tuples of integers
            
              NOTE:
                 - the ordering of the indices corresponds to the atom ordering
                     in the query. For example, the first index is for the atom in
                     this molecule that matches the first atom in the query.
            
        
        """
    def SetNumThreads(self, numThreads: int) -> None:
        """
            Sets the number of threads to be used to enumerate resonance
            structures (defaults to 1; 0 selects the number of concurrent
            threads supported by the hardware; negative values are added
            to the number of concurrent threads supported by the hardware).
            
        
        """
    def SetProgressCallback(self, callback: typing.Any) -> None:
        """
            Pass an instance of a class derived from
            ResonanceMolSupplierCallback, which must implement the
            __call__() method.
            
        
        """
    def WasCanceled(self) -> bool:
        """
            Returns True if the resonance structure generation was canceled.
            
        
        """
    def __getitem__(self, idx: int) -> Mol:
        """
        """
    def __init__(self, mol: Mol, flags: int = 0, maxStructs: int = 1000) -> None:
        ...
    def __iter__(self) -> ResonanceMolSupplier:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __next__(self) -> Mol:
        """
            Returns the next resonance structure in the supplier. Raises _StopIteration_ on end.
            
        
        """
    def atEnd(self) -> bool:
        """
            Returns whether or not we have hit the end of the resonance structure supplier.
            
        
        """
    def reset(self) -> None:
        """
            Resets our position in the resonance structure supplier to the beginning.
            
        
        """
class ResonanceMolSupplierCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
        implement the __call__() method.
        The __call__() method is called at each iteration of the
        algorithm, and provides a mechanism to monitor or stop
        its progress.
    
        To have your callback called, pass an instance of your
        derived class to ResonanceMolSupplier.SetProgressCallback()
    """
    __instance_size__: typing.ClassVar[int] = 88
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetMaxStructures(self) -> int:
        """
            Get the number of conjugated groups this molecule has.
            
        
        """
    def GetNumConjGrps(self) -> int:
        """
            Returns the number of individual conjugated groups in the molecule.
            
        
        """
    def GetNumDiverseStructures(self, conjGrpIdx: int) -> int:
        """
            Get the number of non-degenrate resonance structures generated so far for the passed conjugated group index.
            
        
        """
    def GetNumStructures(self, conjGrpIdx: int) -> int:
        """
            Get the number of resonance structures generated so far for the passed conjugated group index.
            
        
        """
    @typing.overload
    def __call__(self) -> bool:
        """
            This must be implemented in the derived class. Return True if the resonance structure generation should continue; False if the resonance structure generation should stop.
            
        
        """
    @typing.overload
    def __call__(self) -> None:
        """
        """
    def __init__(self) -> None:
        ...
class RingInfo(Boost.Python.instance):
    """
    contains information about a molecule's rings
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
    def AddRing(self, atomIds: typing.Any, bondIds: typing.Any) -> None:
        """
            Adds a ring to the set. Be very careful with this operation.
        
        """
    def AreAtomsInSameRing(self, idx1: int, idx2: int) -> bool:
        """
        """
    def AreAtomsInSameRingOfSize(self, idx1: int, idx2: int, size: int) -> bool:
        """
        """
    def AreBondsInSameRing(self, idx1: int, idx2: int) -> bool:
        """
        """
    def AreBondsInSameRingOfSize(self, idx1: int, idx2: int, size: int) -> bool:
        """
        """
    def AreRingFamiliesInitialized(self) -> bool:
        """
        """
    def AreRingsFused(self, ring1Idx: int, ring2Idx: int) -> bool:
        """
        """
    def AtomMembers(self, idx: int) -> typing.Any:
        """
        """
    def AtomRingFamilies(self) -> typing.Any:
        """
        """
    def AtomRingSizes(self, idx: int) -> typing.Any:
        """
        """
    def AtomRings(self) -> typing.Any:
        """
        """
    def BondMembers(self, idx: int) -> typing.Any:
        """
        """
    def BondRingFamilies(self) -> typing.Any:
        """
        """
    def BondRingSizes(self, idx: int) -> typing.Any:
        """
        """
    def BondRings(self) -> typing.Any:
        """
        """
    def IsAtomInRingOfSize(self, idx: int, size: int) -> bool:
        """
        """
    def IsBondInRingOfSize(self, idx: int, size: int) -> bool:
        """
        """
    def IsRingFused(self, ringIdx: int) -> bool:
        """
        """
    def MinAtomRingSize(self, idx: int) -> int:
        """
        """
    def MinBondRingSize(self, idx: int) -> int:
        """
        """
    def NumAtomRings(self, idx: int) -> int:
        """
        """
    def NumBondRings(self, idx: int) -> int:
        """
        """
    def NumFusedBonds(self, ringIdx: int) -> int:
        """
        """
    def NumRelevantCycles(self) -> int:
        """
        """
    def NumRingFamilies(self) -> int:
        """
        """
    def NumRings(self) -> int:
        """
        """
class StereoDescriptor(Boost.Python.enum):
    Bond_Cis: typing.ClassVar[StereoDescriptor]  # value = rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis
    Bond_Trans: typing.ClassVar[StereoDescriptor]  # value = rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans
    NoValue: typing.ClassVar[StereoDescriptor]  # value = rdkit.Chem.rdchem.StereoDescriptor.NoValue
    Tet_CCW: typing.ClassVar[StereoDescriptor]  # value = rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW
    Tet_CW: typing.ClassVar[StereoDescriptor]  # value = rdkit.Chem.rdchem.StereoDescriptor.Tet_CW
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NoValue': rdkit.Chem.rdchem.StereoDescriptor.NoValue, 'Tet_CW': rdkit.Chem.rdchem.StereoDescriptor.Tet_CW, 'Tet_CCW': rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW, 'Bond_Cis': rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis, 'Bond_Trans': rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.StereoDescriptor.NoValue, 1: rdkit.Chem.rdchem.StereoDescriptor.Tet_CW, 2: rdkit.Chem.rdchem.StereoDescriptor.Tet_CCW, 3: rdkit.Chem.rdchem.StereoDescriptor.Bond_Cis, 4: rdkit.Chem.rdchem.StereoDescriptor.Bond_Trans}
class StereoGroup(Boost.Python.instance):
    """
    A collection of atoms with a defined stereochemical relationship.
    
    Used to help represent a sample with unknown stereochemistry, or that is a mix
    of diastereomers.
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
    def GetAtoms(self) -> typing.Any:
        """
            access the atoms in the StereoGroup.
            
        
        """
    def GetBonds(self) -> typing.Any:
        """
            access the bonds in the StereoGroup.
            
        
        """
    def GetGroupType(self) -> StereoGroupType:
        """
            Returns the StereoGroupType.
            
        
        """
    def GetReadId(self) -> int:
        """
            return the StereoGroup's original ID.
            Note that the ID only makes sense for AND/OR groups.
            
        
        """
    def GetWriteId(self) -> int:
        """
            return the StereoGroup's ID that will be exported.
            Note that the ID only makes sense for AND/OR groups.
            
        
        """
    def SetWriteId(self, id: int) -> None:
        """
            return the StereoGroup's ID that will be exported.
            Note that the ID only makes sense for AND/OR groups.
            
        
        """
class StereoGroupType(Boost.Python.enum):
    STEREO_ABSOLUTE: typing.ClassVar[StereoGroupType]  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
    STEREO_AND: typing.ClassVar[StereoGroupType]  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
    STEREO_OR: typing.ClassVar[StereoGroupType]  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'STEREO_ABSOLUTE': rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE, 'STEREO_OR': rdkit.Chem.rdchem.StereoGroupType.STEREO_OR, 'STEREO_AND': rdkit.Chem.rdchem.StereoGroupType.STEREO_AND}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE, 1: rdkit.Chem.rdchem.StereoGroupType.STEREO_OR, 2: rdkit.Chem.rdchem.StereoGroupType.STEREO_AND}
class StereoGroup_vect(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::StereoGroup*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::StereoGroup, std::__1::allocator<RDKit::StereoGroup>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class StereoInfo(Boost.Python.instance):
    """
    Class describing stereochemistry
    """
    NOATOM: typing.ClassVar[int] = 4294967295
    __instance_size__: typing.ClassVar[int] = 72
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def centeredOn(self) -> int:
        """index of the item the stereo concerns (default: 4294967295)"""
    @centeredOn.setter
    def centeredOn(self, value: int) -> None: ...
    @property
    def controllingAtoms(*args, **kwargs):
        """
        indices of the atoms controlling the stereo
        """
    @property
    def descriptor(*args, **kwargs):
        """
        stereo descriptor
        """
    @descriptor.setter
    def descriptor(*args, **kwargs):
        ...
    @property
    def permutation(self) -> int:
        """permutation index (used for non-tetrahedral chirality) (default: 0)"""
    @permutation.setter
    def permutation(self, value: int) -> None: ...
    @property
    def specified(*args, **kwargs):
        """
        whether or not it is specified
        """
    @specified.setter
    def specified(*args, **kwargs):
        ...
    @property
    def type(*args, **kwargs):
        """
        the type of stereo
        """
    @type.setter
    def type(*args, **kwargs):
        ...
class StereoSpecified(Boost.Python.enum):
    Specified: typing.ClassVar[StereoSpecified]  # value = rdkit.Chem.rdchem.StereoSpecified.Specified
    Unknown: typing.ClassVar[StereoSpecified]  # value = rdkit.Chem.rdchem.StereoSpecified.Unknown
    Unspecified: typing.ClassVar[StereoSpecified]  # value = rdkit.Chem.rdchem.StereoSpecified.Unspecified
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Unspecified': rdkit.Chem.rdchem.StereoSpecified.Unspecified, 'Specified': rdkit.Chem.rdchem.StereoSpecified.Specified, 'Unknown': rdkit.Chem.rdchem.StereoSpecified.Unknown}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.StereoSpecified.Unspecified, 1: rdkit.Chem.rdchem.StereoSpecified.Specified, 2: rdkit.Chem.rdchem.StereoSpecified.Unknown}
class StereoType(Boost.Python.enum):
    Atom_Octahedral: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Atom_Octahedral
    Atom_SquarePlanar: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar
    Atom_Tetrahedral: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral
    Atom_TrigonalBipyramidal: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal
    Bond_Atropisomer: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Bond_Atropisomer
    Bond_Cumulene_Even: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even
    Bond_Double: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Bond_Double
    Unspecified: typing.ClassVar[StereoType]  # value = rdkit.Chem.rdchem.StereoType.Unspecified
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Unspecified': rdkit.Chem.rdchem.StereoType.Unspecified, 'Atom_Tetrahedral': rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral, 'Atom_SquarePlanar': rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar, 'Atom_TrigonalBipyramidal': rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal, 'Atom_Octahedral': rdkit.Chem.rdchem.StereoType.Atom_Octahedral, 'Bond_Double': rdkit.Chem.rdchem.StereoType.Bond_Double, 'Bond_Cumulene_Even': rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even, 'Bond_Atropisomer': rdkit.Chem.rdchem.StereoType.Bond_Atropisomer}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.StereoType.Unspecified, 1: rdkit.Chem.rdchem.StereoType.Atom_Tetrahedral, 2: rdkit.Chem.rdchem.StereoType.Atom_SquarePlanar, 3: rdkit.Chem.rdchem.StereoType.Atom_TrigonalBipyramidal, 4: rdkit.Chem.rdchem.StereoType.Atom_Octahedral, 5: rdkit.Chem.rdchem.StereoType.Bond_Double, 6: rdkit.Chem.rdchem.StereoType.Bond_Cumulene_Even, 7: rdkit.Chem.rdchem.StereoType.Bond_Atropisomer}
class SubstanceGroup(Boost.Python.instance):
    """
    A collection of atoms and bonds with associated properties
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
    def AddAtomWithBookmark(self, mark: int) -> None:
        """
        """
    def AddAtomWithIdx(self, idx: int) -> None:
        """
        """
    def AddAttachPoint(self, aIdx: int, lvIdx: int, idStr: str) -> None:
        """
        """
    def AddBondWithBookmark(self, mark: int) -> None:
        """
        """
    def AddBondWithIdx(self, idx: int) -> None:
        """
        """
    def AddBracket(self, pts: typing.Any) -> None:
        """
        """
    def AddCState(self, bondIdx: int, vector: Point3D) -> None:
        """
        """
    def AddParentAtomWithBookmark(self, mark: int) -> None:
        """
        """
    def AddParentAtomWithIdx(self, idx: int) -> None:
        """
        """
    def ClearAttachPoints(self) -> None:
        """
        """
    def ClearBrackets(self) -> None:
        """
        """
    def ClearCStates(self) -> None:
        """
        """
    def ClearProp(self, key: typing.Any) -> None:
        """
            Removes a particular property (does nothing if not set).
            
            
        
        """
    def GetAtoms(self) -> typing.Sequence[int]:
        """
            returns a list of the indices of the atoms in this SubstanceGroup
        
        """
    def GetAttachPoints(self) -> tuple:
        """
        """
    def GetBonds(self) -> typing.Sequence[int]:
        """
            returns a list of the indices of the bonds in this SubstanceGroup
        
        """
    @typing.overload
    def GetBoolProp(self, key: typing.Any) -> bool:
        """
            returns the value of a particular property
        
        """
    @typing.overload
    def GetBoolProp(self, key: str, default: bool) -> typing.Any:
        """
            returns the value of a particular property, or default if not present
        
        """
    def GetBrackets(self) -> tuple:
        """
        """
    def GetCStates(self) -> tuple:
        """
        """
    @typing.overload
    def GetDoubleProp(self, key: typing.Any) -> float:
        """
            returns the value of a particular property
        
        """
    @typing.overload
    def GetDoubleProp(self, key: str, default: float) -> typing.Any:
        """
            returns the value of a particular property, or default if not present
        
        """
    def GetIndexInMol(self) -> int:
        """
            returns the index of this SubstanceGroup in the owning molecule's list.
        
        """
    @typing.overload
    def GetIntProp(self, key: typing.Any) -> int:
        """
            returns the value of a particular property
        
        """
    @typing.overload
    def GetIntProp(self, key: str, default: int) -> typing.Any:
        """
            returns the value of a particular property, or default if not present
        
        """
    def GetOwningMol(self) -> Mol:
        """
            returns the molecule owning this SubstanceGroup
        
        """
    def GetParentAtoms(self) -> typing.Sequence[int]:
        """
            returns a list of the indices of the parent atoms in this SubstanceGroup
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
              RETURNS: a string
            
              NOTE:
                - If the property has not been set, a KeyError exception will be raised.
            
        
        """
    @typing.overload
    def GetProp(self, key: str, autoConvert: bool = False, default: typing.Any = ...) -> typing.Any:
        """
            Returns the value of the property.
            
              ARGUMENTS:
                - key: the name of the property to return (a string).
            
                - autoConvert: if True attempt to convert the property into a python object
            
                - default: value to return if the property is not present.
            
              RETURNS: the property value, or default if the property is not present.
            
        
        """
    def GetPropNames(self, includePrivate: bool = False, includeComputed: bool = False) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Returns a list of the properties set on the SubstanceGroup.
            
            
        
        """
    def GetPropsAsDict(self, includePrivate: bool = True, includeComputed: bool = True, autoConvertStrings: bool = True) -> dict:
        """
            Returns a dictionary of the properties set on the SubstanceGroup.
             n.b. some properties cannot be converted to python types.
            
        
        """
    def GetStringVectProp(self, key: typing.Any) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            returns the value of a particular property
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: typing.Any) -> int:
        """
            returns the value of a particular property
        
        """
    @typing.overload
    def GetUnsignedProp(self, key: str, default: int) -> typing.Any:
        """
            returns the value of a particular property, or default if not present
        
        """
    def GetUnsignedVectProp(self, key: typing.Any) -> typing.Sequence[int]:
        """
            returns the value of a particular property
        
        """
    def HasProp(self, key: typing.Any) -> bool:
        """
            returns whether or not a particular property exists
        
        """
    def SetAtoms(self, iterable: typing.Any) -> None:
        """
            Set the list of the indices of the atoms in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.
        
        """
    def SetBonds(self, iterable: typing.Any) -> None:
        """
            Set the list of the indices of the bonds in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.
        
        """
    def SetBoolProp(self, key: typing.Any, val: bool, computed: bool = False) -> None:
        """
            sets the value of a particular property
        
        """
    def SetDoubleProp(self, key: typing.Any, val: float, computed: bool = False) -> None:
        """
            sets the value of a particular property
        
        """
    def SetIntProp(self, key: typing.Any, val: int, computed: bool = False) -> None:
        """
            sets the value of a particular property
        
        """
    def SetParentAtoms(self, iterable: typing.Any) -> None:
        """
            Set the list of the indices of the parent atoms in this SubstanceGroup.
            Note that this does not update properties, CStates or Attachment Points.
        
        """
    def SetProp(self, key: typing.Any, val: str, computed: bool = False) -> None:
        """
            sets the value of a particular property
        
        """
    def SetUnsignedProp(self, key: typing.Any, val: int, computed: bool = False) -> None:
        """
            sets the value of a particular property
        
        """
class SubstanceGroupAttach(Boost.Python.instance):
    """
    AttachPoint for a SubstanceGroup
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def aIdx(self) -> int:
        """attachment index (default: 0)"""
    @property
    def id(self) -> str:
        """attachment id (default: '')"""
    @property
    def lvIdx(self) -> int:
        """leaving atom or index (0 for implied) (default: 0)"""
class SubstanceGroupCState(Boost.Python.instance):
    """
    CSTATE for a SubstanceGroup
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def bondIdx(self) -> int:
        """default: 0"""
    @property
    def vector(*args, **kwargs):
        ...
class SubstanceGroup_VECT(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::SubstanceGroup*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::SubstanceGroup, std::__1::allocator<RDKit::SubstanceGroup>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class SubstructMatchParameters(Boost.Python.instance):
    """
    Parameters controlling substructure matching
    """
    __instance_size__: typing.ClassVar[int] = 208
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @typing.overload
    def setExtraAtomCheckFunc(self, func: typing.Any) -> None:
        """
            allows you to provide a function that will be called
                       for each atom pair that matches during substructure searching,
                       after all other comparisons have passed.
                       The function should return true or false indicating whether or not
                       that atom-match should be accepted.
        
        """
    @typing.overload
    def setExtraAtomCheckFunc(self, atomCoordsMatcher: AtomCoordsMatcher) -> None:
        """
            allows you to provide an AtomCoordsMatcher that will be called
                       for each atom pair that matches during substructure searching,
                       after all other comparisons have passed.
        
        """
    def setExtraBondCheckFunc(self, func: typing.Any) -> None:
        """
            allows you to provide a function that will be called
                       for each bond pair that matches during substructure searching,
                       after all other comparisons have passed.
                       The function should return true or false indicating whether or not
                       that bond-match should be accepted.
        
        """
    def setExtraFinalCheck(self, func: typing.Any) -> None:
        """
            allows you to provide a function that will be called
                           with the molecule
                       and a vector of atom IDs containing a potential match.
                       The function should return true or false indicating whether or not
                       that match should be accepted.
        
        """
    @property
    def aromaticMatchesConjugated(self) -> bool:
        """aromatic and conjugated bonds match each other (default: False)"""
    @aromaticMatchesConjugated.setter
    def aromaticMatchesConjugated(self, value: bool) -> None: ...
    @property
    def aromaticMatchesSingleOrDouble(self) -> bool:
        """aromatic and single or double bonds match each other (default: False)"""
    @aromaticMatchesSingleOrDouble.setter
    def aromaticMatchesSingleOrDouble(self, value: bool) -> None: ...
    @property
    def atomProperties(*args, **kwargs):
        """
        atom properties that must be equivalent in order to match.
        """
    @atomProperties.setter
    def atomProperties(*args, **kwargs):
        ...
    @property
    def bondProperties(*args, **kwargs):
        """
        bond properties that must be equivalent in order to match.
        """
    @bondProperties.setter
    def bondProperties(*args, **kwargs):
        ...
    @property
    def extraAtomCheckOverridesDefaultCheck(self) -> bool:
        """if set, only the extraAtomCheck will be used to determine whether or not atoms match (default: False)"""
    @extraAtomCheckOverridesDefaultCheck.setter
    def extraAtomCheckOverridesDefaultCheck(self, value: bool) -> None: ...
    @property
    def extraBondCheckOverridesDefaultCheck(self) -> bool:
        """if set, only the extraBondCheck will be used to determine whether or not bonds match (default: False)"""
    @extraBondCheckOverridesDefaultCheck.setter
    def extraBondCheckOverridesDefaultCheck(self, value: bool) -> None: ...
    @property
    def maxMatches(self) -> int:
        """maximum number of matches to return (default: 1000)"""
    @maxMatches.setter
    def maxMatches(self, value: int) -> None: ...
    @property
    def maxRecursiveMatches(self) -> int:
        """maximum number of recursive matches to find (default: 1000)"""
    @maxRecursiveMatches.setter
    def maxRecursiveMatches(self, value: int) -> None: ...
    @property
    def numThreads(self) -> int:
        """number of threads to use when multi-threading is possible.0 selects the number of concurrent threads supported by thehardware. negative values are added to the number of concurrentthreads supported by the hardware. (default: 1)"""
    @numThreads.setter
    def numThreads(self, value: int) -> None: ...
    @property
    def recursionPossible(self) -> bool:
        """Allow recursive queries (default: True)"""
    @recursionPossible.setter
    def recursionPossible(self, value: bool) -> None: ...
    @property
    def specifiedStereoQueryMatchesUnspecified(self) -> bool:
        """If set, query atoms and bonds with specified stereochemistry will match atoms and bonds with unspecified stereochemistry. (default: False)"""
    @specifiedStereoQueryMatchesUnspecified.setter
    def specifiedStereoQueryMatchesUnspecified(self, value: bool) -> None: ...
    @property
    def uniquify(self) -> bool:
        """uniquify (by atom index) match results (default: True)"""
    @uniquify.setter
    def uniquify(self, value: bool) -> None: ...
    @property
    def useChirality(self) -> bool:
        """Use chirality in determining whether or not atoms/bonds match (default: False)"""
    @useChirality.setter
    def useChirality(self, value: bool) -> None: ...
    @property
    def useEnhancedStereo(self) -> bool:
        """take enhanced stereochemistry into account while doing the match. This only has an effect if useChirality is also True. (default: False)"""
    @useEnhancedStereo.setter
    def useEnhancedStereo(self, value: bool) -> None: ...
    @property
    def useGenericMatchers(self) -> bool:
        """use generic groups (=homology groups) as a post-filtering step (if any are present in the molecule) (default: False)"""
    @useGenericMatchers.setter
    def useGenericMatchers(self, value: bool) -> None: ...
    @property
    def useQueryQueryMatches(self) -> bool:
        """Consider query-query matches, not just simple matches (default: False)"""
    @useQueryQueryMatches.setter
    def useQueryQueryMatches(self, value: bool) -> None: ...
class ValenceType(Boost.Python.enum):
    EXPLICIT: typing.ClassVar[ValenceType]  # value = rdkit.Chem.rdchem.ValenceType.EXPLICIT
    IMPLICIT: typing.ClassVar[ValenceType]  # value = rdkit.Chem.rdchem.ValenceType.IMPLICIT
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'IMPLICIT': rdkit.Chem.rdchem.ValenceType.IMPLICIT, 'EXPLICIT': rdkit.Chem.rdchem.ValenceType.EXPLICIT}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdchem.ValenceType.IMPLICIT, 1: rdkit.Chem.rdchem.ValenceType.EXPLICIT}
class _ROConformerSeq(Boost.Python.instance):
    """
    Read-only sequence of conformers, not constructible from Python.
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
    def __getitem__(self, i: int) -> Conformer:
        """
        """
    def __iter__(self) -> typing.Sequence[rdkit.Chem.Conformer]:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __next__(self) -> Conformer:
        """
        """
class _ROQAtomSeq(Boost.Python.instance):
    """
    Read-only sequence of atoms matching a query, not constructible from Python.
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
    def __getitem__(self, which: int) -> Atom:
        """
        """
    def __iter__(self) -> typing.Sequence[rdkit.Chem.QueryAtom]:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __next__(self) -> Atom:
        """
        """
class _cppAtomKekulizeException(_cppMolSanitizeException):
    """
    exception arising from sanitization
    """
    @staticmethod
    def GetAtomIndices(arg1: _cppAtomKekulizeException) -> tuple:
        """
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
class _cppAtomSanitizeException(_cppMolSanitizeException):
    """
    exception arising from sanitization
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
    def GetAtomIdx(self) -> int:
        """
        """
class _cppAtomValenceException(_cppAtomSanitizeException):
    """
    exception arising from sanitization
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
class _cppMolSanitizeException(Boost.Python.instance):
    """
    exception arising from sanitization
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
    def GetType(self) -> str:
        """
        """
    def Message(self) -> str:
        """
        """
class _listN5boost10shared_ptrIN5RDKit9ConformerEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__list_iterator<boost::shared_ptr<RDKit::Conformer>, void*>> __iter__(boost::python::back_reference<std::__1::list<boost::shared_ptr<RDKit::Conformer>, std::__1::allocator<boost::shared_ptr<RDKit::Conformer>>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _listPN5RDKit4AtomE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__list_iterator<RDKit::Atom*, void*>> __iter__(boost::python::back_reference<std::__1::list<RDKit::Atom*, std::__1::allocator<RDKit::Atom*>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _listPN5RDKit4BondE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
            C++ signature :
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__list_iterator<RDKit::Bond*, void*>> __iter__(boost::python::back_reference<std::__1::list<RDKit::Bond*, std::__1::allocator<RDKit::Bond*>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
def AddMolSubstanceGroup(mol: Mol, sgroup: SubstanceGroup) -> SubstanceGroup:
    """
        adds a copy of a SubstanceGroup to a molecule, returns the new SubstanceGroup
    
    """
def ClearMolSubstanceGroups(mol: Mol) -> None:
    """
        removes all SubstanceGroups from a molecule (if any)
    
    """
def CreateMolDataSubstanceGroup(mol: Mol, fieldName: str, value: str) -> SubstanceGroup:
    """
        creates a new DATA SubstanceGroup associated with a molecule, returns the new SubstanceGroup
    
    """
def CreateMolSubstanceGroup(mol: Mol, type: str) -> SubstanceGroup:
    """
        creates a new SubstanceGroup associated with a molecule, returns the new SubstanceGroup
    
    """
def CreateStereoGroup(stereoGroupType: StereoGroupType, mol: Mol, atomIds: typing.Any = [], bondIds: typing.Any = [], readId: int = 0) -> StereoGroup:
    """
        creates a StereoGroup associated with a molecule from a list of atom Ids
    
    """
def ForwardStereoGroupIds(mol: Mol) -> None:
    """
        Forward the original Stereo Group IDs when exporting the Mol.
    
    """
def GetAtomAlias(atom: Atom) -> str:
    """
        Returns the atom's MDL alias text
    
    """
def GetAtomRLabel(atom: Atom) -> int:
    """
        Returns the atom's MDL AtomRLabel (this is an integer from 0 to 99)
    
    """
def GetAtomValue(atom: Atom) -> str:
    """
        Returns the atom's MDL alias text
    
    """
def GetDefaultPickleProperties() -> int:
    """
        Get the current global mol pickler options.
    
    """
def GetMolSubstanceGroupWithIdx(mol: Mol, idx: int) -> SubstanceGroup:
    """
        returns a particular SubstanceGroup from the molecule
    
    """
def GetMolSubstanceGroups(mol: Mol) -> SubstanceGroup_VECT:
    """
        returns a copy of the molecule's SubstanceGroups (if any)
    
    """
def GetNumPiElectrons(atom: Atom) -> int:
    """
        Returns the number of electrons an atom is using for pi bonding
    
    """
def GetPeriodicTable() -> PeriodicTable:
    """
        Returns the application's PeriodicTable instance.
        
        
    
    """
def GetSupplementalSmilesLabel(atom: Atom) -> str:
    """
        Gets the supplemental smiles label on an atom, returns an empty string if not present.
    
    """
def MolBundleCanSerialize() -> bool:
    """
        Returns True if the MolBundle is serializable (requires boost serialization
    
    """
def SetAtomAlias(atom: Atom, rlabel: str) -> None:
    """
        Sets the atom's MDL alias text.
        Setting to an empty string clears the alias.
    
    """
def SetAtomRLabel(atom: Atom, rlabel: int) -> None:
    """
        Sets the atom's MDL RLabel (this is an integer from 0 to 99).
        Setting to 0 clears the rlabel.
    
    """
def SetAtomValue(atom: Atom, rlabel: str) -> None:
    """
        Sets the atom's MDL alias text.
        Setting to an empty string clears the alias.
    
    """
def SetDefaultPickleProperties(arg1: int) -> None:
    """
        Set the current global mol pickler options.
    
    """
def SetSupplementalSmilesLabel(atom: Atom, label: str) -> None:
    """
        Sets a supplemental label on an atom that is written to the smiles string.
        
        >>> m = Chem.MolFromSmiles("C")
        >>> Chem.SetSupplementalSmilesLabel(m.GetAtomWithIdx(0), '<xxx>')
        >>> Chem.MolToSmiles(m)
        'C<xxx>'
        
    
    """
def _HasSubstructMatchStr(pkl: str, query: Mol, recursionPossible: bool = True, useChirality: bool = False, useQueryQueryMatches: bool = False) -> bool:
    """
        This function is included to speed substructure queries from databases, 
        it's probably not of
        general interest.
        
          ARGUMENTS:
            - pkl: a Molecule pickle
        
            - query: a Molecule
        
            - recursionPossible: (optional)
        
            - useChirality: (optional)
        
            - useQueryQueryMatches: use query-query matching logic
        
          RETURNS: True or False
        
    
    """
def tossit() -> None:
    """
    """
ALLOW_CHARGE_SEPARATION: ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_CHARGE_SEPARATION
ALLOW_INCOMPLETE_OCTETS: ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.ALLOW_INCOMPLETE_OCTETS
AllProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AllProps
AtomProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.AtomProps
BondProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.BondProps
CHI_ALLENE: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_ALLENE
CHI_OCTAHEDRAL: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OCTAHEDRAL
CHI_OTHER: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_OTHER
CHI_SQUAREPLANAR: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_SQUAREPLANAR
CHI_TETRAHEDRAL: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL
CHI_TETRAHEDRAL_CCW: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
CHI_TETRAHEDRAL_CW: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
CHI_TRIGONALBIPYRAMIDAL: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_TRIGONALBIPYRAMIDAL
CHI_UNSPECIFIED: ChiralType  # value = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
COMPOSITE_AND: CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_AND
COMPOSITE_OR: CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_OR
COMPOSITE_XOR: CompositeQueryType  # value = rdkit.Chem.rdchem.CompositeQueryType.COMPOSITE_XOR
ComputedProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.ComputedProps
CoordsAsDouble: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.CoordsAsDouble
EXPLICIT: ValenceType  # value = rdkit.Chem.rdchem.ValenceType.EXPLICIT
IMPLICIT: ValenceType  # value = rdkit.Chem.rdchem.ValenceType.IMPLICIT
KEKULE_ALL: ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.KEKULE_ALL
MolProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.MolProps
NoConformers: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoConformers
NoProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.NoProps
PrivateProps: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.PrivateProps
QueryAtomData: PropertyPickleOptions  # value = rdkit.Chem.rdchem.PropertyPickleOptions.QueryAtomData
STEREO_ABSOLUTE: StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_ABSOLUTE
STEREO_AND: StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_AND
STEREO_OR: StereoGroupType  # value = rdkit.Chem.rdchem.StereoGroupType.STEREO_OR
UNCONSTRAINED_ANIONS: ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_ANIONS
UNCONSTRAINED_CATIONS: ResonanceFlags  # value = rdkit.Chem.rdchem.ResonanceFlags.UNCONSTRAINED_CATIONS
