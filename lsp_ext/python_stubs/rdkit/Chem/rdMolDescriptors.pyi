# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing functions to compute molecular descriptors
"""
from __future__ import annotations
import rdkit.rdBase
import typing
__all__: list[str] = ['AtomPairsParameters', 'BCUT2D', 'CalcAUTOCORR2D', 'CalcAUTOCORR3D', 'CalcAsphericity', 'CalcChi0n', 'CalcChi0v', 'CalcChi1n', 'CalcChi1v', 'CalcChi2n', 'CalcChi2v', 'CalcChi3n', 'CalcChi3v', 'CalcChi4n', 'CalcChi4v', 'CalcChiNn', 'CalcChiNv', 'CalcCoulombMat', 'CalcCrippenDescriptors', 'CalcEEMcharges', 'CalcEccentricity', 'CalcExactMolWt', 'CalcFractionCSP3', 'CalcGETAWAY', 'CalcHallKierAlpha', 'CalcInertialShapeFactor', 'CalcKappa1', 'CalcKappa2', 'CalcKappa3', 'CalcLabuteASA', 'CalcMORSE', 'CalcMolFormula', 'CalcNPR1', 'CalcNPR2', 'CalcNumAliphaticCarbocycles', 'CalcNumAliphaticHeterocycles', 'CalcNumAliphaticRings', 'CalcNumAmideBonds', 'CalcNumAromaticCarbocycles', 'CalcNumAromaticHeterocycles', 'CalcNumAromaticRings', 'CalcNumAtomStereoCenters', 'CalcNumAtoms', 'CalcNumBridgeheadAtoms', 'CalcNumHBA', 'CalcNumHBD', 'CalcNumHeavyAtoms', 'CalcNumHeteroatoms', 'CalcNumHeterocycles', 'CalcNumLipinskiHBA', 'CalcNumLipinskiHBD', 'CalcNumRings', 'CalcNumRotatableBonds', 'CalcNumSaturatedCarbocycles', 'CalcNumSaturatedHeterocycles', 'CalcNumSaturatedRings', 'CalcNumSpiroAtoms', 'CalcNumUnspecifiedAtomStereoCenters', 'CalcOxidationNumbers', 'CalcPBF', 'CalcPMI1', 'CalcPMI2', 'CalcPMI3', 'CalcPhi', 'CalcRDF', 'CalcRadiusOfGyration', 'CalcSpherocityIndex', 'CalcTPSA', 'CalcWHIM', 'CustomProp_VSA_', 'DoubleCubicLatticeVolume', 'GetAtomFeatures', 'GetAtomPairAtomCode', 'GetAtomPairCode', 'GetAtomPairFingerprint', 'GetConnectivityInvariants', 'GetFeatureInvariants', 'GetHashedAtomPairFingerprint', 'GetHashedAtomPairFingerprintAsBitVect', 'GetHashedMorganFingerprint', 'GetHashedTopologicalTorsionFingerprint', 'GetHashedTopologicalTorsionFingerprintAsBitVect', 'GetMACCSKeysFingerprint', 'GetMorganFingerprint', 'GetMorganFingerprintAsBitVect', 'GetTopologicalTorsionFingerprint', 'GetUSR', 'GetUSRCAT', 'GetUSRDistributions', 'GetUSRDistributionsFromPoints', 'GetUSRFromDistributions', 'GetUSRScore', 'MQNs_', 'MakePropertyRangeQuery', 'NumRotatableBondsOptions', 'PEOE_VSA_', 'Properties', 'PropertyFunctor', 'PropertyRangeQuery', 'PythonPropertyFunctor', 'SMR_VSA_', 'SlogP_VSA_']
class AtomPairsParameters(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    atomTypes: typing.ClassVar[rdkit.rdBase._vectj]  # value = <rdkit.rdBase._vectj object>
    codeSize: typing.ClassVar[int] = 9
    numAtomPairFingerprintBits: typing.ClassVar[int] = 23
    numBranchBits: typing.ClassVar[int] = 3
    numChiralBits: typing.ClassVar[int] = 2
    numPathBits: typing.ClassVar[int] = 5
    numPiBits: typing.ClassVar[int] = 2
    numTypeBits: typing.ClassVar[int] = 4
    version: typing.ClassVar[str] = '1.1.0'
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: AtomPairsParameters, arg2: str, arg3: AtomPairsParameters) -> None:
        """
        """
    def __init__(self) -> None:
        ...
class DoubleCubicLatticeVolume(Boost.Python.instance):
    """
    Class for the Double Cubic Lattice Volume method
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def GetAtomSurfaceArea(arg1: DoubleCubicLatticeVolume, atom_idx: int) -> float:
        """
            Get the surface area of atom with atom_idx
        
        """
    @staticmethod
    def GetAtomVolume(arg1: DoubleCubicLatticeVolume, atomIdx: int, solventRadius: float) -> float:
        """
            Get the volume atom of atom_idx with volume for specified Probe Radius
        
        """
    @staticmethod
    def GetCompactness(arg1: DoubleCubicLatticeVolume) -> float:
        """
            Get the Compactness of the Protein
        
        """
    @staticmethod
    def GetPackingDensity(arg1: DoubleCubicLatticeVolume) -> float:
        """
            Get the PackingDensity of the Protein
        
        """
    @staticmethod
    def GetPartialSurfaceArea(arg1: DoubleCubicLatticeVolume, atomIndices: AtomPairsParameters) -> float:
        """
            Get the Partial Surface Area of the Molecule or Protein for specified subset of atoms
        
        """
    @staticmethod
    def GetPartialVolume(arg1: DoubleCubicLatticeVolume, atomIdx: AtomPairsParameters) -> float:
        """
            Get the Partial Volume of the Molecule or Protein for specified subset of atoms
        
        """
    @staticmethod
    def GetPolarSurfaceArea(arg1: DoubleCubicLatticeVolume, includeSandP: bool = False, includeHs: bool = False) -> float:
        """
            Get the Polar Surface Area of the Molecule or Protein
        
        """
    @staticmethod
    def GetPolarVolume(arg1: DoubleCubicLatticeVolume, includeSandP: bool = False, includeHs: bool = False) -> float:
        """
            Get the Polar Volume of the Molecule or Protein
        
        """
    @staticmethod
    def GetSurfacePoints(arg1: DoubleCubicLatticeVolume) -> dict:
        """
            Get the set of points representing the surface
        
        """
    @staticmethod
    def GetVDWVolume(arg1: DoubleCubicLatticeVolume) -> float:
        """
            Get the van der Waals Volume of the Molecule or Protein
        
        """
    @staticmethod
    def GetVolume(arg1: DoubleCubicLatticeVolume) -> float:
        """
            Get the Total Volume of the Molecule or Protein
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetSurfaceArea(self) -> float:
        """
            Get the Surface Area of the Molecule or Protein
        
        """
    @typing.overload
    def __init__(self, mol: Mol, isProtein: bool = False, includeLigand: bool = True, probeRadius: float = 1.4, confId: int = -1) -> None:
        ...
    @typing.overload
    def __init__(self: AtomPairsParameters, mol: Mol, radii: list, isProtein: bool = False, includeLigand: bool = True, probeRadius: float = 1.4, confId: int = -1) -> typing.Any:
        ...
class NumRotatableBondsOptions(Boost.Python.enum):
    """
    Options for generating rotatable bonds
     NonStrict - standard loose definitions
     Strict - stricter definition excluding amides, esters, etc
     StrictLinkages - adds rotors between rotatable bonds
     Default - Current RDKit default
    """
    Default: typing.ClassVar[NumRotatableBondsOptions]  # value = rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default
    NonStrict: typing.ClassVar[NumRotatableBondsOptions]  # value = rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.NonStrict
    Strict: typing.ClassVar[NumRotatableBondsOptions]  # value = rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Strict
    StrictLinkages: typing.ClassVar[NumRotatableBondsOptions]  # value = rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NonStrict': rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.NonStrict, 'Strict': rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Strict, 'StrictLinkages': rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages, 'Default': rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.NonStrict, 1: rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Strict, 2: rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.StrictLinkages, -1: rdkit.Chem.rdMolDescriptors.NumRotatableBondsOptions.Default}
class Properties(Boost.Python.instance):
    """
    Property computation and registry system.  To compute all registered properties:
    mol = Chem.MolFromSmiles('c1ccccc1')
    properties = rdMolDescriptors.Properties()
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
      print(name, value)
    
    To compute a subset
    properties = rdMolDescriptors.Properties(['exactmw', 'lipinskiHBA'])
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(mol)):
      print(name, value)
    
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def GetAvailableProperties() -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Return all available property names that can be computed
        
        """
    @staticmethod
    def GetProperty(propName: str) -> PropertyFunctor:
        """
            Return the named property if it exists
        
        """
    @staticmethod
    def RegisterProperty(propertyFunctor: AtomPairsParameters) -> int:
        """
            Register a new property object (not thread safe)
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AnnotateProperties(self, mol: Mol) -> None:
        """
            Annotate the molecule with the computed properties.  These properties will be available as SDData or from mol.GetProp(prop)
        
        """
    def ComputeProperties(self, mol: Mol, annotateMol: bool = False) -> typing.Sequence[double]:
        """
            Return a list of computed properties, if annotateMol==True, annotate the molecule with the computed properties.
        
        """
    def GetPropertyNames(self) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Return the property names computed by this instance
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, propNames: _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE) -> None:
        ...
class PropertyFunctor(Boost.Python.instance):
    """
    Property computation class stored in the property registry.
    See rdkit.Chem.rdMolDescriptor.Properties.GetProperty and 
    rdkit.Chem.Descriptor.Properties.PropertyFunctor for creating new ones
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
    def GetName(self) -> str:
        """
            Return the name of the property to calculate
        
        """
    def GetVersion(self) -> str:
        """
            Return the version of the calculated property
        
        """
    def __call__(self, mol: Mol) -> float:
        """
            Compute the property for the specified molecule
        
        """
class PropertyRangeQuery(Boost.Python.instance):
    """
    Property Range Query for a molecule.  Match(mol) -> true if in range
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
    def Match(self, what: Mol) -> bool:
        """
        """
class PythonPropertyFunctor(PropertyFunctor):
    """
    """
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, mol: Mol) -> float:
        """
            Compute the property for the specified molecule
        
        """
    def __init__(self, callback: typing.Any, name: str, version: str) -> None:
        ...
@typing.overload
def BCUT2D(mol: Mol) -> list:
    """
        Implements BCUT descriptors From J. Chem. Inf. Comput. Sci., Vol. 39, No. 1, 1999Diagonal elements are (currently) atomic mass, gasteiger charge,crippen logP and crippen MRReturns the 2D BCUT2D descriptors vector as described in
        returns [mass eigen value high, mass eigen value low,
                 gasteiger charge eigenvalue high, gasteiger charge low,
                 crippen lowgp  eigenvalue high, crippen lowgp  low,
                 crippen mr eigenvalue high, crippen mr low]
        
    
    """
@typing.overload
def BCUT2D(mol: Mol, atom_props: list) -> tuple:
    """
        Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule and the specified atom props
         atom_props must be a list or tuple of floats equal in size to the number of atoms in mol
    
    """
@typing.overload
def BCUT2D(mol: Mol, atom_props: tuple) -> tuple:
    """
        Returns a 2D BCUT (eigen value hi, eigenvalue low) given the molecule and the specified atom props
         atom_props must be a list or tuple of floats equal in size to the number of atoms in mol
    
    """
@typing.overload
def BCUT2D(mol: Mol, atom_propname: str) -> tuple:
    """
        Returns a 2D BCUT (eigen value high, eigen value low) given the molecule and the specified atom prop name
        atom_propname must exist on each atom and be convertible to a float
    
    """
def CalcAUTOCORR2D(mol: Mol, CustomAtomProperty: str = '') -> list:
    """
        Returns 2D Autocorrelation descriptors vector
    
    """
def CalcAUTOCORR3D(mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
        Returns 3D Autocorrelation descriptors vector
    
    """
def CalcAsphericity(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcChi0n(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi0v(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi1n(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi1v(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi2n(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi2v(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi3n(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi3v(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi4n(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChi4v(mol: Mol, force: bool = False) -> float:
    """
    """
def CalcChiNn(mol: Mol, n: int, force: bool = False) -> float:
    """
    """
def CalcChiNv(mol: Mol, n: int, force: bool = False) -> float:
    """
    """
def CalcCoulombMat(mol: Mol, confId: int = -1) -> tuple:
    """
        Returns severals Coulomb randomized matrices
    
    """
def CalcCrippenDescriptors(mol: Mol, includeHs: bool = True, force: bool = False) -> tuple:
    """
        returns a 2-tuple with the Wildman-Crippen logp,mr values
    
    """
def CalcEEMcharges(mol: Mol, confId: int = -1) -> list:
    """
        Returns EEM atomic partial charges
    
    """
def CalcEccentricity(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcExactMolWt(mol: Mol, onlyHeavy: bool = False) -> float:
    """
        returns the molecule's exact molecular weight
    
    """
def CalcFractionCSP3(mol: Mol) -> float:
    """
        returns the fraction of C atoms that are SP3 hybridized
    
    """
def CalcGETAWAY(mol: Mol, confId: int = -1, precision: float = 2, CustomAtomProperty: str = '') -> list:
    """
        Returns the GETAWAY descriptors vector
    
    """
def CalcHallKierAlpha(mol: Mol, atomContribs: AtomPairsParameters = None) -> float:
    """
    """
def CalcInertialShapeFactor(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcKappa1(mol: Mol) -> float:
    """
    """
def CalcKappa2(mol: Mol) -> float:
    """
    """
def CalcKappa3(mol: Mol) -> float:
    """
    """
def CalcLabuteASA(mol: Mol, includeHs: bool = True, force: bool = False) -> float:
    """
        returns the Labute ASA value for a molecule
    
    """
def CalcMORSE(mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
        Returns Molecule Representation of Structures based on Electron diffraction descriptors
    
    """
def CalcMolFormula(mol: Mol, separateIsotopes: bool = False, abbreviateHIsotopes: bool = True) -> str:
    """
        returns the molecule's formula
    
    """
def CalcNPR1(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcNPR2(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcNumAliphaticCarbocycles(mol: Mol) -> int:
    """
        returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule
    
    """
def CalcNumAliphaticHeterocycles(mol: Mol) -> int:
    """
        returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule
    
    """
def CalcNumAliphaticRings(mol: Mol) -> int:
    """
        returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule
    
    """
def CalcNumAmideBonds(mol: Mol) -> int:
    """
        returns the number of amide bonds in a molecule
    
    """
def CalcNumAromaticCarbocycles(mol: Mol) -> int:
    """
        returns the number of aromatic carbocycles for a molecule
    
    """
def CalcNumAromaticHeterocycles(mol: Mol) -> int:
    """
        returns the number of aromatic heterocycles for a molecule
    
    """
def CalcNumAromaticRings(mol: Mol) -> int:
    """
        returns the number of aromatic rings for a molecule
    
    """
def CalcNumAtomStereoCenters(mol: Mol) -> int:
    """
        Returns the total number of atomic stereocenters (specified and unspecified)
    
    """
def CalcNumAtoms(mol: Mol) -> int:
    """
        returns the total number of atoms for a molecule
    
    """
def CalcNumBridgeheadAtoms(mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
        Returns the number of bridgehead atoms (atoms shared between rings that share at least two bonds)
    
    """
def CalcNumHBA(mol: Mol) -> int:
    """
        returns the number of H-bond acceptors for a molecule
    
    """
def CalcNumHBD(mol: Mol) -> int:
    """
        returns the number of H-bond donors for a molecule
    
    """
def CalcNumHeavyAtoms(mol: Mol) -> int:
    """
        returns the number of heavy atoms for a molecule
    
    """
def CalcNumHeteroatoms(mol: Mol) -> int:
    """
        returns the number of heteroatoms for a molecule
    
    """
def CalcNumHeterocycles(mol: Mol) -> int:
    """
        returns the number of heterocycles for a molecule
    
    """
def CalcNumLipinskiHBA(mol: Mol) -> int:
    """
        returns the number of Lipinski H-bond acceptors for a molecule
    
    """
def CalcNumLipinskiHBD(mol: Mol) -> int:
    """
        returns the number of Lipinski H-bond donors for a molecule
    
    """
def CalcNumRings(mol: Mol) -> int:
    """
        returns the number of rings for a molecule
    
    """
@typing.overload
def CalcNumRotatableBonds(mol: Mol, strict: bool) -> int:
    """
        returns the number of rotatable bonds for a molecule.
           strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
           strict = NumRotatableBondsOptions.Strict - (default) does not count things like
                    amide or ester bonds
           strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
              systems.
              - Single bonds between aliphatic ring Cs are always rotatable. This
                means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now 
                considered rotatable; it was not before
              - Heteroatoms in the linked rings no longer affect whether or not
                the linking bond is rotatable
              - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
                 considered non-rotatable
    
    """
@typing.overload
def CalcNumRotatableBonds(mol: Mol, strict: NumRotatableBondsOptions = ...) -> int:
    """
        returns the number of rotatable bonds for a molecule.
           strict = NumRotatableBondsOptions.NonStrict - Simple rotatable bond definition.
           strict = NumRotatableBondsOptions.Strict - (default) does not count things like
                    amide or ester bonds
           strict = NumRotatableBondsOptions.StrictLinkages - handles linkages between ring
              systems.
              - Single bonds between aliphatic ring Cs are always rotatable. This
                means that the central bond in CC1CCCC(C)C1-C1C(C)CCCC1C is now 
                considered rotatable; it was not before
              - Heteroatoms in the linked rings no longer affect whether or not
                the linking bond is rotatable
              - the linking bond in systems like Cc1cccc(C)c1-c1c(C)cccc1 is now
                 considered non-rotatable
    
    """
def CalcNumSaturatedCarbocycles(mol: Mol) -> int:
    """
        returns the number of saturated carbocycles for a molecule
    
    """
def CalcNumSaturatedHeterocycles(mol: Mol) -> int:
    """
        returns the number of saturated heterocycles for a molecule
    
    """
def CalcNumSaturatedRings(mol: Mol) -> int:
    """
        returns the number of saturated rings for a molecule
    
    """
def CalcNumSpiroAtoms(mol: Mol, atoms: AtomPairsParameters = None) -> int:
    """
        Returns the number of spiro atoms (atoms shared between rings that share exactly one atom)
    
    """
def CalcNumUnspecifiedAtomStereoCenters(mol: Mol) -> int:
    """
        Returns the number of unspecified atomic stereocenters
    
    """
def CalcOxidationNumbers(mol: Mol) -> None:
    """
        Adds the oxidation number/state to the atoms of a molecule as property OxidationNumber on each atom.  Use Pauling electronegativities.  This is experimental code, still under development.
    
    """
def CalcPBF(mol: Mol, confId: int = -1) -> float:
    """
        Returns the PBF (plane of best fit) descriptor (https://doi.org/10.1021/ci300293f)
    
    """
def CalcPMI1(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcPMI2(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcPMI3(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcPhi(mol: Mol) -> float:
    """
    """
def CalcRDF(mol: Mol, confId: int = -1, CustomAtomProperty: str = '') -> list:
    """
        Returns radial distribution fonction descriptors (RDF)
    
    """
def CalcRadiusOfGyration(mol: Mol, confId: int = -1, useAtomicMasses: bool = True, force: bool = True) -> float:
    """
    """
def CalcSpherocityIndex(mol: Mol, confId: int = -1, force: bool = True) -> float:
    """
    """
def CalcTPSA(mol: Mol, force: bool = False, includeSandP: bool = False) -> float:
    """
        returns the TPSA value for a molecule
    
    """
def CalcWHIM(mol: Mol, confId: int = -1, thresh: float = 0.001, CustomAtomProperty: str = '') -> list:
    """
        Returns the WHIM descriptors vector
    
    """
def CustomProp_VSA_(mol: Mol, customPropName: str, bins: AtomPairsParameters, force: bool = False) -> list:
    """
    """
def GetAtomFeatures(mol: Mol, atomid: int, addchiral: bool = False) -> list:
    """
        Returns the Atom Features vector
    
    """
def GetAtomPairAtomCode(atom: Atom, branchSubtract: int = 0, includeChirality: bool = False) -> int:
    """
        Returns the atom code (hash) for an atom
    
    """
def GetAtomPairCode(atom1Code: int, atom2Code: int, distance: int, includeChirality: bool = False) -> int:
    """
        Returns the atom-pair code (hash) for a pair of atoms separated by a certain number of bonds
    
    """
def GetAtomPairFingerprint(mol: Mol, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
        Returns the atom-pair fingerprint for a molecule as an IntSparseIntVect
    
    """
def GetConnectivityInvariants(mol: Mol, includeRingMembership: bool = True) -> list:
    """
        Returns connectivity invariants (ECFP-like) for a molecule.
    
    """
def GetFeatureInvariants(mol: Mol) -> list:
    """
        Returns feature invariants (FCFP-like) for a molecule.
    
    """
def GetHashedAtomPairFingerprint(mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> IntSparseIntVect:
    """
        Returns the hashed atom-pair fingerprint for a molecule as an IntSparseIntVect
    
    """
def GetHashedAtomPairFingerprintAsBitVect(mol: Mol, nBits: int = 2048, minLength: int = 1, maxLength: int = 30, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False, use2D: bool = True, confId: int = -1) -> ExplicitBitVect:
    """
        Returns the atom-pair fingerprint for a molecule as an ExplicitBitVect
    
    """
def GetHashedMorganFingerprint(mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect:
    """
        Returns a hashed Morgan fingerprint for a molecule
    
    """
def GetHashedTopologicalTorsionFingerprint(mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
        Returns the hashed topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
    
    """
def GetHashedTopologicalTorsionFingerprintAsBitVect(mol: Mol, nBits: int = 2048, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, nBitsPerEntry: int = 4, includeChirality: bool = False) -> ExplicitBitVect:
    """
        Returns the topological-torsion fingerprint for a molecule as an ExplicitBitVect
    
    """
def GetMACCSKeysFingerprint(mol: Mol) -> ExplicitBitVect:
    """
        Returns the MACCS keys for a molecule as an ExplicitBitVect
    
    """
def GetMorganFingerprint(mol: Mol, radius: int, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, useCounts: bool = True, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> UIntSparseIntVect:
    """
        Returns a Morgan fingerprint for a molecule
    
    """
def GetMorganFingerprintAsBitVect(mol: Mol, radius: int, nBits: int = 2048, invariants: AtomPairsParameters = [], fromAtoms: AtomPairsParameters = [], useChirality: bool = False, useBondTypes: bool = True, useFeatures: bool = False, bitInfo: AtomPairsParameters = None, includeRedundantEnvironments: bool = False) -> ExplicitBitVect:
    """
        Returns a Morgan fingerprint for a molecule as a bit vector
    
    """
def GetTopologicalTorsionFingerprint(mol: Mol, targetSize: int = 4, fromAtoms: AtomPairsParameters = 0, ignoreAtoms: AtomPairsParameters = 0, atomInvariants: AtomPairsParameters = 0, includeChirality: bool = False) -> LongSparseIntVect:
    """
        Returns the topological-torsion fingerprint for a molecule as a LongIntSparseIntVect
    
    """
def GetUSR(mol: Mol, confId: int = -1) -> list:
    """
        Returns a USR descriptor for one conformer of a molecule
    
    """
def GetUSRCAT(mol: Mol, atomSelections: AtomPairsParameters = None, confId: int = -1) -> list:
    """
        Returns a USRCAT descriptor for one conformer of a molecule
    
    """
def GetUSRDistributions(coords: AtomPairsParameters, points: AtomPairsParameters = None) -> list:
    """
        Returns the four USR distance distributions for a set of coordinates
    
    """
def GetUSRDistributionsFromPoints(coords: AtomPairsParameters, points: AtomPairsParameters) -> list:
    """
        Returns the USR distance distributions for a set of coordinates and points
    
    """
def GetUSRFromDistributions(distances: AtomPairsParameters) -> list:
    """
        Returns the USR descriptor from a set of distance distributions
    
    """
def GetUSRScore(descriptor1: AtomPairsParameters, descriptor2: AtomPairsParameters, weights: AtomPairsParameters = []) -> float:
    """
        Returns the USR score for two USR or USRCAT descriptors
    
    """
def MQNs_(mol: Mol, force: bool = False) -> list:
    """
    """
def MakePropertyRangeQuery(name: str, min: float, max: float) -> PropertyRangeQuery:
    """
        Generates a Range property for the specified property, between min and max
        query = MakePropertyRangeQuery('exactmw', 0, 500)
        query.Match( mol )
    
    """
def PEOE_VSA_(mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    """
def SMR_VSA_(mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    """
def SlogP_VSA_(mol: Mol, bins: AtomPairsParameters = [], force: bool = False) -> list:
    """
    """
def _CalcCrippenContribs(mol: Mol, force: bool = False, atomTypes: list = [], atomTypeLabels: list = []) -> list:
    """
        returns (as a list of 2-tuples) the contributions of each atom to
        the Wildman-Cripppen logp and mr value
    
    """
def _CalcLabuteASAContribs(mol: Mol, includeHs: bool = True, force: bool = False) -> tuple:
    """
        returns a list of atomic contributions to the Labute ASA
    
    """
def _CalcMolWt(mol: Mol, onlyHeavy: bool = False) -> float:
    """
        returns the molecule's molecular weight
    
    """
def _CalcTPSAContribs(mol: Mol, force: bool = False, includeSandP: bool = False) -> tuple:
    """
        returns a list of atomic contributions to the TPSA
    
    """
_BCUT2D_version: str = '1.0.0'
_CalcAUTOCORR2D_version: str = '1.0.0'
_CalcAUTOCORR3D_version: str = '1.0.0'
_CalcAsphericity_version: str = '1.0.0'
_CalcChi0n_version: str = '1.2.0'
_CalcChi0v_version: str = '1.2.0'
_CalcChi1n_version: str = '1.2.0'
_CalcChi1v_version: str = '1.2.0'
_CalcChi2n_version: str = '1.2.0'
_CalcChi2v_version: str = '1.2.0'
_CalcChi3n_version: str = '1.2.0'
_CalcChi3v_version: str = '1.2.0'
_CalcChi4n_version: str = '1.2.0'
_CalcChi4v_version: str = '1.2.0'
_CalcChiNn_version: str = '1.2.0'
_CalcChiNv_version: str = '1.2.0'
_CalcCoulombMat_version: str = '1.0.0'
_CalcCrippenDescriptors_version: str = '1.2.0'
_CalcEMMcharges_version: str = '1.0.0'
_CalcEccentricity_version: str = '1.0.0'
_CalcExactMolWt_version: str = '1.0.0'
_CalcFractionCSP3_version: str = '1.0.0'
_CalcGETAWAY_version: str = '1.0.0'
_CalcHallKierAlpha_version: str = '1.2.0'
_CalcInertialShapeFactor_version: str = '1.0.0'
_CalcKappa1_version: str = '1.1.0'
_CalcKappa2_version: str = '1.1.0'
_CalcKappa3_version: str = '1.1.0'
_CalcLabuteASA_version: str = '1.0.2'
_CalcMORSE_version: str = '1.0.0'
_CalcMolFormula_version: str = '1.3.0'
_CalcMolWt_version: str = '1.0.0'
_CalcNPR1_version: str = '1.0.0'
_CalcNPR2_version: str = '1.0.0'
_CalcNumAliphaticCarbocycles_version: str = '1.0.0'
_CalcNumAliphaticHeterocycles_version: str = '1.0.0'
_CalcNumAliphaticRings_version: str = '1.0.0'
_CalcNumAmideBonds_version: str = '1.0.0'
_CalcNumAromaticCarbocycles_version: str = '1.0.0'
_CalcNumAromaticHeterocycles_version: str = '1.0.0'
_CalcNumAromaticRings_version: str = '1.0.0'
_CalcNumAtomStereoCenters_version: str = '1.0.1'
_CalcNumAtoms_version: str = '1.0.0'
_CalcNumBridgeheadAtoms_version: str = '2.0.0'
_CalcNumHBA_version: str = '2.0.2'
_CalcNumHBD_version: str = '2.0.1'
_CalcNumHeavyAtoms_version: str = '1.0.0'
_CalcNumHeteroatoms_version: str = '1.0.1'
_CalcNumHeterocycles_version: str = '1.0.0'
_CalcNumLipinskiHBA_version: str = '1.0.0'
_CalcNumLipinskiHBD_version: str = '2.0.0'
_CalcNumRings_version: str = '1.0.1'
_CalcNumRotatableBonds_version: str = '3.2.0'
_CalcNumSaturatedCarbocycles_version: str = '1.0.0'
_CalcNumSaturatedHeterocycles_version: str = '1.0.0'
_CalcNumSaturatedRings_version: str = '1.0.0'
_CalcNumSpiroAtoms_version: str = '1.0.0'
_CalcNumUnspecifiedAtomStereoCenters_version: str = '1.0.1'
_CalcPBF_version: str = '1.0.0'
_CalcPMI1_version: str = '1.0.0'
_CalcPMI2_version: str = '1.0.0'
_CalcPMI3_version: str = '1.0.0'
_CalcPhi_version: str = '1.0.0'
_CalcRDF_version: str = '1.0.0'
_CalcRadiusOfGyration_version: str = '1.0.0'
_CalcSpherocityIndex_version: str = '1.0.0'
_CalcTPSA_version: str = '2.0.0'
_CalcWHIM_version: str = '1.0.0'
_ConnectivityInvariants_version: str = '1.0.0'
_FeatureInvariants_version: str = '0.1.0'
_GetAtomFeatures_version: str = '1.0.0'
_MorganFingerprint_version: str = '1.0.0'
