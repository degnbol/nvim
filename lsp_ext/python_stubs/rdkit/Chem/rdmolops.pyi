# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing RDKit functionality for manipulating molecules.
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['ADJUST_IGNOREALL', 'ADJUST_IGNORECHAINS', 'ADJUST_IGNOREDUMMIES', 'ADJUST_IGNOREMAPPED', 'ADJUST_IGNORENONDUMMIES', 'ADJUST_IGNORENONE', 'ADJUST_IGNORERINGS', 'AROMATICITY_CUSTOM', 'AROMATICITY_DEFAULT', 'AROMATICITY_MDL', 'AROMATICITY_MMFF94', 'AROMATICITY_RDKIT', 'AROMATICITY_SIMPLE', 'AddHs', 'AddHsParameters', 'AddRecursiveQuery', 'AddStereoAnnotations', 'AddWavyBondsForStereoAny', 'AdjustQueryParameters', 'AdjustQueryProperties', 'AdjustQueryPropertiesWithGenericGroups', 'AdjustQueryWhichFlags', 'AromaticityModel', 'AssignAtomChiralTagsFromMolParity', 'AssignAtomChiralTagsFromStructure', 'AssignChiralTypesFromBondDirs', 'AssignRadicals', 'AssignStereochemistry', 'AssignStereochemistryFrom3D', 'AtomHasConjugatedBond', 'BondWedgingParameters', 'BoolVector', 'CanonicalizeStereoGroups', 'Cleanup', 'CleanupAtropisomers', 'CleanupChirality', 'CleanupOrganometallics', 'CleanupStereoGroups', 'CollapseAttachmentPoints', 'CombineMols', 'ComputeAtomCIPRanks', 'ConvertGenericQueriesToSubstanceGroups', 'CopyMolSubset', 'CountAtomElec', 'DativeBondsToHaptic', 'DeleteSubstructs', 'DetectBondStereoChemistry', 'DetectBondStereochemistry', 'DetectChemistryProblems', 'ExpandAttachmentPoints', 'FastFindRings', 'FindAllPathsOfLengthN', 'FindAllSubgraphsOfLengthMToN', 'FindAllSubgraphsOfLengthN', 'FindAtomEnvironmentOfRadiusN', 'FindMesoCenters', 'FindPotentialStereo', 'FindPotentialStereoBonds', 'FindRingFamilies', 'FindUniqueSubgraphsOfLengthN', 'FragmentOnBRICSBonds', 'FragmentOnBonds', 'FragmentOnSomeBonds', 'Get3DDistanceMatrix', 'GetAdjacencyMatrix', 'GetAllowNontetrahedralChirality', 'GetDistanceMatrix', 'GetFormalCharge', 'GetMolFrags', 'GetMostSubstitutedCoreMatch', 'GetSSSR', 'GetShortestPath', 'GetSymmSSSR', 'GetUseLegacyStereoPerception', 'HapticBondsToDative', 'HasQueryHs', 'Kekulize', 'KekulizeIfPossible', 'LayeredFingerprint', 'LayeredFingerprint_substructLayers', 'MergeQueryHs', 'MolAddRecursiveQueries', 'MolzipLabel', 'MolzipParams', 'MurckoDecompose', 'NeedsHs', 'ParseMolQueryDefFile', 'PathToSubmol', 'PatternFingerprint', 'RDKFingerprint', 'ReapplyMolBlockWedging', 'RemoveAllHs', 'RemoveHs', 'RemoveHsParameters', 'RemoveNonExplicit3DChirality', 'RemoveStereochemistry', 'RenumberAtoms', 'ReplaceCore', 'ReplaceSidechains', 'ReplaceSubstructs', 'SANITIZE_ADJUSTHS', 'SANITIZE_ALL', 'SANITIZE_CLEANUP', 'SANITIZE_CLEANUPATROPISOMERS', 'SANITIZE_CLEANUPCHIRALITY', 'SANITIZE_CLEANUP_ORGANOMETALLICS', 'SANITIZE_FINDRADICALS', 'SANITIZE_KEKULIZE', 'SANITIZE_NONE', 'SANITIZE_PROPERTIES', 'SANITIZE_SETAROMATICITY', 'SANITIZE_SETCONJUGATION', 'SANITIZE_SETHYBRIDIZATION', 'SANITIZE_SYMMRINGS', 'SanitizeFlags', 'SanitizeMol', 'SetAllowNontetrahedralChirality', 'SetAromaticity', 'SetBondStereoFromDirections', 'SetConjugation', 'SetDoubleBondNeighborDirections', 'SetGenericQueriesFromProperties', 'SetHybridization', 'SetTerminalAtomCoords', 'SetUseLegacyStereoPerception', 'SimplifyEnhancedStereo', 'SortMatchesByDegreeOfCoreSubstitution', 'SplitMolByPDBChainId', 'SplitMolByPDBResidues', 'StereoBondThresholds', 'StereoGroupAbsOptions', 'SubsetInfo', 'SubsetMethod', 'SubsetOptions', 'TranslateChiralFlagToStereoGroups', 'UIntUIntMap', 'UnfoldedRDKFingerprintCountBased', 'WedgeBond', 'WedgeMolBonds', 'map_indexing_suite_UIntUIntMap_entry', 'molzip', 'molzipFragments']
class AddHsParameters(Boost.Python.instance):
    """
    Parameters controlling H addition.
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def addCoords(self) -> bool:
        """add coordinates for the Hs (default: False)"""
    @addCoords.setter
    def addCoords(self, value: bool) -> None: ...
    @property
    def addResidueInfo(self) -> bool:
        """add residue info to the Hs (default: False)"""
    @addResidueInfo.setter
    def addResidueInfo(self, value: bool) -> None: ...
    @property
    def explicitOnly(self) -> bool:
        """only add explict Hs (default: False)"""
    @explicitOnly.setter
    def explicitOnly(self, value: bool) -> None: ...
    @property
    def skipQueries(self) -> bool:
        """do not add Hs to query atoms or atoms with query bonds (default: False)"""
    @skipQueries.setter
    def skipQueries(self, value: bool) -> None: ...
class AdjustQueryParameters(Boost.Python.instance):
    """
    Parameters controlling which components of the query atoms/bonds are adjusted.
    
    Note that some of the options here are either directly contradictory or make
      no sense when combined with each other. We generally assume that client code
      is doing something sensible and don't attempt to detect possible conflicts or
      problems.
    
    A note on the flags controlling which atoms/bonds are modified: 
       These generally limit the set of atoms/bonds to be modified.
       For example:
           - ADJUST_IGNORERINGS atoms/bonds in rings will not be modified.
           - ADJUST_IGNORENONE causes all atoms/bonds to be modified
           - ADJUST_IGNOREALL no atoms/bonds will be modified
       Some of the options obviously make no sense for bonds
    """
    __instance_size__: typing.ClassVar[int] = 80
    @staticmethod
    def NoAdjustments() -> AdjustQueryParameters:
        """
            Returns an AdjustQueryParameters object with all parameters set to false
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def adjustConjugatedFiveRings(self) -> bool:
        """set bond queries in conjugated five-rings to SINGLE|DOUBLE|AROMATIC (default: False)"""
    @adjustConjugatedFiveRings.setter
    def adjustConjugatedFiveRings(self, value: bool) -> None: ...
    @property
    def adjustDegree(self) -> bool:
        """add degree queries (default: True)"""
    @adjustDegree.setter
    def adjustDegree(self, value: bool) -> None: ...
    @property
    def adjustDegreeFlags(self) -> int:
        """controls which atoms have their degree queries changed (default: 3)"""
    @adjustDegreeFlags.setter
    def adjustDegreeFlags(self, value: int) -> None: ...
    @property
    def adjustHeavyDegree(self) -> bool:
        """adjust the heavy-atom degree (default: False)"""
    @adjustHeavyDegree.setter
    def adjustHeavyDegree(self, value: bool) -> None: ...
    @property
    def adjustHeavyDegreeFlags(self) -> int:
        """controls which atoms have their heavy-atom degree queries changed (default: 3)"""
    @adjustHeavyDegreeFlags.setter
    def adjustHeavyDegreeFlags(self, value: int) -> None: ...
    @property
    def adjustRingChain(self) -> bool:
        """add ring-chain queries to atoms (default: False)"""
    @adjustRingChain.setter
    def adjustRingChain(self, value: bool) -> None: ...
    @property
    def adjustRingChainFlags(self) -> int:
        """controls which atoms have ring-chain queries added (default: 0)"""
    @adjustRingChainFlags.setter
    def adjustRingChainFlags(self, value: int) -> None: ...
    @property
    def adjustRingCount(self) -> bool:
        """add ring-count queries (default: False)"""
    @adjustRingCount.setter
    def adjustRingCount(self, value: bool) -> None: ...
    @property
    def adjustRingCountFlags(self) -> int:
        """controls which atoms have ring-count queries added (default: 3)"""
    @adjustRingCountFlags.setter
    def adjustRingCountFlags(self, value: int) -> None: ...
    @property
    def adjustSingleBondsBetweenAromaticAtoms(self) -> bool:
        """sets non-ring single bonds between two aromatic or conjugated atoms to SINGLE|AROMATIC (default: False)"""
    @adjustSingleBondsBetweenAromaticAtoms.setter
    def adjustSingleBondsBetweenAromaticAtoms(self, value: bool) -> None: ...
    @property
    def adjustSingleBondsToDegreeOneNeighbors(self) -> bool:
        """set single bonds bewteen aromatic or conjugated atoms and degree-one neighbors to SINGLE|AROMATIC (default: False)"""
    @adjustSingleBondsToDegreeOneNeighbors.setter
    def adjustSingleBondsToDegreeOneNeighbors(self, value: bool) -> None: ...
    @property
    def aromatizeIfPossible(self) -> bool:
        """perceive and set aromaticity (default: True)"""
    @aromatizeIfPossible.setter
    def aromatizeIfPossible(self, value: bool) -> None: ...
    @property
    def makeAtomsGeneric(self) -> bool:
        """convert atoms to generic queries (any atoms) (default: False)"""
    @makeAtomsGeneric.setter
    def makeAtomsGeneric(self, value: bool) -> None: ...
    @property
    def makeAtomsGenericFlags(self) -> int:
        """controls which atoms are converted to generic queries (default: 0)"""
    @makeAtomsGenericFlags.setter
    def makeAtomsGenericFlags(self, value: int) -> None: ...
    @property
    def makeBondsGeneric(self) -> bool:
        """converts bonds to generic queries (any bonds) (default: False)"""
    @makeBondsGeneric.setter
    def makeBondsGeneric(self, value: bool) -> None: ...
    @property
    def makeBondsGenericFlags(self) -> int:
        """controls which bonds are converted to generic queries (default: 0)"""
    @makeBondsGenericFlags.setter
    def makeBondsGenericFlags(self, value: int) -> None: ...
    @property
    def makeDummiesQueries(self) -> bool:
        """convert dummy atoms without isotope labels to any-atom queries (default: True)"""
    @makeDummiesQueries.setter
    def makeDummiesQueries(self, value: bool) -> None: ...
    @property
    def setMDLFiveRingAromaticity(self) -> bool:
        """uses the 5-ring aromaticity behavior of the (former) MDL software as documented in the Chemical Representation Guide (default: False)"""
    @setMDLFiveRingAromaticity.setter
    def setMDLFiveRingAromaticity(self, value: bool) -> None: ...
    @property
    def useStereoCareForBonds(self) -> bool:
        """if this is set sterochemistry information will be removed from double bonds that do not have the stereoCare property set (default: False)"""
    @useStereoCareForBonds.setter
    def useStereoCareForBonds(self, value: bool) -> None: ...
class AdjustQueryWhichFlags(Boost.Python.enum):
    ADJUST_IGNOREALL: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL
    ADJUST_IGNORECHAINS: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
    ADJUST_IGNOREDUMMIES: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
    ADJUST_IGNOREMAPPED: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED
    ADJUST_IGNORENONDUMMIES: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES
    ADJUST_IGNORENONE: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE
    ADJUST_IGNORERINGS: typing.ClassVar[AdjustQueryWhichFlags]  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'ADJUST_IGNORENONE': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE, 'ADJUST_IGNORECHAINS': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS, 'ADJUST_IGNORERINGS': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS, 'ADJUST_IGNOREDUMMIES': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES, 'ADJUST_IGNORENONDUMMIES': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES, 'ADJUST_IGNOREMAPPED': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED, 'ADJUST_IGNOREALL': rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE, 1: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS, 4: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS, 2: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES, 8: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES, 16: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED, 268435455: rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL}
class AromaticityModel(Boost.Python.enum):
    AROMATICITY_CUSTOM: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM
    AROMATICITY_DEFAULT: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT
    AROMATICITY_MDL: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL
    AROMATICITY_MMFF94: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94
    AROMATICITY_RDKIT: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT
    AROMATICITY_SIMPLE: typing.ClassVar[AromaticityModel]  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'AROMATICITY_DEFAULT': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT, 'AROMATICITY_RDKIT': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT, 'AROMATICITY_SIMPLE': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE, 'AROMATICITY_MDL': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL, 'AROMATICITY_MMFF94': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94, 'AROMATICITY_CUSTOM': rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT, 1: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT, 2: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE, 4: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL, 8: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94, 268435455: rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM}
class BondWedgingParameters(Boost.Python.instance):
    """
    Parameters controlling how bond wedging is done.
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def wedgeTwoBondsIfPossible(self) -> bool:
        """If this is enabled then two bonds will be wedged at chiral
          centers subject to the following constraints:
            1. ring bonds will not be wedged
            2. bonds to chiral centers will not be wedged
            3. bonds separated by more than 120 degrees will not be
                wedged (default: False)"""
    @wedgeTwoBondsIfPossible.setter
    def wedgeTwoBondsIfPossible(self, value: bool) -> None: ...
class BoolVector(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__bit_iterator<std::__1::vector<bool, std::__1::allocator<bool>>, false, 0ul>> __iter__(boost::python::back_reference<std::__1::vector<bool, std::__1::allocator<bool>>&>)
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
class MolzipLabel(Boost.Python.enum):
    AtomMapNumber: typing.ClassVar[MolzipLabel]  # value = rdkit.Chem.rdmolops.MolzipLabel.AtomMapNumber
    AtomType: typing.ClassVar[MolzipLabel]  # value = rdkit.Chem.rdmolops.MolzipLabel.AtomType
    FragmentOnBonds: typing.ClassVar[MolzipLabel]  # value = rdkit.Chem.rdmolops.MolzipLabel.FragmentOnBonds
    Isotope: typing.ClassVar[MolzipLabel]  # value = rdkit.Chem.rdmolops.MolzipLabel.Isotope
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'AtomMapNumber': rdkit.Chem.rdmolops.MolzipLabel.AtomMapNumber, 'Isotope': rdkit.Chem.rdmolops.MolzipLabel.Isotope, 'FragmentOnBonds': rdkit.Chem.rdmolops.MolzipLabel.FragmentOnBonds, 'AtomType': rdkit.Chem.rdmolops.MolzipLabel.AtomType}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.MolzipLabel.AtomMapNumber, 1: rdkit.Chem.rdmolops.MolzipLabel.Isotope, 2: rdkit.Chem.rdmolops.MolzipLabel.FragmentOnBonds, 3: rdkit.Chem.rdmolops.MolzipLabel.AtomType}
class MolzipParams(Boost.Python.instance):
    """
    Parameters controlling how to zip molecules together
    
      OPTIONS:
          label : set the MolzipLabel option [default MolzipLabel.AtomMapNumber]
    
      MolzipLabel.AtomMapNumber: atom maps are on dummy atoms, zip together the corresponding
         attached atoms, i.e.  zip 'C[*:1]' 'N[*:1]' results in 'CN'
    
      MolzipLabel.Isotope: isotope labels are on dummy atoms, zip together the corresponding
         attached atoms, i.e.  zip 'C[1*]' 'N[1*]' results in 'CN'
    
      MolzipLabel.FragmentOnBonds: zip together molecules generated by fragment on bonds.
        Note the atom indices cannot change or be reordered from the output of fragmentOnBonds
    
      MolzipLabel.AtomTypes: choose the atom types to act as matching dummy atoms.
        i.e.  'C[V]' and 'N[Xe]' with atoms pairs [('V', 'Xe')] results in 'CN'
    """
    __instance_size__: typing.ClassVar[int] = 88
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    def setAtomSymbols(self, symbols: typing.Any) -> None:
        """
            Set the atom symbols used to zip mols together when using AtomType labeling
        
        """
    @property
    def alignCoordinates(self) -> bool:
        """if true and the input fragments have coordinates, the fragments
        will be aligned along connection vectors in the output molecule (default: False)"""
    @alignCoordinates.setter
    def alignCoordinates(self, value: bool) -> None: ...
    @property
    def enforceValenceRules(self) -> bool:
        """If true (default) enforce valences after zipping
        Setting this to false allows assembling chemically incorrect fragments. (default: True)"""
    @enforceValenceRules.setter
    def enforceValenceRules(self, value: bool) -> None: ...
    @property
    def generateCoordinates(self) -> bool:
        """If true will add depiction coordinates to input molecules and
        zipped molecule (for molzipFragments only) (default: False)"""
    @generateCoordinates.setter
    def generateCoordinates(self, value: bool) -> None: ...
    @property
    def label(*args, **kwargs):
        """
        Set the atom labeling system to zip together
        """
    @label.setter
    def label(*args, **kwargs):
        ...
class RemoveHsParameters(Boost.Python.instance):
    """
    Parameters controlling which Hs are removed.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def removeAndTrackIsotopes(self) -> bool:
        """hydrogens with non-default isotopes and store them in the _isotopicHs atom property such that AddHs() can add the same isotope at a later stage (default: False)"""
    @removeAndTrackIsotopes.setter
    def removeAndTrackIsotopes(self, value: bool) -> None: ...
    @property
    def removeDefiningBondStereo(self) -> bool:
        """hydrogens defining bond stereochemistry (default: False)"""
    @removeDefiningBondStereo.setter
    def removeDefiningBondStereo(self, value: bool) -> None: ...
    @property
    def removeDegreeZero(self) -> bool:
        """hydrogens that have no bonds (default: False)"""
    @removeDegreeZero.setter
    def removeDegreeZero(self, value: bool) -> None: ...
    @property
    def removeDummyNeighbors(self) -> bool:
        """hydrogens with at least one dummy-atom neighbor (default: False)"""
    @removeDummyNeighbors.setter
    def removeDummyNeighbors(self, value: bool) -> None: ...
    @property
    def removeHigherDegrees(self) -> bool:
        """hydrogens with two (or more) bonds (default: False)"""
    @removeHigherDegrees.setter
    def removeHigherDegrees(self, value: bool) -> None: ...
    @property
    def removeHydrides(self) -> bool:
        """hydrogens with formal charge -1 (default: False)"""
    @removeHydrides.setter
    def removeHydrides(self, value: bool) -> None: ...
    @property
    def removeInSGroups(self) -> bool:
        """hydrogens involved in SubstanceGroups (default: True)"""
    @removeInSGroups.setter
    def removeInSGroups(self, value: bool) -> None: ...
    @property
    def removeIsotopes(self) -> bool:
        """hydrogens with non-default isotopes (default: False)"""
    @removeIsotopes.setter
    def removeIsotopes(self, value: bool) -> None: ...
    @property
    def removeMapped(self) -> bool:
        """mapped hydrogens (default: True)"""
    @removeMapped.setter
    def removeMapped(self, value: bool) -> None: ...
    @property
    def removeNonimplicit(self) -> bool:
        """DEPRECATED (default: True)"""
    @removeNonimplicit.setter
    def removeNonimplicit(self, value: bool) -> None: ...
    @property
    def removeNontetrahedralNeighbors(self) -> bool:
        """hydrogens with neighbors that have non-tetrahedral stereochemistry (default: False)"""
    @removeNontetrahedralNeighbors.setter
    def removeNontetrahedralNeighbors(self, value: bool) -> None: ...
    @property
    def removeOnlyHNeighbors(self) -> bool:
        """hydrogens with bonds only to other hydrogens (default: False)"""
    @removeOnlyHNeighbors.setter
    def removeOnlyHNeighbors(self, value: bool) -> None: ...
    @property
    def removeWithQuery(self) -> bool:
        """hydrogens with queries defined (default: False)"""
    @removeWithQuery.setter
    def removeWithQuery(self, value: bool) -> None: ...
    @property
    def removeWithWedgedBond(self) -> bool:
        """hydrogens with wedged bonds to them (default: True)"""
    @removeWithWedgedBond.setter
    def removeWithWedgedBond(self, value: bool) -> None: ...
    @property
    def showWarnings(self) -> bool:
        """display warning messages for some classes of removed Hs (default: True)"""
    @showWarnings.setter
    def showWarnings(self, value: bool) -> None: ...
    @property
    def updateExplicitCount(self) -> bool:
        """DEPRECATED (default: False)"""
    @updateExplicitCount.setter
    def updateExplicitCount(self, value: bool) -> None: ...
class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUSTHS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
    SANITIZE_ALL: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL
    SANITIZE_CLEANUP: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
    SANITIZE_CLEANUPATROPISOMERS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS
    SANITIZE_CLEANUPCHIRALITY: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
    SANITIZE_CLEANUP_ORGANOMETALLICS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
    SANITIZE_FINDRADICALS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
    SANITIZE_KEKULIZE: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
    SANITIZE_NONE: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    SANITIZE_PROPERTIES: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
    SANITIZE_SETAROMATICITY: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
    SANITIZE_SETCONJUGATION: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
    SANITIZE_SETHYBRIDIZATION: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
    SANITIZE_SYMMRINGS: typing.ClassVar[SanitizeFlags]  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'SANITIZE_NONE': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, 'SANITIZE_CLEANUP': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP, 'SANITIZE_PROPERTIES': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES, 'SANITIZE_SYMMRINGS': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS, 'SANITIZE_KEKULIZE': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE, 'SANITIZE_FINDRADICALS': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS, 'SANITIZE_SETAROMATICITY': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY, 'SANITIZE_SETCONJUGATION': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION, 'SANITIZE_SETHYBRIDIZATION': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION, 'SANITIZE_CLEANUPCHIRALITY': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY, 'SANITIZE_CLEANUPATROPISOMERS': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS, 'SANITIZE_ADJUSTHS': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS, 'SANITIZE_CLEANUP_ORGANOMETALLICS': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS, 'SANITIZE_ALL': rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE, 1: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP, 2: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES, 4: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS, 8: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE, 16: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS, 32: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY, 64: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION, 128: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION, 256: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY, 2048: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS, 512: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS, 1024: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS, 268435455: rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL}
class StereoBondThresholds(Boost.Python.instance):
    """
    Constants used to set the thresholds for which single bonds can be made wavy.
    """
    CHIRAL_ATOM: typing.ClassVar[int] = 100000
    DBL_BOND_NO_STEREO: typing.ClassVar[int] = 1000
    DBL_BOND_SPECIFIED_STEREO: typing.ClassVar[int] = 10000
    DIRECTION_SET: typing.ClassVar[int] = 1000000
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
class StereoGroupAbsOptions(Boost.Python.enum):
    AlwaysInclude: typing.ClassVar[StereoGroupAbsOptions]  # value = rdkit.Chem.rdmolops.StereoGroupAbsOptions.AlwaysInclude
    NeverInclude: typing.ClassVar[StereoGroupAbsOptions]  # value = rdkit.Chem.rdmolops.StereoGroupAbsOptions.NeverInclude
    OnlyIncludeWhenOtherGroupsExist: typing.ClassVar[StereoGroupAbsOptions]  # value = rdkit.Chem.rdmolops.StereoGroupAbsOptions.OnlyIncludeWhenOtherGroupsExist
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'OnlyIncludeWhenOtherGroupsExist': rdkit.Chem.rdmolops.StereoGroupAbsOptions.OnlyIncludeWhenOtherGroupsExist, 'NeverInclude': rdkit.Chem.rdmolops.StereoGroupAbsOptions.NeverInclude, 'AlwaysInclude': rdkit.Chem.rdmolops.StereoGroupAbsOptions.AlwaysInclude}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.StereoGroupAbsOptions.OnlyIncludeWhenOtherGroupsExist, 1: rdkit.Chem.rdmolops.StereoGroupAbsOptions.NeverInclude, 2: rdkit.Chem.rdmolops.StereoGroupAbsOptions.AlwaysInclude}
class SubsetInfo(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 136
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def atomMapping(*args, **kwargs):
        """
        mapping from the original atom index to the subset atom index
        """
    @atomMapping.setter
    def atomMapping(*args, **kwargs):
        ...
    @property
    def bondMapping(*args, **kwargs):
        """
         mapping from the original bond index to the subset bond index
        """
    @bondMapping.setter
    def bondMapping(*args, **kwargs):
        ...
class SubsetMethod(Boost.Python.enum):
    BONDS: typing.ClassVar[SubsetMethod]  # value = rdkit.Chem.rdmolops.SubsetMethod.BONDS
    BONDS_BETWEEN_ATOMS: typing.ClassVar[SubsetMethod]  # value = rdkit.Chem.rdmolops.SubsetMethod.BONDS_BETWEEN_ATOMS
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'BONDS_BETWEEN_ATOMS': rdkit.Chem.rdmolops.SubsetMethod.BONDS_BETWEEN_ATOMS, 'BONDS': rdkit.Chem.rdmolops.SubsetMethod.BONDS}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdmolops.SubsetMethod.BONDS_BETWEEN_ATOMS, 1: rdkit.Chem.rdmolops.SubsetMethod.BONDS}
class SubsetOptions(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def clearComputedProps(self) -> bool:
        """clear all computed props on the subsetted molecule (default: False)"""
    @clearComputedProps.setter
    def clearComputedProps(self, value: bool) -> None: ...
    @property
    def conformerIdx(self) -> int:
        """What conformer idx to use for the coordinates default is -1 (default: 4294967295)"""
    @conformerIdx.setter
    def conformerIdx(self, value: int) -> None: ...
    @property
    def copyAsQuery(self) -> bool:
        """Return the subset as a query (default: False)"""
    @copyAsQuery.setter
    def copyAsQuery(self, value: bool) -> None: ...
    @property
    def copyCoordinates(self) -> bool:
        """Copy the active coordinates from the molecule (default: True)"""
    @copyCoordinates.setter
    def copyCoordinates(self, value: bool) -> None: ...
    @property
    def method(*args, **kwargs):
        """
        Subsetting method to use
        """
    @method.setter
    def method(*args, **kwargs):
        ...
    @property
    def sanitize(self) -> bool:
        """Sanitize the resulting subset (default: False)"""
    @sanitize.setter
    def sanitize(self, value: bool) -> None: ...
class UIntUIntMap(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__map_iterator<std::__1::__tree_iterator<std::__1::__value_type<unsigned int, unsigned int>, std::__1::__tree_node<std::__1::__value_type<unsigned int, unsigned int>, void*>*, long>>> __iter__(boost::python::back_reference<std::__1::map<unsigned int, unsigned int, std::__1::less<unsigned int>, std::__1::allocator<std::__1::pair<unsigned int const, unsigned int>>>&>)
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _vectN5RDKit9Chirality10StereoInfoE(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::Chirality::StereoInfo*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::Chirality::StereoInfo, std::__1::allocator<RDKit::Chirality::StereoInfo>>&>)
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
class map_indexing_suite_UIntUIntMap_entry(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __repr__(arg1: map_indexing_suite_UIntUIntMap_entry) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def data(self) -> int:
        """
        """
    def key(self) -> int:
        """
        """
@typing.overload
def AddHs(mol: Mol, params: AddHsParameters, onlyOnAtoms: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Adds hydrogens to the graph of a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - params: AddHsParameters object controlling the addition.
        
            - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
              considered to have Hs added to them
        
          RETURNS: a new molecule with added Hs
        
          NOTES:
        
            - The original molecule is *not* modified.
        
            - Much of the code assumes that Hs are not included in the molecular
              topology, so be *very* careful with the molecule that comes back from
              this function.\\n
    
    """
@typing.overload
def AddHs(mol: Mol, explicitOnly: bool = False, addCoords: bool = False, onlyOnAtoms: typing.Any = None, addResidueInfo: bool = False) -> rdkit.Chem.Mol:
    """
        Adds hydrogens to the graph of a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - explicitOnly: (optional) if this toggle is set, only explicit Hs will
              be added to the molecule.  Default value is 0 (add implicit and explicit Hs).
        
            - addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates
              set.  Default value is 0 (no 3D coords).
        
            - onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
              considered to have Hs added to them
        
            - addResidueInfo: (optional) if this is true, add residue info to
              hydrogen atoms (useful for PDB files).
        
          RETURNS: a new molecule with added Hs
        
          NOTES:
        
            - The original molecule is *not* modified.
        
            - Much of the code assumes that Hs are not included in the molecular
              topology, so be *very* careful with the molecule that comes back from
              this function.
        
        
    
    """
def AddRecursiveQuery(mol: Mol, query: Mol, atomIdx: int, preserveExistingQuery: bool = True) -> None:
    """
        Adds a recursive query to an atom
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as the recursive query (this will be copied)
        
            - atomIdx: the atom to modify
        
            - preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved
        
          RETURNS: None
        
        
    
    """
def AddStereoAnnotations(*args, **kwargs) -> None:
    """
        add R/S, relative stereo, and E/Z annotations to atoms and bonds
        
          Arguments:
           - mol: molecule to modify
           - absLabel: label for atoms in an ABS stereo group
           - orLabel: label for atoms in an OR stereo group
           - andLabel: label for atoms in an AND stereo group
           - cipLabel: label for chiral atoms that aren't in a stereo group.
           - bondLabel: label for CIP stereochemistry on bonds
        
         If any label is empty, the corresponding annotations will not be added.
        
         The labels can contain the following placeholders:
           - {id} - the stereo group's index
           - {cip} - the atom or bond's CIP stereochemistry
        
         Note that CIP labels will only be added if CIP stereochemistry has been
         assigned to the molecule.
        
    
    """
def AddWavyBondsForStereoAny(mol: Mol, clearDoubleBondFlags: bool = True, addWhenImpossible: int = 1000) -> None:
    """
        set wavy bonds around double bonds with STEREOANY stereo
          ARGUMENTS :
            - molecule : the molecule to update\\n -
            - conformer : the conformer to use to determine wedge direction
        
    
    """
def AdjustQueryProperties(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Returns a new molecule where the query properties of atoms have been modified.
    
    """
def AdjustQueryPropertiesWithGenericGroups(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Returns a new molecule where the query properties of atoms have been modified and generic group queries have been prepared.
    
    """
def AssignAtomChiralTagsFromMolParity(mol: Mol, replaceExistingTags: bool = True) -> None:
    """
        Sets the chiral tags on a molecule's atoms based on
          the molParity atom property.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - replaceExistingTags: if True, existing stereochemistry information will be cleared
            before running the calculation.
        
        
    
    """
def AssignAtomChiralTagsFromStructure(mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
        Sets the chiral tags on a molecule's atoms based on
          a 3D conformation.
          NOTE that this does not check to see if atoms are chiral centers (i.e. all
          substituents are different), it merely sets the chiral type flags based on the
          coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you
          want chiral flags only on actual stereocenters.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: the conformer id to use, -1 for the default 
            - replaceExistingTags: if True, existing stereochemistry information will be cleared
            before running the calculation.
        
        
    
    """
def AssignChiralTypesFromBondDirs(mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
        Uses bond directions to assign ChiralTypes to a molecule's atoms.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: (optional) the conformation to use 
            - replaceExistingTags: (optional) replace any existing information about stereochemistry
        
        
    
    """
def AssignRadicals(mol: Mol) -> None:
    """
        Assigns radical counts to atoms
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        
    
    """
def AssignStereochemistry(mol: Mol, cleanIt: bool = False, force: bool = False, flagPossibleStereoCenters: bool = False) -> None:
    """
        Assign stereochemistry tags to atoms and bonds.
          If useLegacyStereoPerception is true, it also does the CIP stereochemistry
          assignment for the molecule's atoms (R/S) and double bonds (Z/E).
          This assignment is based on legacy code which is fast, but is
          known to incorrectly assign CIP labels in some cases.
          instead, to assign CIP labels based on an accurate, though slower,
          implementation of the CIP rules, call CIPLabeler::assignCIPLabels().
          Chiral atoms will have a property '_CIPCode' indicating their chiral code.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - cleanIt: (optional) if provided, any existing values of the property `_CIPCode`
                will be cleared, atoms with a chiral specifier that aren't
              actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
              etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with 
              STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP 
              atom ranks will be marked STEREONONE. 
            - force: (optional) causes the calculation to be repeated, even if it has already
              been done
            - flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
              atoms that are possible stereocenters
        
    
    """
def AssignStereochemistryFrom3D(mol: Mol, confId: int = -1, replaceExistingTags: bool = True) -> None:
    """
        Uses a conformer (should be 3D) to assign ChiralTypes to a molecule's atoms
                and stereo flags to its bonds
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - confId: (optional) the conformation to use 
            - replaceExistingTags: (optional) replace any existing information about stereochemistry
        
        
    
    """
def AtomHasConjugatedBond(atom: Atom) -> bool:
    """
        returns whether or not the atom is involved in a conjugated bond
    
    """
def CanonicalizeStereoGroups(mol: Mol, outputAbsoluteGroups: StereoGroupAbsOptions = ..., maxStereoGroups: int = 12) -> rdkit.Chem.Mol:
    """
        Rationalize Enhanced Stereo indications to a canonical form 
                
                  ARGUMENTS:
                
                    - molecule: the molecule to update
                    -StereoGroupAbsOptions outputAbsoluteGroups: controls output of abs groups: 
                      one of: OnlyIncludeWhenOtherGroupsExist, NeverInclude, AlwaysInclude 
                     maxStereoGroups: maximm number of OR or AND stereo groups to process (default is 12): 
                
                
         
    
    """
def Cleanup(mol: Mol) -> None:
    """
        cleans up certain common bad functionalities in the molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
    
    """
def CleanupAtropisomers(mol: Mol) -> None:
    """
        removes bogus atropisomeric markers (e.g. those without sp2 begin and end atoms)
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
    
    """
def CleanupChirality(mol: Mol) -> None:
    """
        removes bogus chirality markers (e.g. tetrahedral flags on non-sp3 centers)
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
    
    """
def CleanupOrganometallics(mol: Mol) -> None:
    """
        cleans up certain common bad functionalities in the organometallic molecule
        
          Note that this function is experimental and may either change in behavior
          or be replaced with something else in future releases.
        
                ARGUMENTS :
        
         - mol : the molecule to use
        
         NOTES :
        
         - The molecule is modified in place.
        
         
    
    """
def CleanupStereoGroups(mol: Mol) -> None:
    """
        removes atoms without specified chirality from stereo groups
    
    """
def CollapseAttachmentPoints(mol: Mol, markedOnly: bool = True) -> None:
    """
        dummy atoms in the graph are removed and replaced with attachment point annotations on the attached atoms
        
          Arguments:
           - mol: molecule to be modified
           - markedOnly: if true, only dummy atoms with the _fromAttachPoint
             property will be collapsed
        
          In order for a dummy atom to be considered for collapsing it must have:
           - degree 1 with a single or unspecified bond
           - the bond to it can not be wedged
           - either no query or be an AtomNullQuery
        
    
    """
def CombineMols(mol1: Mol, mol2: Mol, offset: Point3D = ...) -> rdkit.Chem.Mol:
    """
        Combine the atoms from two molecules to produce a third
    
    """
def ComputeAtomCIPRanks(mol: Mol) -> tuple:
    """
        Computes the CIP ranks for the atoms in a molecule.
          The ranks are stored as an atom property '_CIPRank' and returned as a tuple.
        
    
    """
def ConvertGenericQueriesToSubstanceGroups(mol: Mol) -> None:
    """
        documentation
    
    """
@typing.overload
def CopyMolSubset(mol: Mol, atomIndices: typing.Any, bondIndices: typing.Any, options: SubsetOptions = ...) -> rdkit.Chem.Mol:
    """
        Extract a subgraph from an ROMol. Bonds, atoms, substance groups and 
        stereo groups are only extracted to the subgraph if all participant entities 
        are contained within the given atoms and bonds. 
        
        ARGUMENTS:
         - mol - starting mol 
         - atoms - indices atoms to extract 
         - bonds - indices bonds to extract 
         - subsetInfo - optional subsetInfo to record the atoms and bonds used 
         - options - optional subset options, note the method is ignored since all the atoms and bonds are sp 
        
        
    
    """
@typing.overload
def CopyMolSubset(mol: Mol, atomIndices: typing.Any, bondIndices: typing.Any, subsetInfo: SubsetInfo, options: SubsetOptions = ...) -> rdkit.Chem.Mol:
    """
        Extract a subgraph from an ROMol. Bonds, atoms, substance groups and 
        stereo groups are only extracted to the subgraph if all participant entities 
        are contained within the given atoms and bonds. 
        
        ARGUMENTS:
         - mol - starting mol 
         - atoms - indices atoms to extract 
         - bonds - indices bonds to extract 
         - subsetInfo - optional subsetInfo to record the atoms and bonds used 
         - options - optional subset options, note the method is ignored since all the atoms and bonds are sp 
        
        
    
    """
@typing.overload
def CopyMolSubset(mol: Mol, path: typing.Any, options: SubsetOptions = ...) -> rdkit.Chem.Mol:
    """
        copy a subset of a molecule
    
    """
@typing.overload
def CopyMolSubset(mol: Mol, path: typing.Any, subsetInfo: SubsetInfo, options: SubsetOptions = ...) -> rdkit.Chem.Mol:
    """
        copy a subset of a molecule
    
    """
def CountAtomElec(atom: Atom) -> int:
    """
        returns the number of electrons available on an atom to donate for aromaticity
    
    """
def DativeBondsToHaptic(mol: Mol) -> rdkit.Chem.Mol:
    """
        Does the reverse of hapticBondsToDative.  If there are multiple
        contiguous atoms attached by dative bonds to an atom (probably a metal
        atom), the dative bonds will be replaced by a dummy atom in their
        centre attached to the (metal) atom by a dative bond, which is
        labelled with ENDPTS of the atoms that had the original dative bonds.
        
        ARGUMENTS:
        
          - mol: the molecule to use
        
        RETURNS:
          a modified copy of the molecule
    
    """
def DeleteSubstructs(mol: Mol, query: Mol, onlyFrags: bool = False, useChirality: bool = False) -> rdkit.Chem.Mol:
    """
        Removes atoms matching a substructure query from a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as a substructure query
        
            - onlyFrags: (optional) if this toggle is set, atoms will only be removed if
              the entire fragment in which they are found is matched by the query.
              See below for examples.
              Default value is 0 (remove the atoms whether or not the entire fragment matches)
        
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a new molecule with the substructure removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - DeleteSubstructs('CCOC','OC') -> 'CC'
        
            - DeleteSubstructs('CCOC','OC',1) -> 'CCOC'
        
            - DeleteSubstructs('CCOCCl.Cl','Cl',1) -> 'CCOCCl'
        
            - DeleteSubstructs('CCOCCl.Cl','Cl') -> 'CCOC'
        
        
    
    """
def DetectBondStereoChemistry(mol: Mol, conformer: Conformer) -> None:
    """
        Assign stereochemistry to bonds based on coordinates and a conformer.
                DEPRECATED
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
            - conformer: Conformer providing the coordinates
        
        
    
    """
def DetectBondStereochemistry(mol: Mol, confId: int = -1) -> None:
    """
        DEPRECATED
            - mol: the molecule to be modified
            - confId: Conformer to use for the coordinates
        
        
    
    """
def DetectChemistryProblems(mol: Mol, sanitizeOps: int = ...) -> tuple:
    """
        checks for chemistry problems
    
    """
def ExpandAttachmentPoints(mol: Mol, addAsQueries: bool = True, addCoords: bool = True) -> None:
    """
        attachment points encoded as attachPt properties are added to the graph as dummy atoms
        
          Arguments:
           - mol: molecule to be modified
           - addAsQueries: if true, the dummy atoms will be added as null queries
                (i.e. they will match any atom in a substructure search)
           - addCoords: if true and the molecule has one or more conformers, 
                positions for the attachment points will be added to the conformer(s)
        
    
    """
def FastFindRings(mol: Mol) -> None:
    """
        Does a non-SSSR ring finding for a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use.
        
          RETURNS: Nothing
        
        
    
    """
def FindAllPathsOfLengthN(mol: Mol, length: int, useBonds: bool = True, useHs: bool = False, rootedAtAtom: int = -1, onlyShortestPaths: bool = False) -> _listNSt3__16vectorIiNS_9allocatorIiEEEE:
    """
        Finds all paths of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target length for the paths.
        
            - useBonds: (optional) toggles the use of bond indices in the paths.
              Otherwise atom indices are used.  *Note* this behavior is different
              from that for subgraphs.
              Defaults to 1.
        
            - rootedAtAtom: (optional) if nonzero, only paths from the specified
              atom will be returned.
        
            - onlyShortestPaths: (optional) if set then only paths which are <= the shortest
              path between the begin and end atoms will be included in the results
        
          RETURNS: a tuple of tuples with IDs for the bonds.
        
          NOTES: 
        
           - Difference between _subgraphs_ and _paths_ :: 
        
               Subgraphs are potentially branched, whereas paths (in our 
               terminology at least) cannot be.  So, the following graph: 
        
                    C--0--C--1--C--3--C
                          |
                          2
                          |
                          C
        
               has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
               but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
        
        
    
    """
def FindAllSubgraphsOfLengthMToN(mol: Mol, min: int, max: int, useHs: bool = False, rootedAtAtom: int = -1) -> typing.Any:
    """
        Finds all subgraphs of a particular length in a molecule
          See documentation for FindAllSubgraphsOfLengthN for definitions
        
        
    
    """
def FindAllSubgraphsOfLengthN(mol: Mol, length: int, useHs: bool = False, rootedAtAtom: int = -1) -> _listNSt3__16vectorIiNS_9allocatorIiEEEE:
    """
        Finds all subgraphs of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target number of bonds for the subgraphs.
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
              atom will be returned.
        
          RETURNS: a tuple of 2-tuples with bond IDs
        
          NOTES: 
        
           - Difference between _subgraphs_ and _paths_ :: 
        
               Subgraphs are potentially branched, whereas paths (in our 
               terminology at least) cannot be.  So, the following graph: 
        
                    C--0--C--1--C--3--C
                          |
                          2
                          |
                          C
          has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
          but only 2 _paths_ of length 3: (0,1,3),(2,1,3)
        
        
    
    """
def FindAtomEnvironmentOfRadiusN(mol: Mol, radius: int, rootedAtAtom: int, useHs: bool = False, enforceSize: bool = True, atomMap: typing.Any = None) -> typing.Sequence[int]:
    """
        Find bonds of a particular radius around an atom. 
                 Return empty result if there is no bond at the requested radius.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - radius: an integer with the target radius for the environment.
        
            - rootedAtAtom: the atom to consider
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - enforceSize (optional) If set to False, all bonds within the requested radius is 
              collected. Defaults to 1. 
        
            - atomMap: (optional) If provided, it will measure the minimum distance of the atom 
              from the rooted atom (start with 0 from the rooted atom). The result is a pair of 
              the atom ID and the distance. 
        
          RETURNS: a vector of bond IDs
        
        
    
    """
def FindMesoCenters(mol: Mol, includeIsotopes: bool = True, includeAtomMaps: bool = False) -> typing.Any:
    """
        returns the meso centers in a molecule (if any).
            
          ARGUMENTS:
            
            - mol: the molecule to use
            - includeIsotopes: (optional) toggles whether or not isotopes should be included in the
              calculation.
            - includeAtomMaps: (optional) toggles whether or not atom maps should be included in the
              calculation.
            
    
    """
def FindPotentialStereo(mol: Mol, cleanIt: bool = False, flagPossible: bool = True) -> typing.Sequence[rdkit.Chem.StereoInfo]:
    """
        find potential stereo elements in a molecule and returns them as StereoInfo objects
        Note that this function is still somewhat experimental and the API
        and results may change in a future release.
    
    """
def FindPotentialStereoBonds(mol: Mol, cleanIt: bool = False) -> None:
    """
        Find bonds than can be cis/trans in a molecule and mark them as 'any'.
                 This function finds any double bonds that can potentially be part
                 of a cis/trans system. No attempt is made here to mark them cis or trans
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - cleanIt: (optional) if this option is set to true, any previous marking of _CIPCode
                       on the bond is cleared - otherwise it is left untouched
        
        
    
    """
def FindRingFamilies(mol: Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False) -> None:
    """
        Generate Unique Ring Families.
        
          ARGUMENTS:
        
            - mol: the molecule to use.
            - includeDativeBonds: whether or not dative bonds should be included in the ring families finding.
            - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring families .
        
          RETURNS: Nothing
        
        
    
    """
def FindUniqueSubgraphsOfLengthN(mol: Mol, length: int, useHs: bool = False, useBO: bool = True, rootedAtAtom: int = -1) -> _listNSt3__16vectorIiNS_9allocatorIiEEEE:
    """
        Finds unique subgraphs of a particular length in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - length: an integer with the target number of bonds for the subgraphs.
        
            - useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
              should be included in the results.
              Defaults to 0.
        
            - useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from
              another.
              Defaults to 1.
        
            - rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
              atom will be returned.
        
          RETURNS: a tuple of tuples with bond IDs
        
        
        
    
    """
def FragmentOnBRICSBonds(mol: Mol) -> rdkit.Chem.Mol:
    """
        Return a new molecule with all BRICS bonds broken
    
    """
def FragmentOnBonds(mol: Mol, bondIndices: typing.Any, addDummies: bool = True, dummyLabels: typing.Any = None, bondTypes: typing.Any = None, cutsPerAtom: list = []) -> rdkit.Chem.Mol:
    """
        Return a new molecule with all specified bonds broken
        
          ARGUMENTS:
        
              - mol            - the molecule to be modified
              - bondIndices    - indices of the bonds to be broken
              - addDummies  - toggles addition of dummy atoms to indicate where 
                bonds were broken
              - dummyLabels - used to provide the labels to be used for the dummies.
                the first element in each pair is the label for the dummy
                that replaces the bond's beginAtom, the second is for the 
                dummy that replaces the bond's endAtom. If not provided, the
                dummies are labeled with atom indices.
              - bondTypes - used to provide the bond type to use between the
                fragments and the dummy atoms. If not provided, defaults to single. 
              - cutsPerAtom - used to return the number of cuts made at each atom. 
        
          RETURNS:
              a new Mol with the modifications
        
    
    """
def FragmentOnSomeBonds(mol: Mol, bondIndices: typing.Any, numToBreak: int = 1, addDummies: bool = True, dummyLabels: typing.Any = None, bondTypes: typing.Any = None, returnCutsPerAtom: bool = False) -> tuple:
    """
        fragment on some bonds
    
    """
def Get3DDistanceMatrix(mol: Mol, confId: int = -1, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> typing.Any:
    """
        Returns the molecule's 3D distance matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - confId: (optional) chooses the conformer Id to use
              Default value is -1.
        
            - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
              matrix (to return a "Balaban" distance matrix).
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats with the distance matrix
        
        
    
    """
def GetAdjacencyMatrix(mol: Mol, useBO: bool = False, emptyVal: int = 0, force: bool = False, prefix: str = '') -> typing.Any:
    """
        Returns the molecule's adjacency matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - useBO: (optional) toggles use of bond orders in calculating the matrix.
              Default value is 0.
        
            - emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats containing the adjacency matrix
        
        
    
    """
def GetAllowNontetrahedralChirality() -> bool:
    """
        returns whether or not recognition of non-tetrahedral chirality from 3D structures is enabled
    
    """
def GetDistanceMatrix(mol: Mol, useBO: bool = False, useAtomWts: bool = False, force: bool = False, prefix: str = '') -> typing.Any:
    """
        Returns the molecule's topological distance matrix.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - useBO: (optional) toggles use of bond orders in calculating the distance matrix.
              Default value is 0.
        
            - useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
              matrix (to return a "Balaban" distance matrix).
              Default value is 0.
        
            - force: (optional) forces the calculation to proceed, even if there is a cached value.
              Default value is 0.
        
            - prefix: (optional, internal use) sets the prefix used in the property cache
              Default value is .
        
          RETURNS: a Numeric array of floats with the distance matrix
        
        
    
    """
def GetFormalCharge(mol: Mol) -> int:
    """
        Returns the formal charge for the molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
        
    
    """
def GetMolFrags(mol: Mol, asMols: bool = False, sanitizeFrags: bool = True, frags: typing.Any = None, fragsMolAtomMapping: typing.Any = None) -> tuple:
    """
        Finds the disconnected fragments from a molecule.
        
          For example, for the molecule 'CC(=O)[O-].[NH3+]C' GetMolFrags() returns
          ((0, 1, 2, 3), (4, 5))
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - asMols: (optional) if this is provided and true, the fragments
              will be returned as molecules instead of atom ids.
            - sanitizeFrags: (optional) if this is provided and true, the fragments
              molecules will be sanitized before returning them.
            - frags: (optional, defaults to None) if asMols is true and this is provided
               as an empty list, the result will be mol.GetNumAtoms() long on return and
               will contain the fragment assignment for each Atom
            - fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this
              is provided as an empty list, the result will be numFrags long on 
              return, and each entry will contain the indices of the Atoms in that fragment:
              [(0, 1, 2, 3), (4, 5)]
        
          RETURNS: a tuple of tuples with IDs for the atoms in each fragment
                   or a tuple of molecules.
        
        
    
    """
def GetMostSubstitutedCoreMatch(mol: Mol, core: Mol, matches: typing.Any) -> typing.Any:
    """
        Postprocesses the results of a mol.GetSubstructMatches(core) call 
        where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
        It returns the match with the largest number of non-hydrogen matches to 
        the terminal dummy atoms.
        
          ARGUMENTS:
        
            - mol: the molecule GetSubstructMatches was run on
        
            - core: the molecule used as a substructure query
        
            - matches: the result returned by GetSubstructMatches
        
          RETURNS: the tuple where terminal dummy atoms in the core match the largest 
                   number of non-hydrogen atoms in mol
        
    
    """
def GetSSSR(mol: Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False) -> _vectNSt3__16vectorIiNS_9allocatorIiEEEE:
    """
        Get the smallest set of simple rings for a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to use.
            - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
            - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.
        
          RETURNS: a sequence of sequences containing the rings found as atom ids
                 The length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.
        
        
    
    """
def GetShortestPath(mol: Mol, aid1: int, aid2: int) -> tuple:
    """
        Find the shortest path between two atoms using the Bellman-Ford algorithm.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - idx1: index of the first atom
            - idx2: index of the second atom
        
        
    
    """
def GetSymmSSSR(mol: Mol, includeDativeBonds: bool = False, includeHydrogenBonds: bool = False) -> _vectNSt3__16vectorIiNS_9allocatorIiEEEE:
    """
        Get a symmetrized SSSR for a molecule.
        
          The symmetrized SSSR is at least as large as the SSSR for a molecule.
          In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be
          a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).
        
          ARGUMENTS:
        
            - mol: the molecule to use.
            - includeDativeBonds: whether or not dative bonds should be included in the ring finding.
            - includeHydrogenBonds: whether or not hydrogen bonds should be included in the ring finding.
        
          RETURNS: a sequence of sequences containing the rings found as atom ids
        
        
    
    """
def GetUseLegacyStereoPerception() -> bool:
    """
        returns whether or not the legacy stereo perception code is being used
    
    """
def HapticBondsToDative(mol: Mol) -> rdkit.Chem.Mol:
    """
        One way of showing haptic bonds (such as cyclopentadiene to
        iron in ferrocene) is to use a dummy atom with a dative bond to the
        iron atom with the bond labelled with the atoms involved in the
        organic end of the bond.  Another way is to have explicit dative
        bonds from the atoms of the haptic group to the metal atom.  This
        function converts the former representation to the latter.
        
        ARGUMENTS:
        
          - mol: the molecule to use
        
        RETURNS:
          a modified copy of the molecule
    
    """
def HasQueryHs(mol: Mol) -> tuple:
    """
        Check to see if the molecule has query Hs, this is normally used on query molecules
        such as those returned from MolFromSmarts
        Example: 
              (hasQueryHs, hasUnmergeableQueryHs) = HasQueryHs(mol)
        
        if hasUnmergeableQueryHs, these query hs cannot be removed by calling
        MergeQueryHs
    
    """
def Kekulize(mol: Mol, clearAromaticFlags: bool = False, canonical: bool = True) -> None:
    """
        Kekulizes the molecule
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
              molecule will be marked non-aromatic following the kekulization.
              Default value is False.
        
            - canonical: (optional) if true, uses canonical atom ranking so
              that the kekulization result is independent of the atom ordering in the
              molecule.  Set to false to skip the ranking step for better performance
              when deterministic output is not required (e.g. during sanitization).
              Note, this "canonical" order only really makes sense when the molecule's
              chemistry is sane, like after sanitization. If stereochemistry hasn't been
              perceived, the chemistry of the molecule is inconsistent, and
              "canonical" atom ranks are only a technical artifact.
        
          NOTES:
        
            - The molecule is modified in place.
        
            - this does not modify query bonds which have bond type queries (like those
              which come from SMARTS) or rings containing them.
        
            - even if clearAromaticFlags is False the BondType for all modified
              aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
              Kekulization.
        
        
    
    """
def KekulizeIfPossible(mol: Mol, clearAromaticFlags: bool = False, canonical: bool = True) -> None:
    """
        Kekulizes the molecule if possible. Otherwise the molecule is not modified
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the 
              molecule will be marked non-aromatic if the kekulization succeds.
              Default value is False.
        
            - canonical: (optional) if true  uses canonical atom ranking so
              that the kekulization result is independent of the atom ordering in the
              molecule.  Set to false to skip the ranking step for better performance
              when deterministic output is not required (e.g. during sanitization).
              Note, this "canonical" order only really makes sense when the molecule's
              chemistry is sane, like after sanitization. If stereochemistry hasn't been
              perceived, the chemistry of the molecule is inconsistent, and
              "canonical" atom ranks are only a technical artifact.
        
        \\n\\
            - canonical: (optional) if True, uses canonical atom ranking so that the\\n\\
              kekulization result is independent of the atom ordering in the molecule.\\n\\
              Default value is False.\\n\\
        \\n\\
          NOTES:\\n\\
        \\n\\
            - The molecule is modified in place.\\n\\
            
    
    """
def LayeredFingerprint(mol: Mol, layerFlags: int = 4294967295, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, branchedPaths: bool = True, fromAtoms: typing.Any = 0) -> ExplicitBitVect:
    """
        Returns a layered fingerprint for a molecule
        
          NOTE: This function is experimental. The API or results may change from
            release to release.
        
          Explanation of the algorithm below.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - layerFlags: (optional) which layers to include in the fingerprint
              See below for definitions. Defaults to all.
        
            - minPath: (optional) minimum number of bonds to include in the subgraphs
              Defaults to 1.
        
            - maxPath: (optional) maximum number of bonds to include in the subgraphs
              Defaults to 7.
        
            - fpSize: (optional) number of bits in the fingerprint
              Defaults to 2048.
        
            - atomCounts: (optional) 
              if provided, this should be a list at least as long as the number of atoms
              in the molecule. It will be used to provide the count of the number 
              of paths that set bits each atom is involved in.
              NOTE: the list is not zeroed out here.
        
            - setOnlyBits: (optional) 
              if provided, only bits that are set in this bit vector will be set
              in the result. This is essentially the same as doing:
                   res &= setOnlyBits
              but also has an impact on the atomCounts (if being used)
        
            - branchedPaths: (optional) if set both branched and unbranched paths will be
              used in the fingerprint.
              Defaults to True.
        
            - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
              starting from these atoms will be used.
              Defaults to empty.
        
          RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits
        
          Layer definitions:
             - 0x01: pure topology
             - 0x02: bond order
             - 0x04: atom types
             - 0x08: presence of rings
             - 0x10: ring sizes
             - 0x20: aromaticity
        
        
        
    
    """
def MergeQueryHs(mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False) -> rdkit.Chem.Mol:
    """
        merges hydrogens into their neighboring atoms as queries
    
    """
def MolAddRecursiveQueries(mol: Mol, queries: dict, propName: str) -> None:
    """
        Adds named recursive queries to atoms
        
    
    """
def MurckoDecompose(mol: Mol) -> rdkit.Chem.Mol:
    """
        Do a Murcko decomposition and return the scaffold
    
    """
def NeedsHs(mol: Mol) -> bool:
    """
        returns whether or not the molecule needs to have Hs added
    
    """
def ParseMolQueryDefFile(fileobj: typing.Any, standardize: bool = True, delimiter: str = '\t', comment: str = '//', nameColumn: int = 0, smartsColumn: int = 1) -> dict:
    """
        reads query definitions from a simply formatted file
        
    
    """
def PathToSubmol(mol: Mol, path: typing.Any, useQuery: bool = False, atomMap: typing.Any = None) -> rdkit.Chem.Mol:
    """
    """
@typing.overload
def PatternFingerprint(mol: Mol, fpSize: int = 2048, atomCounts: list = [], setOnlyBits: ExplicitBitVect = None, tautomerFingerprints: bool = False) -> ExplicitBitVect:
    """
        A fingerprint using SMARTS patterns 
        
          NOTE: This function is experimental. The API or results may change from
            release to release.
        
    
    """
@typing.overload
def PatternFingerprint(mol: MolBundle, fpSize: int = 2048, setOnlyBits: ExplicitBitVect = None, tautomerFingerprints: bool = False) -> ExplicitBitVect:
    """
        A fingerprint using SMARTS patterns 
        
          NOTE: This function is experimental. The API or results may change from
            release to release.
        
    
    """
def RDKFingerprint(mol: Mol, minPath: int = 1, maxPath: int = 7, fpSize: int = 2048, nBitsPerHash: int = 2, useHs: bool = True, tgtDensity: float = 0.0, minSize: int = 128, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: typing.Any = 0, fromAtoms: typing.Any = 0, atomBits: typing.Any = None, bitInfo: typing.Any = None) -> ExplicitBitVect:
    """
        Returns an RDKit topological fingerprint for a molecule
        
          Explanation of the algorithm below.
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
            - minPath: (optional) minimum number of bonds to include in the subgraphs
              Defaults to 1.
        
            - maxPath: (optional) maximum number of bonds to include in the subgraphs
              Defaults to 7.
        
            - fpSize: (optional) number of bits in the fingerprint
              Defaults to 2048.
        
            - nBitsPerHash: (optional) number of bits to set per path
              Defaults to 2.
        
            - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
              has explicit Hs.
              Defaults to True.
        
            - tgtDensity: (optional) fold the fingerprint until this minimum density has
              been reached
              Defaults to 0.
        
            - minSize: (optional) the minimum size the fingerprint will be folded to when
              trying to reach tgtDensity
              Defaults to 128.
        
            - branchedPaths: (optional) if set both branched and unbranched paths will be
              used in the fingerprint.
              Defaults to True.
        
            - useBondOrder: (optional) if set both bond orders will be used in the path hashes
              Defaults to True.
        
            - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
              Defaults to empty.
        
            - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
              starting from these atoms will be used.
              Defaults to empty.
        
            - atomBits: (optional) an empty list. If provided, the result will contain a list 
              containing the bits each atom sets.
              Defaults to empty.
        
            - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
              with bits as keys and corresponding bond paths as values.
              Defaults to empty.
        
          RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits
        
          ALGORITHM:
        
           This algorithm functions by find all subgraphs between minPath and maxPath in
           length.  For each subgraph:
        
             1) A hash is calculated.
        
             2) The hash is used to seed a random-number generator
        
             3) _nBitsPerHash_ random numbers are generated and used to set the corresponding
                bits in the fingerprint
        
        
        
    
    """
def ReapplyMolBlockWedging(mol: Mol, allBondTypes: bool = True, verify: bool = False) -> None:
    """
        Set the wedging to that which was read from the original
             MolBlock, over-riding anything that was originally there.
        
                  ARGUMENTS:
                
                    - molecule: the molecule to update
                    - allBondTypes: reapply the wedging also on bonds other
                      than single and aromatic ones
                    - verify: if true, the function will check that the wedges are only
                      applied in sensible places (i.e.single bonds connected to chiral
                      centers or atropisomeric bonds)
                
                
        
    
    """
def RemoveAllHs(mol: Mol, sanitize: bool = True) -> rdkit.Chem.Mol:
    """
        Returns a copy of the molecule with all Hs removed.
    
    """
@typing.overload
def RemoveHs(mol: Mol, implicitOnly: bool = False, updateExplicitCount: bool = False, sanitize: bool = True) -> rdkit.Chem.Mol:
    """
        Removes any hydrogens from the graph of a molecule.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - implicitOnly: (optional) if this toggle is set, only implicit Hs will
              be removed from the graph.  Default value is 0 (remove implicit and explicit Hs).
        
            - updateExplicitCount: (optional) if this toggle is set, the explicit H count on atoms with 
              Hs will be updated. Default value is 0 (do not update explicit H count).
        
            - sanitize: (optional) if this toggle is set, the molecule will be sanitized after the Hs
              are removed. Default value is 1 (do sanitize).
        
          RETURNS: a new molecule with the Hs removed
        
          NOTES:
        
            - The original molecule is *not* modified.
            - Hydrogens which aren't connected to a heavy atom will not be
              removed.  This prevents molecules like [H][H] from having
              all atoms removed.
            - Labelled hydrogen (e.g. atoms with atomic number=1, but isotope > 1),
              will not be removed.
            - two coordinate Hs, like the central H in C[H-]C, will not be removed
            - Hs connected to dummy atoms will not be removed
            - Hs that are part of the definition of double bond Stereochemistry
              will not be removed
            - Hs that are not connected to anything else will not be removed
        
         
    
    """
@typing.overload
def RemoveHs(mol: Mol, params: RemoveHsParameters, sanitize: bool = True) -> rdkit.Chem.Mol:
    """
        Returns a copy of the molecule with Hs removed. Which Hs are removed is controlled by the params argument
    
    """
def RemoveNonExplicit3DChirality(mol: Mol) -> None:
    """
        Remove chiral markings that were derived from a 3D mol but were not 
                explicity marked in the mol block. (wedge bond or CFG indication
                
                  ARGUMENTS:
                
                    - molecule: the molecule to update
                
                
        
    
    """
def RemoveStereochemistry(mol: Mol) -> None:
    """
        Removes all stereochemistry info from the molecule.
        
        
    
    """
def RenumberAtoms(mol: Mol, newOrder: typing.Any) -> rdkit.Chem.Mol:
    """
        Returns a copy of a molecule with renumbered atoms
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - newOrder: the new ordering the atoms (should be numAtoms long)
              for example: if newOrder is [3,2,0,1], then atom 3 in the original 
              molecule will be atom 0 in the new one
        
        
        
    
    """
@typing.overload
def ReplaceCore(mol: Mol, core: Mol, matches: typing.Any, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False) -> rdkit.Chem.Mol:
    """
        Removes the core of a molecule and labels the sidechains with dummy atoms based on
        The matches indices given in the matching vector matches.
        Calling:
          ReplaceCore(mol,core,mol.GetSubstructMatch(core))
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - coreQuery: the molecule to be used as a substructure query for recognizing the core
        
            - matches: a matching vector of the type returned by mol.GetSubstructMatch(...)
        
            - replaceDummies: toggles replacement of atoms that match dummies in the query
        
            - labelByIndex: toggles labeling the attachment point dummy atoms with 
              the index of the core atom they're attached to.
        
            - requireDummyMatch: if the molecule has side chains that attach at points not
              flagged with a dummy, it will be rejected (None is returned)
        
          RETURNS: a new molecule with the core removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        EXAMPLES:
        
            >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, ReplaceCore
            >>> mol = MolFromSmiles('C1ONNCC1')
            >>> core = MolFromSmiles('NN')
        
            >>> MolToSmiles(ReplaceCore(mol, core, mol.GetSubstructMatch(core)))
            '[1*]OCCC[2*]'
        
            Since NN is symmetric, we should actually get two matches here if we don't
            uniquify the matches.
        
            >>> [MolToSmiles(ReplaceCore(mol, core, match))
            ...     for match in mol.GetSubstructMatches(core, uniquify=False)]
            ['[1*]OCCC[2*]', '[1*]CCCO[2*]']
        
        
    
    """
@typing.overload
def ReplaceCore(mol: Mol, coreQuery: Mol, replaceDummies: bool = True, labelByIndex: bool = False, requireDummyMatch: bool = False, useChirality: bool = False) -> rdkit.Chem.Mol:
    """
        Removes the core of a molecule and labels the sidechains with dummy atoms.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - coreQuery: the molecule to be used as a substructure query for recognizing the core
        
            - replaceDummies: toggles replacement of atoms that match dummies in the query
        
            - labelByIndex: toggles labeling the attachment point dummy atoms with 
              the index of the core atom they're attached to.
        
            - requireDummyMatch: if the molecule has side chains that attach at points not
              flagged with a dummy, it will be rejected (None is returned)
        
            - useChirality: use chirality matching in the coreQuery
        
          RETURNS: a new molecule with the core removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        
          EXAMPLES:
        
           >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFromSmarts, ReplaceCore
        
           Basic usage: remove a core as specified by SMILES (or another molecule).
           To get the atom labels which are stored as an isotope of the matched atom, 
           the output must be written as isomeric smiles.  
           A small confusion is that atom isotopes of 0 aren't shown in smiles strings.
        
           Here we remove a ring and leave the decoration (r-group) behind.
        
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCCC1CCC1'),MolFromSmiles('C1CCC1')),
           ...             isomericSmiles=True)
           '[1*]CCC'
        
           The isotope label by default is matched by the first connection found. In order to
           indicate which atom the decoration is attached in the core query, use labelByIndex=True.
           Here the attachment is from the third atom in the smiles string, which is indexed by 3
           in the core, like all good computer scientists expect, atoms indices start at 0.
        
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCN1CCC1'),MolFromSmiles('C1CCN1'),
           ...                         labelByIndex=True),
           ...   isomericSmiles=True)
           '[3*]CC'
        
           Non-core matches just return None
        
           >>> ReplaceCore(MolFromSmiles('CCC1CC1'),MolFromSmiles('C1CCC1'))
        
           The bond between atoms are considered part of the core and are removed as well
        
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CC2C1CCC2'),MolFromSmiles('C1CCC1')),
           ...             isomericSmiles=True)
           '[1*]CCC[2*]'
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmiles('N')),
           ...             isomericSmiles=True)
           '[1*]CCCC[2*]'
        
           When using dummy atoms, cores should be read in as SMARTS.  When read as SMILES
           dummy atoms only match other dummy atoms.
           The replaceDummies flag indicates whether matches to the dummy atoms should be considered as part
           of the core or as part of the decoration (r-group)
        
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
           ...                         replaceDummies=True),
           ...             isomericSmiles=True)
           '[1*]CC[2*]'
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
           ...                         replaceDummies=False),
           ...             isomericSmiles=True)
           '[1*]CCCC[2*]'
        
        
           >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CCC1CN'),MolFromSmarts('C1CCC1[*]'),
           ...                         replaceDummies=False),
           ...             isomericSmiles=True)
           '[1*]CN'
        
        
        
    
    """
def ReplaceSidechains(mol: Mol, coreQuery: Mol, useChirality: bool = False) -> rdkit.Chem.Mol:
    """
        Replaces sidechains in a molecule with dummy atoms for their attachment points.
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - coreQuery: the molecule to be used as a substructure query for recognizing the core
        
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a new molecule with the sidechains removed
        
          NOTES:
        
            - The original molecule is *not* modified.
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - ReplaceSidechains('CCC1CCC1','C1CCC1') -> '[Xa]C1CCC1'
        
            - ReplaceSidechains('CCC1CC1','C1CCC1') -> ''
        
            - ReplaceSidechains('C1CC2C1CCC2','C1CCC1') -> '[Xa]C1CCC1[Xb]'
        
        
    
    """
def ReplaceSubstructs(mol: Mol, query: Mol, replacement: Mol, replaceAll: bool = False, replacementConnectionPoint: int = 0, useChirality: bool = False) -> typing.Any:
    """
        Replaces atoms matching a substructure query in a molecule
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
        
            - query: the molecule to be used as a substructure query
        
            - replacement: the molecule to be used as the replacement
        
            - replaceAll: (optional) if this toggle is set, all substructures matching
              the query will be replaced in a single result, otherwise each result will
              contain a separate replacement.
              Default value is False (return multiple replacements)
            - replacementConnectionPoint: (optional) index of the atom in the replacement that
              the bond should be made to.
            - useChirality: (optional) match the substructure query using chirality
        
          RETURNS: a tuple of new molecules with the substructures replaced removed
        
          NOTES:
        
            - The original molecule is *not* modified.
            - A bond is only formed to the remaining atoms, if any, that were bonded 
              to the first atom in the substructure query. (For finer control over
              substructure replacement, consider using ChemicalReaction.)
        
          EXAMPLES:
        
           The following examples substitute SMILES/SMARTS strings for molecules, you'd have
           to actually use molecules:
        
            - ReplaceSubstructs('CCOC','O[CH3]','NC') -> ('CCNC',)
        
            - ReplaceSubstructs('COCCOC','O[CH3]','NC') -> ('COCCNC','CNCCOC')
        
            - ReplaceSubstructs('COCCOC','O[CH3]','NC',True) -> ('CNCCNC',)
        
            - ReplaceSubstructs('COCCOC','O[CH3]','CN',True,1) -> ('CNCCNC',)
        
            - ReplaceSubstructs('CCOC','[CH3]O','NC') -> ('CC.CN',)
        
        
    
    """
def SanitizeMol(mol: Mol, sanitizeOps: int = ..., catchErrors: bool = False) -> SanitizeFlags:
    """
        Kekulize, check valencies, set aromaticity, conjugation and hybridization
        
            - The molecule is modified in place.
        
            - If sanitization fails, an exception will be thrown unless catchErrors is set
        
          ARGUMENTS:
        
            - mol: the molecule to be modified
            - sanitizeOps: (optional) sanitization operations to be carried out
              these should be constructed by or'ing together the
              operations in rdkit.Chem.SanitizeFlags
            - catchErrors: (optional) if provided, instead of raising an exception
              when sanitization fails (the default behavior), the 
              first operation that failed (as defined in rdkit.Chem.SanitizeFlags)
              is returned. Zero is returned on success.
        
        
    
    """
def SetAllowNontetrahedralChirality(val: bool) -> None:
    """
        toggles recognition of non-tetrahedral chirality from 3D structures
    
    """
def SetAromaticity(mol: Mol, model: AromaticityModel = ...) -> None:
    """
        does aromaticity perception
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - model: the model to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        
    
    """
def SetBondStereoFromDirections(mol: Mol) -> None:
    """
        Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
        
        
    
    """
def SetConjugation(mol: Mol) -> None:
    """
        finds conjugated bonds
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        
    
    """
def SetDoubleBondNeighborDirections(mol: Mol, conf: typing.Any = None) -> None:
    """
        Uses the stereo info on double bonds to set the directions of neighboring single bonds
                
          ARGUMENTS:
          
            - mol: the molecule to be modified
        
        
    
    """
def SetGenericQueriesFromProperties(mol: Mol, useAtomLabels: bool = True, useSGroups: bool = True) -> None:
    """
        documentation
    
    """
def SetHybridization(mol: Mol) -> None:
    """
        Assigns hybridization states to atoms
        
          ARGUMENTS:
        
            - mol: the molecule to use
        
          NOTES:
        
            - The molecule is modified in place.
        
        
    
    """
def SetTerminalAtomCoords(mol: Mol, idx: int, otherIdx: int) -> None:
    """
        Sets Cartesian coordinates for a terminal atom.
        
          Useful for growing an atom off a molecule with sensible 
          coordinates based on the geometry of the neighbor.
        
          NOTE: this sets the appropriate coordinates in all of the molecule's conformers 
          ARGUMENTS:
        
            - mol: the molecule the atoms belong to.
            - idx: index of the terminal atom whose coordinates are set.
            - mol: index of the bonded neighbor atom.
        
          RETURNS: Nothing
        
        
    
    """
def SetUseLegacyStereoPerception(val: bool) -> None:
    """
        sets usage of the legacy stereo perception code
    
    """
def SimplifyEnhancedStereo(mol: Mol, removeAffectedStereoGroups: bool = True) -> None:
    """
        Simplifies the stereochemical representation of a molecule where all
        specified stereocenters are in the same StereoGroup
        
          Arguments:
           - mol: molecule to modify
           - removeAffectedStereoGroups: if set then the affected StereoGroups will be removed
        
        If all specified stereocenters are in the same AND or OR stereogroup, a
        moleculeNote property will be set on the molecule with the value "AND
        enantiomer" or "OR enantiomer". CIP labels, if present, are removed.
        
    
    """
def SortMatchesByDegreeOfCoreSubstitution(mol: Mol, core: Mol, matches: typing.Any) -> typing.Any:
    """
        Postprocesses the results of a mol.GetSubstructMatches(core) call 
        where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups). 
        It returns a copy of matches sorted by decreasing number of non-hydrogen matches 
        to the terminal dummy atoms.
        
          ARGUMENTS:
        
            - mol: the molecule GetSubstructMatches was run on
        
            - core: the molecule used as a substructure query
        
            - matches: the result returned by GetSubstructMatches
        
          RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches 
                   to the terminal dummy atoms
        
    
    """
def SplitMolByPDBChainId(mol: Mol, whiteList: typing.Any = None, negateList: bool = False) -> dict:
    """
        Splits a molecule into pieces based on PDB chain information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by chain id with molecules as the values
        
        
    
    """
def SplitMolByPDBResidues(mol: Mol, whiteList: typing.Any = None, negateList: bool = False) -> dict:
    """
        Splits a molecule into pieces based on PDB residue information.
        
          ARGUMENTS:
        
            - mol: the molecule to use
            - whiteList: only residues in this list will be returned
            - negateList: if set, negates the white list inclusion logic
        
          RETURNS: a dictionary keyed by residue name with molecules as the values
        
        
    
    """
def TranslateChiralFlagToStereoGroups(mol: Mol, zeroFlagGroupType: StereoGroupType = ...) -> None:
    """
        Generate enhanced stereo groups based on the status of the chiral flag property.
        
          Arguments:
           - mol: molecule to be modified
           - zeroFlagGroupType: how to handle non-grouped stereo centers when the
                  chiral flag is set to zero
        
          If the chiral flag is set to a value of 1 then all specified tetrahedral
          chiral centers which are not already in StereoGroups will be added to an
          ABS StereoGroup.
        
          If the chiral flag is set to a value of 0 then all specified tetrahedral
          chiral centers will be added to a StereoGroup of the type zeroFlagGroupType
        
          If there is no chiral flag set (i.e. the property is not present), the
          molecule will not be modified.
    
    """
def UnfoldedRDKFingerprintCountBased(mol: Mol, minPath: int = 1, maxPath: int = 7, useHs: bool = True, branchedPaths: bool = True, useBondOrder: bool = True, atomInvariants: typing.Any = 0, fromAtoms: typing.Any = 0, atomBits: typing.Any = None, bitInfo: typing.Any = None) -> ULongSparseIntVect:
    """
        Returns an unfolded count-based version of the RDKit fingerprint for a molecule
        
        ARGUMENTS:
            
                - mol: the molecule to use
            
                - minPath: (optional) minimum number of bonds to include in the subgraphs
                  Defaults to 1.
            
                - maxPath: (optional) maximum number of bonds to include in the subgraphs
                  Defaults to 7.
            
                - useHs: (optional) include paths involving Hs in the fingerprint if the molecule
                  has explicit Hs.
                  Defaults to True.
            
                - branchedPaths: (optional) if set both branched and unbranched paths will be
                  used in the fingerprint.
                  Defaults to True.
            
                - useBondOrder: (optional) if set both bond orders will be used in the path hashes
                  Defaults to True.
            
                - atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
                  Defaults to empty.
            
                - fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs 
                  starting from these atoms will be used.
                  Defaults to empty.
            
                - atomBits: (optional) an empty list. If provided, the result will contain a list 
                  containing the bits each atom sets.
                  Defaults to empty.
            
                - bitInfo: (optional) an empty dict. If provided, the result will contain a dict 
                  with bits as keys and corresponding bond paths as values.
                  Defaults to empty.
             
             
        
    
    """
def WedgeBond(bond: Bond, fromAtomIdx: int, conf: Conformer) -> None:
    """
        Set the wedging on an individual bond from a molecule.
           The wedging scheme used is that from Mol files.
          ARGUMENTS:
            - bond: the bond to update
            - atom ID: the atom from which to do the wedging
            - conformer: the conformer to use to determine wedge direction
        
    
    """
def WedgeMolBonds(mol: Mol, conformer: Conformer, params: BondWedgingParameters = None) -> None:
    """
        Set the wedging on single bonds in a molecule.
           The wedging scheme used is that from Mol files.
        
          ARGUMENTS:
        
            - molecule: the molecule to update
            - conformer: the conformer to use to determine wedge direction
        
        
        
    
    """
def _TestSetProps(mol: Mol) -> None:
    """
    """
@typing.overload
def molzip(a: Mol, b: Mol, params: MolzipParams = ...) -> rdkit.Chem.Mol:
    """
        zip together two molecules using the given matching parameters
    
    """
@typing.overload
def molzip(a: Mol, params: MolzipParams = ...) -> rdkit.Chem.Mol:
    """
        zip together multiple molecules within a combined molecule using the given matching parameters
    
    """
@typing.overload
def molzip(row: dict, params: MolzipParams = ...) -> rdkit.Chem.Mol:
    """
        zip an RGroupRow together to recreate the original molecule.  This correctly handles
        broken cycles that can occur in decompositions.
         example:
        
          >>> from rdkit import Chem
          >>> from rdkit.Chem import rdRGroupDecomposition as rgd
          >>> core = Chem.MolFromSmiles('CO')
          >>> mols = [Chem.MolFromSmiles('C1NNO1')]
          >>> rgroups, unmatched = rgd.RGroupDecompose(core, mols)
          >>> for rgroup in rgroups:
          ...     mol = rgd.molzip(rgroup)
        
        
    
    """
def molzipFragments(mols: typing.Any, params: MolzipParams = ...) -> rdkit.Chem.Mol:
    """
        zip together multiple molecules from an R group decomposition 
        using the given matching parameters.  The first molecule in the list
        must be the core
    
    """
ADJUST_IGNOREALL: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREALL
ADJUST_IGNORECHAINS: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORECHAINS
ADJUST_IGNOREDUMMIES: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREDUMMIES
ADJUST_IGNOREMAPPED: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNOREMAPPED
ADJUST_IGNORENONDUMMIES: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONDUMMIES
ADJUST_IGNORENONE: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORENONE
ADJUST_IGNORERINGS: AdjustQueryWhichFlags  # value = rdkit.Chem.rdmolops.AdjustQueryWhichFlags.ADJUST_IGNORERINGS
AROMATICITY_CUSTOM: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_CUSTOM
AROMATICITY_DEFAULT: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT
AROMATICITY_MDL: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MDL
AROMATICITY_MMFF94: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_MMFF94
AROMATICITY_RDKIT: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_RDKIT
AROMATICITY_SIMPLE: AromaticityModel  # value = rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_SIMPLE
LayeredFingerprint_substructLayers: int = 7
SANITIZE_ADJUSTHS: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS
SANITIZE_ALL: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL
SANITIZE_CLEANUP: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP
SANITIZE_CLEANUPATROPISOMERS: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPATROPISOMERS
SANITIZE_CLEANUPCHIRALITY: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
SANITIZE_CLEANUP_ORGANOMETALLICS: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP_ORGANOMETALLICS
SANITIZE_FINDRADICALS: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
SANITIZE_KEKULIZE: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
SANITIZE_NONE: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
SANITIZE_PROPERTIES: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
SANITIZE_SETAROMATICITY: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
SANITIZE_SETCONJUGATION: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
SANITIZE_SETHYBRIDIZATION: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
SANITIZE_SYMMRINGS: SanitizeFlags  # value = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
_LayeredFingerprint_version: str = '0.7.0'
_PatternFingerprint_version: str = '1.0.0'
_RDKFingerprint_version: str = '2.0.0'
