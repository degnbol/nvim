# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing RGroupDecomposition classes and functions.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['AtomIndexLabels', 'AtomMap', 'AtomMapLabels', 'AutoDetect', 'DummyAtomLabels', 'Exhaustive', 'FingerprintVariance', 'GA', 'Greedy', 'GreedyChunks', 'Isotope', 'IsotopeLabels', 'MCS', 'MDLRGroup', 'MDLRGroupLabels', 'MOL_SPTR_VECT', 'Match', 'NoAlignment', 'NoSymmetrization', 'None', 'RGroupCoreAlignment', 'RGroupDecompose', 'RGroupDecomposition', 'RGroupDecompositionParameters', 'RGroupLabelling', 'RGroupLabels', 'RGroupMatching', 'RGroupScore', 'RelabelDuplicateLabels', 'RelabelMappedDummies']
class MOL_SPTR_VECT(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__wrap_iter<boost::shared_ptr<RDKit::ROMol>*>> __iter__(boost::python::back_reference<std::__1::vector<boost::shared_ptr<RDKit::ROMol>, std::__1::allocator<boost::shared_ptr<RDKit::ROMol>>>&>)
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
class RGroupCoreAlignment(Boost.Python.enum):
    MCS: typing.ClassVar[RGroupCoreAlignment]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS
    NoAlignment: typing.ClassVar[RGroupCoreAlignment]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment
#     None: typing.ClassVar[RGroupCoreAlignment]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'None': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None, 'NoAlignment': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment, 'MCS': rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment, 1: rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS}
class RGroupDecomposition(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def Add(self, mol: Mol) -> int:
        """
        """
    def GetMatchingCoreIdx(self, mol: Mol, matches: typing.Any = None) -> int:
        """
        """
    def GetRGroupLabels(self) -> list:
        """
            Return the current list of found rgroups.
            Note, Process() should be called first
        
        """
    def GetRGroupsAsColumns(self, asSmiles: bool = False) -> dict:
        """
            Return the rgroups as columns (note: can be fed directly into a pandas datatable)
              ARGUMENTS:
               - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
                Column structure:
                   columns[rgroup_label] = [ mols_or_smiles ]
            
        
        """
    def GetRGroupsAsRows(self, asSmiles: bool = False) -> list:
        """
            Return the rgroups as rows (note: can be fed directly into a pandas datatable)
              ARGUMENTS:
               - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
                Row structure:
                   rows[idx] = {rgroup_label: molecule_or_smiles}
            
        
        """
    def Process(self) -> bool:
        """
            Process the rgroups (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
        
        """
    def ProcessAndScore(self) -> tuple:
        """
            Process the rgroups and returns the score (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)
        
        """
    @typing.overload
    def __init__(self, cores: typing.Any) -> None:
        ...
    @typing.overload
    def __init__(self, cores: typing.Any, params: RGroupDecompositionParameters) -> None:
        ...
class RGroupDecompositionParameters(Boost.Python.instance):
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
      OPTIONS:
        - RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCS
                               If set to MCS, cores labels are mapped to each other using their
                               Maximum common substructure overlap.
        - RGroupLabels: optionally set where the rgroup labels to use are encoded.
                         RGroupLabels.IsotopeLabels - labels are stored on isotopes
                         RGroupLabels.AtomMapLabels - labels are stored on atommaps
                         RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
                         RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
                         RGroupLabels.AtomIndexLabels - use the atom index as the label
                         RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
                         RGroupLabels.AutoDetect - auto detect the label [default]
           Note: in all cases, any rgroups found on unlabelled atoms will be automatically
                  labelled.
        - RGroupLabelling: choose where the rlabels are stored on the decomposition
                            RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
                            RGroupLabelling.Isotope - store rgroups on the isotope
                            RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)
                           default: AtomMap | MDLRGroup
        - onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
        - removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
        - removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
        - removeHydrogensPostMatch: remove all hydrogens from the output molecules
        - allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more
        - doTautomers: match all tautomers of a core against each input structure
        - doEnumeration: expand input cores into enumerated mol bundles
        - allowMultipleRGroupsOnUnlabelled: permit more than one rgroup to be attached to an unlabelled core atom
        - allowMultipleCoresInSameMol: permit a core to match more than once in the same molecule if the sets of matched atoms are not equal (default=False)
    """
    __instance_size__: typing.ClassVar[int] = 288
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
    def alignment(self) -> int:
        """default: 1"""
    @alignment.setter
    def alignment(self, value: int) -> None: ...
    @property
    def allowMultipleCoresInSameMol(self) -> bool:
        """default: False"""
    @allowMultipleCoresInSameMol.setter
    def allowMultipleCoresInSameMol(self, value: bool) -> None: ...
    @property
    def allowMultipleRGroupsOnUnlabelled(self) -> bool:
        """default: False"""
    @allowMultipleRGroupsOnUnlabelled.setter
    def allowMultipleRGroupsOnUnlabelled(self, value: bool) -> None: ...
    @property
    def allowNonTerminalRGroups(self) -> bool:
        """default: False"""
    @allowNonTerminalRGroups.setter
    def allowNonTerminalRGroups(self, value: bool) -> None: ...
    @property
    def chunkSize(self) -> int:
        """default: 5"""
    @chunkSize.setter
    def chunkSize(self, value: int) -> None: ...
    @property
    def doEnumeration(self) -> bool:
        """default: False"""
    @doEnumeration.setter
    def doEnumeration(self, value: bool) -> None: ...
    @property
    def doTautomers(self) -> bool:
        """default: False"""
    @doTautomers.setter
    def doTautomers(self, value: bool) -> None: ...
    @property
    def gaMaximumOperations(self) -> int:
        """default: -1"""
    @gaMaximumOperations.setter
    def gaMaximumOperations(self, value: int) -> None: ...
    @property
    def gaNumberOperationsWithoutImprovement(self) -> int:
        """default: -1"""
    @gaNumberOperationsWithoutImprovement.setter
    def gaNumberOperationsWithoutImprovement(self, value: int) -> None: ...
    @property
    def gaNumberRuns(self) -> int:
        """default: 1"""
    @gaNumberRuns.setter
    def gaNumberRuns(self, value: int) -> None: ...
    @property
    def gaParallelRuns(self) -> bool:
        """default: True"""
    @gaParallelRuns.setter
    def gaParallelRuns(self, value: bool) -> None: ...
    @property
    def gaPopulationSize(self) -> int:
        """default: -1"""
    @gaPopulationSize.setter
    def gaPopulationSize(self, value: int) -> None: ...
    @property
    def gaRandomSeed(self) -> int:
        """default: -1"""
    @gaRandomSeed.setter
    def gaRandomSeed(self, value: int) -> None: ...
    @property
    def includeTargetMolInResults(self) -> bool:
        """default: False"""
    @includeTargetMolInResults.setter
    def includeTargetMolInResults(self, value: bool) -> None: ...
    @property
    def labels(self) -> int:
        """default: 255"""
    @labels.setter
    def labels(self, value: int) -> None: ...
    @property
    def matchingStrategy(self) -> int:
        """default: 2"""
    @matchingStrategy.setter
    def matchingStrategy(self, value: int) -> None: ...
    @property
    def onlyMatchAtRGroups(self) -> bool:
        """default: False"""
    @onlyMatchAtRGroups.setter
    def onlyMatchAtRGroups(self, value: bool) -> None: ...
    @property
    def removeAllHydrogenRGroups(self) -> bool:
        """default: True"""
    @removeAllHydrogenRGroups.setter
    def removeAllHydrogenRGroups(self, value: bool) -> None: ...
    @property
    def removeAllHydrogenRGroupsAndLabels(self) -> bool:
        """default: True"""
    @removeAllHydrogenRGroupsAndLabels.setter
    def removeAllHydrogenRGroupsAndLabels(self, value: bool) -> None: ...
    @property
    def removeHydrogensPostMatch(self) -> bool:
        """default: True"""
    @removeHydrogensPostMatch.setter
    def removeHydrogensPostMatch(self, value: bool) -> None: ...
    @property
    def rgroupLabelling(self) -> int:
        """default: 5"""
    @rgroupLabelling.setter
    def rgroupLabelling(self, value: int) -> None: ...
    @property
    def scoreMethod(self) -> int:
        """default: 1"""
    @scoreMethod.setter
    def scoreMethod(self, value: int) -> None: ...
    @property
    def substructMatchParams(*args, **kwargs):
        ...
    @property
    def timeout(self) -> float:
        """default: -1.0"""
    @timeout.setter
    def timeout(self, value: float) -> None: ...
class RGroupLabelling(Boost.Python.enum):
    AtomMap: typing.ClassVar[RGroupLabelling]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap
    Isotope: typing.ClassVar[RGroupLabelling]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope
    MDLRGroup: typing.ClassVar[RGroupLabelling]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'AtomMap': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap, 'Isotope': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope, 'MDLRGroup': rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap, 2: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope, 4: rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup}
class RGroupLabels(Boost.Python.enum):
    AtomIndexLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels
    AtomMapLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels
    AutoDetect: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect
    DummyAtomLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels
    IsotopeLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels
    MDLRGroupLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels
    RelabelDuplicateLabels: typing.ClassVar[RGroupLabels]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'IsotopeLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels, 'AtomMapLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels, 'AtomIndexLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels, 'RelabelDuplicateLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels, 'MDLRGroupLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels, 'DummyAtomLabels': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels, 'AutoDetect': rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels, 2: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels, 4: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels, 8: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels, 16: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels, 32: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels, 255: rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect}
class RGroupMatching(Boost.Python.enum):
    Exhaustive: typing.ClassVar[RGroupMatching]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive
    GA: typing.ClassVar[RGroupMatching]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA
    Greedy: typing.ClassVar[RGroupMatching]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy
    GreedyChunks: typing.ClassVar[RGroupMatching]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks
    NoSymmetrization: typing.ClassVar[RGroupMatching]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Greedy': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy, 'GreedyChunks': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks, 'Exhaustive': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive, 'NoSymmetrization': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization, 'GA': rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy, 2: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks, 4: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive, 8: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization, 16: rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA}
class RGroupScore(Boost.Python.enum):
    FingerprintVariance: typing.ClassVar[RGroupScore]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance
    Match: typing.ClassVar[RGroupScore]  # value = rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Match': rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match, 'FingerprintVariance': rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match, 4: rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance}
def RGroupDecompose(cores: typing.Any, mols: typing.Any, asSmiles: bool = False, asRows: bool = True, options: RGroupDecompositionParameters = ...) -> typing.Any:
    """
        Decompose a collection of molecules into their Rgroups
          ARGUMENTS:
            - cores: a set of cores from most to least specific.
                     See RGroupDecompositionParameters for more details
                     on how the cores can be labelled
            - mols: the molecules to be decomposed
            - asSmiles: if True return smiles strings, otherwise return molecules [default: False]
            - asRows: return the results as rows (default) otherwise return columns
            - options: RGroupDecompositionParameters object that defines the parameters for the decomposition.
                     See RGroupDecompositionParameters for defaults
        
          RETURNS: row_or_column_results, unmatched
        
            Row structure:
               rows[idx] = {rgroup_label: molecule_or_smiles}
            Column structure:
               columns[rgroup_label] = [ mols_or_smiles ]
        
            unmatched is a vector of indices in the input mols that were not matched.
        
    
    """
def RelabelMappedDummies(mol: Mol, inputLabels: int = ..., outputLabels: int = ...) -> None:
    """
        Relabel dummy atoms bearing an R-group mapping (as
        atom map number, isotope or MDLRGroup label) such that
        they will be displayed by the rendering code as R# rather
        than #*, *:#, #*:#, etc. By default, only the MDLRGroup label
        is retained on output; this may be configured through the
        outputLabels parameter.
        In case there are multiple potential R-group mappings,
        the priority on input is Atom map number > Isotope > MDLRGroup.
        The inputLabels parameter allows to configure which mappings
        are taken into consideration.
        
    
    """
AtomIndexLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomIndexLabels
AtomMap: RGroupLabelling  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.AtomMap
AtomMapLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AtomMapLabels
AutoDetect: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.AutoDetect
DummyAtomLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.DummyAtomLabels
Exhaustive: RGroupMatching  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Exhaustive
FingerprintVariance: RGroupScore  # value = rdkit.Chem.rdRGroupDecomposition.RGroupScore.FingerprintVariance
GA: RGroupMatching  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GA
Greedy: RGroupMatching  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.Greedy
GreedyChunks: RGroupMatching  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.GreedyChunks
Isotope: RGroupLabelling  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.Isotope
IsotopeLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.IsotopeLabels
MCS: RGroupCoreAlignment  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.MCS
MDLRGroup: RGroupLabelling  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabelling.MDLRGroup
MDLRGroupLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.MDLRGroupLabels
Match: RGroupScore  # value = rdkit.Chem.rdRGroupDecomposition.RGroupScore.Match
NoAlignment: RGroupCoreAlignment  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.NoAlignment
NoSymmetrization: RGroupMatching  # value = rdkit.Chem.rdRGroupDecomposition.RGroupMatching.NoSymmetrization
# None: RGroupCoreAlignment  # value = rdkit.Chem.rdRGroupDecomposition.RGroupCoreAlignment.None
RelabelDuplicateLabels: RGroupLabels  # value = rdkit.Chem.rdRGroupDecomposition.RGroupLabels.RelabelDuplicateLabels
