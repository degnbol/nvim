# fix_pybind_stubs: rdkit 2026.3.5
from __future__ import annotations
import rdkit.Chem
import typing
from .FilterMatchOps import *
__all__: list[str] = ['ExclusionList', 'FilterCatalog', 'FilterCatalogCanSerialize', 'FilterCatalogEntry', 'FilterCatalogEntryList', 'FilterCatalogListOfEntryList', 'FilterCatalogParams', 'FilterHierarchyMatcher', 'FilterMatch', 'FilterMatchOps', 'FilterMatcherBase', 'GetFlattenedFunctionalGroupHierarchy', 'GetFunctionalGroupHierarchy', 'IntPair', 'MolList', 'PythonFilterMatcher', 'RunFilterCatalog', 'SmartsMatcher', 'VectFilterMatch']
class ExclusionList(FilterMatcherBase):
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddPattern(self, base: FilterMatcherBase) -> None:
        """
            Add a FilterMatcherBase that should not appear in a molecule
        
        """
    def SetExclusionPatterns(self, list: typing.Any) -> None:
        """
            Set a list of FilterMatcherBases that should not appear in a molecule
        
        """
    def __init__(self) -> None:
        ...
class FilterCatalog(Boost.Python.instance):
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 72
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def AddEntry(entry: FilterCatalog, updateFPLength: FilterCatalogEntry = False) -> None:
        """
            Add a FilterCatalogEntry to the catalog
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetEntry(self, idx: int) -> FilterCatalogEntry:
        """
            Return the FilterCatalogEntry at the specified index
        
        """
    def GetEntryWithIdx(self, idx: int) -> FilterCatalogEntry:
        """
            Return the FilterCatalogEntry at the specified index
        
        """
    def GetFilterMatches(self, mol: Mol) -> VectFilterMatch:
        """
            Return every matching filter from all catalog entries that match mol
        
        """
    def GetFirstMatch(self, mol: Mol) -> FilterCatalogEntry:
        """
            Return the first catalog entry that matches mol
        
        """
    def GetMatches(self, mol: Mol) -> FilterCatalogEntryList:
        """
            Return all catalog entries that match mol
        
        """
    def GetNumEntries(self) -> int:
        """
            Returns the number of entries in the catalog
        
        """
    def HasMatch(self, mol: Mol) -> bool:
        """
            Returns True if the catalog has an entry that matches mol
        
        """
    def RemoveEntry(self, obj: typing.Any) -> bool:
        """
            Remove the given entry from the catalog
        
        """
    def Serialize(self) -> typing.Any:
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
    def __init__(self, binStr: str) -> None:
        ...
    @typing.overload
    def __init__(self, params: FilterCatalogParams) -> None:
        ...
    @typing.overload
    def __init__(self, catalogs: FilterCatalogs) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
class FilterCatalogEntry(Boost.Python.instance):
    """
    FilterCatalogEntry
    A filter catalog entry is an entry in a filter catalog.
    Each filter is named and is used to flag a molecule usually for some
    undesirable property.
    
    For example, a PAINS (Pan Assay INterference) catalog entry be appear as
    follows:
    
    >>> from rdkit.Chem.FilterCatalog import *
    >>> params = FilterCatalogParams()
    >>> params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    True
    >>> catalog = FilterCatalog(params)
    >>> mol = Chem.MolFromSmiles('O=C(Cn1cnc2c1c(=O)n(C)c(=O)n2C)N/N=C/c1c(O)ccc2c1cccc2')
    >>> entry = catalog.GetFirstMatch(mol)
    >>> print (entry.GetProp('Scope'))
    PAINS filters (family A)
    >>> print (entry.GetDescription())
    hzone_phenol_A(479)
    
    
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def ClearProp(self, key: typing.Any) -> None:
        """
        """
    def GetDescription(self) -> str:
        """
            Get the description of the catalog entry
        
        """
    def GetFilterMatches(self, mol: Mol) -> VectFilterMatch:
        """
            Retrieve the list of filters that match the molecule
        
        """
    def GetProp(self, key: typing.Any) -> str:
        """
        """
    def GetPropList(self) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
        """
    def HasFilterMatch(self, mol: Mol) -> bool:
        """
            Returns True if the catalog entry contains filters that match the molecule
        
        """
    def IsValid(self) -> bool:
        """
        """
    def Serialize(self) -> typing.Any:
        """
        """
    def SetDescription(self, description: str) -> None:
        """
            Set the description of the catalog entry
        
        """
    def SetProp(self, key: typing.Any, val: str) -> None:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, name: str, matcher: FilterMatcherBase) -> None:
        ...
class FilterCatalogEntryList(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__wrap_iter<boost::shared_ptr<RDKit::FilterCatalogEntry const>*>> __iter__(boost::python::back_reference<std::__1::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::__1::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const>>>&>)
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
class FilterCatalogListOfEntryList(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<std::__1::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::__1::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const>>>*>> __iter__(boost::python::back_reference<std::__1::vector<std::__1::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::__1::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const>>>, std::__1::allocator<std::__1::vector<boost::shared_ptr<RDKit::FilterCatalogEntry const>, std::__1::allocator<boost::shared_ptr<RDKit::FilterCatalogEntry const>>>>>&>)
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
class FilterCatalogParams(Boost.Python.instance):
    class FilterCatalogs(Boost.Python.enum):
        ALL: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.ALL
        BRENK: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.BRENK
        CHEMBL: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL
        CHEMBL_BMS: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_BMS
        CHEMBL_Dundee: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Dundee
        CHEMBL_Glaxo: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Glaxo
        CHEMBL_Inpharmatica: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Inpharmatica
        CHEMBL_LINT: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_LINT
        CHEMBL_MLSMR: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_MLSMR
        CHEMBL_SureChEMBL: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_SureChEMBL
        NIH: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.NIH
        PAINS: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS
        PAINS_A: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_A
        PAINS_B: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_B
        PAINS_C: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_C
        ZINC: typing.ClassVar[FilterCatalogs]  # value = rdkit.Chem.rdfiltercatalog.FilterCatalogs.ZINC
        __slots__: typing.ClassVar[tuple] = tuple()
        names: typing.ClassVar[dict]  # value = {'PAINS_A': rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_A, 'PAINS_B': rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_B, 'PAINS_C': rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_C, 'PAINS': rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS, 'BRENK': rdkit.Chem.rdfiltercatalog.FilterCatalogs.BRENK, 'NIH': rdkit.Chem.rdfiltercatalog.FilterCatalogs.NIH, 'ZINC': rdkit.Chem.rdfiltercatalog.FilterCatalogs.ZINC, 'CHEMBL_Glaxo': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Glaxo, 'CHEMBL_Dundee': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Dundee, 'CHEMBL_BMS': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_BMS, 'CHEMBL_SureChEMBL': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_SureChEMBL, 'CHEMBL_MLSMR': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_MLSMR, 'CHEMBL_Inpharmatica': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Inpharmatica, 'CHEMBL_LINT': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_LINT, 'CHEMBL': rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL, 'ALL': rdkit.Chem.rdfiltercatalog.FilterCatalogs.ALL}
        values: typing.ClassVar[dict]  # value = {2: rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_A, 4: rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_B, 8: rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS_C, 14: rdkit.Chem.rdfiltercatalog.FilterCatalogs.PAINS, 16: rdkit.Chem.rdfiltercatalog.FilterCatalogs.BRENK, 32: rdkit.Chem.rdfiltercatalog.FilterCatalogs.NIH, 64: rdkit.Chem.rdfiltercatalog.FilterCatalogs.ZINC, 128: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Glaxo, 256: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Dundee, 512: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_BMS, 1024: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_SureChEMBL, 2048: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_MLSMR, 4096: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_Inpharmatica, 8192: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL_LINT, 16256: rdkit.Chem.rdfiltercatalog.FilterCatalogs.CHEMBL, 16382: rdkit.Chem.rdfiltercatalog.FilterCatalogs.ALL}
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddCatalog(self, catalogs: FilterCatalogs) -> bool:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, catalogs: FilterCatalogs) -> None:
        ...
class FilterHierarchyMatcher(FilterMatcherBase):
    """
    Hierarchical Filter
     basic constructors: 
       FilterHierarchyMatcher( matcher )
       where can be any FilterMatcherBase (SmartsMatcher, etc)
     FilterHierarchyMatcher's have children and can form matching
      trees.  then GetFilterMatches is called, the most specific (
      i.e. lowest node in a branch) is returned.
    
     n.b. A FilterHierarchicalMatcher of functional groups is returned
      by calling GetFunctionalGroupHierarchy()
    
    >>> from rdkit.Chem import MolFromSmiles
    >>> from rdkit.Chem.FilterCatalog import *
    >>> functionalGroups = GetFunctionalGroupHierarchy()
    >>> [match.filterMatch.GetName() 
    ...     for match in functionalGroups.GetFilterMatches(
    ...         MolFromSmiles('c1ccccc1Cl'))]
    ['Halogen.Aromatic', 'Halogen.NotFluorine.Aromatic']
    
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddChild(self, hierarchy: FilterHierarchyMatcher) -> FilterHierarchyMatcher:
        """
            Add a child node to this hierarchy.
        
        """
    def SetPattern(self, matcher: FilterMatcherBase) -> None:
        """
            Set the filtermatcher pattern for this node.  An empty node is considered a root node and passes along the matches to the children.
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, matcher: FilterMatcherBase) -> None:
        ...
class FilterMatch(Boost.Python.instance):
    """
    Object that holds the result of running FilterMatcherBase::GetMatches
    
     - filterMatch holds the FilterMatchBase that triggered the match
     - atomParis holds the [ (query_atom_idx, target_atom_idx) ] pairs for the matches.
    
    
    Note that some matches may not have atom pairs (especially matches that use FilterMatchOps.Not
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, filter: FilterMatcherBase, atomPairs: MatchTypeVect) -> None:
        ...
    @property
    def atomPairs(*args, **kwargs):
        ...
    @property
    def filterMatch(*args, **kwargs):
        ...
class FilterMatcherBase(Boost.Python.instance):
    """
    Base class for matching molecules to filters.
    
     A FilterMatcherBase supplies the following API 
     - IsValid() returns True if the matcher is valid for use, False otherwise
     - HasMatch(mol) returns True if the molecule matches the filter
     - GetMatches(mol) -> [FilterMatch, FilterMatch] returns all the FilterMatch data
           that matches the molecule
    
    
    print( FilterMatcherBase ) will print user-friendly information about the filter
    Note that a FilterMatcherBase can be combined from may FilterMatcherBases
    This is why GetMatches can return multiple FilterMatcherBases.
    >>> from rdkit.Chem.FilterCatalog import *
    >>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', 0, 1)
    >>> oxygen_matcher = SmartsMatcher('Oxygen', '[#8]', 0, 1)
    >>> co_matcher = FilterMatchOps.Or(carbon_matcher, oxygen_matcher)
    >>> mol = Chem.MolFromSmiles('C')
    >>> matches = co_matcher.GetMatches(mol)
    >>> len(matches)
    1
    >>> print(matches[0].filterMatch)
    Carbon
    
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
    def GetMatches(self, mol: Mol) -> VectFilterMatch:
        """
            Returns the list of matching subfilters mol matches any filter
        
        """
    def GetName(self) -> str:
        """
        """
    def HasMatch(self, mol: Mol) -> bool:
        """
            Returns True if mol matches the filter
        
        """
    def IsValid(self) -> bool:
        """
            Return True if the filter matcher is valid, False otherwise
        
        """
    def __str__(self) -> str:
        """
        """
class IntPair(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __getitem__(self, idx: int) -> int:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, query: int, target: int) -> None:
        ...
    @property
    def query(self) -> int:
        """default: 0"""
    @query.setter
    def query(self, value: int) -> None: ...
    @property
    def target(self) -> int:
        """default: 0"""
    @target.setter
    def target(self, value: int) -> None: ...
class MolList(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_value_policy<boost::python::return_by_value, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::ROMol**>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::ROMol*, std::__1::allocator<RDKit::ROMol*>>&>)
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
class PythonFilterMatcher(FilterMatcherBase):
    __instance_size__: typing.ClassVar[int] = 88
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, callback: typing.Any) -> None:
        ...
class SmartsMatcher(FilterMatcherBase):
    """
    Smarts Matcher Filter
     basic constructors: 
       SmartsMatcher( name, smarts_pattern, minCount=1, maxCount=UINT_MAX )
       SmartsMatcher( name, molecule, minCount=1, maxCount=UINT_MAX )
    
      note: If the supplied smarts pattern is not valid, the IsValid() function will
       return False
    >>> from rdkit.Chem.FilterCatalog import *
    >>> minCount, maxCount = 1,2
    >>> carbon_matcher = SmartsMatcher('Carbon', '[#6]', minCount, maxCount)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CC')))
    True
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    False
    >>> carbon_matcher.SetMinCount(2)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('C')))
    False
    >>> carbon_matcher.SetMaxCount(3)
    >>> print (carbon_matcher.HasMatch(Chem.MolFromSmiles('CCC')))
    True
    
    """
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetMaxCount(self) -> int:
        """
            Get the maximum times pattern can appear for the filter to match
        
        """
    def GetMinCount(self) -> int:
        """
            Get the minimum times pattern must appear for the filter to match
        
        """
    def GetPattern(self) -> rdkit.Chem.Mol:
        """
        """
    def IsValid(self) -> bool:
        """
            Returns True if the SmartsMatcher is valid
        
        """
    def SetMaxCount(self, count: int) -> None:
        """
            Set the maximum times pattern can appear for the filter to match
        
        """
    def SetMinCount(self, count: int) -> None:
        """
            Set the minimum times pattern must appear to match
        
        """
    @typing.overload
    def SetPattern(self, pat: Mol) -> None:
        """
            Set the pattern molecule for the SmartsMatcher
        
        """
    @typing.overload
    def SetPattern(self, pat: str) -> None:
        """
            Set the smarts pattern for the Smarts Matcher (warning: MinimumCount is not reset)
        
        """
    @typing.overload
    def __init__(self, name: str) -> None:
        ...
    @typing.overload
    def __init__(self, rhs: Mol) -> None:
        ...
    @typing.overload
    def __init__(self, name: str, mol: Mol, minCount: int = 1, maxCount: int = 4294967295) -> None:
        ...
    @typing.overload
    def __init__(self, name: str, smarts: str, minCount: int = 1, maxCount: int = 4294967295) -> None:
        ...
class VectFilterMatch(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::FilterMatch*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::FilterMatch, std::__1::allocator<RDKit::FilterMatch>>&>)
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
def FilterCatalogCanSerialize() -> bool:
    """
        Returns True if the FilterCatalog is serializable (requires boost serialization
    
    """
def GetFlattenedFunctionalGroupHierarchy(normalized: bool = False) -> dict:
    """
        Returns the flattened functional group hierarchy as a dictionary  of name:ROMOL_SPTR substructure items
    
    """
def GetFunctionalGroupHierarchy() -> FilterCatalog:
    """
        Returns the functional group hierarchy filter catalog
    
    """
def RunFilterCatalog(filterCatalog: FilterCatalog, smiles: _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE, numThreads: int = 1) -> FilterCatalogListOfEntryList:
    """
        Run the filter catalog on the input list of smiles strings.
        Use numThreads=0 to use all available processors. Returns a vector of vectors.  For each input smiles, a vector of FilterCatalogEntry objects are returned for each matched filter.  If a molecule matches no filter, the vector will be empty. If a smiles string can't be parsed, a 'Bad smiles' entry is returned.
    
    """
