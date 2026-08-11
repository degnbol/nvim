# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing functions for creating a Scaffold Network
"""
from __future__ import annotations
import typing
__all__: list[str] = ['BRICSScaffoldParams', 'CreateScaffoldNetwork', 'EdgeType', 'NetworkEdge', 'NetworkEdge_VECT', 'ScaffoldNetwork', 'ScaffoldNetworkParams', 'UpdateScaffoldNetwork']
class EdgeType(Boost.Python.enum):
    Fragment: typing.ClassVar[EdgeType]  # value = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment
    Generic: typing.ClassVar[EdgeType]  # value = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic
    GenericBond: typing.ClassVar[EdgeType]  # value = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond
    Initialize: typing.ClassVar[EdgeType]  # value = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize
    RemoveAttachment: typing.ClassVar[EdgeType]  # value = rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Fragment': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment, 'Generic': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic, 'GenericBond': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond, 'RemoveAttachment': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment, 'Initialize': rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Fragment, 2: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Generic, 3: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.GenericBond, 4: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.RemoveAttachment, 5: rdkit.Chem.Scaffolds.rdScaffoldNetwork.EdgeType.Initialize}
class NetworkEdge(Boost.Python.instance):
    """
    A scaffold network edge
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
    def __str__(self) -> typing.Any:
        """
        """
    @property
    def beginIdx(*args, **kwargs):
        """
        index of the begin node in node list
        """
    @property
    def endIdx(*args, **kwargs):
        """
        index of the end node in node list
        """
    @property
    def type(*args, **kwargs):
        """
        type of the edge
        """
class NetworkEdge_VECT(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::ScaffoldNetwork::NetworkEdge*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::__1::allocator<RDKit::ScaffoldNetwork::NetworkEdge>>&>)
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
class ScaffoldNetwork(Boost.Python.instance):
    """
    A scaffold network
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 120
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
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
    def __init__(self, pkl: str) -> None:
        ...
    def __setstate__(self, data: tuple) -> None:
        """
        """
    @property
    def counts(*args, **kwargs):
        """
        the number of times each node was encountered while building the network.
        """
    @property
    def edges(*args, **kwargs):
        """
        the sequence of network edges
        """
    @property
    def molCounts(*args, **kwargs):
        """
        the number of moleclues each node was found in.
        """
    @property
    def nodes(*args, **kwargs):
        """
        the sequence of SMILES defining the nodes
        """
class ScaffoldNetworkParams(Boost.Python.instance):
    """
    Scaffold network parameters
    """
    __instance_size__: typing.ClassVar[int] = 64
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, bondBreakerSmartsList: _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE) -> None:
        ...
    @property
    def collectMolCounts(self) -> bool:
        """keep track of the number of molecules each scaffold was found in (default: True)"""
    @collectMolCounts.setter
    def collectMolCounts(self, value: bool) -> None: ...
    @property
    def flattenChirality(self) -> bool:
        """remove chirality and bond stereo when flattening (default: True)"""
    @flattenChirality.setter
    def flattenChirality(self, value: bool) -> None: ...
    @property
    def flattenIsotopes(self) -> bool:
        """remove isotopes when flattening (default: True)"""
    @flattenIsotopes.setter
    def flattenIsotopes(self, value: bool) -> None: ...
    @property
    def flattenKeepLargest(self) -> bool:
        """keep only the largest fragment when doing flattening (default: True)"""
    @flattenKeepLargest.setter
    def flattenKeepLargest(self, value: bool) -> None: ...
    @property
    def includeGenericBondScaffolds(self) -> bool:
        """include scaffolds with all bonds replaced by single bonds (default: False)"""
    @includeGenericBondScaffolds.setter
    def includeGenericBondScaffolds(self, value: bool) -> None: ...
    @property
    def includeGenericScaffolds(self) -> bool:
        """include scaffolds with all atoms replaced by dummies (default: True)"""
    @includeGenericScaffolds.setter
    def includeGenericScaffolds(self, value: bool) -> None: ...
    @property
    def includeNames(self) -> bool:
        """Include molecules names of the input molecules (default: False)"""
    @includeNames.setter
    def includeNames(self, value: bool) -> None: ...
    @property
    def includeScaffoldsWithAttachments(self) -> bool:
        """Include the version of the scaffold with attachment points (default: True)"""
    @includeScaffoldsWithAttachments.setter
    def includeScaffoldsWithAttachments(self, value: bool) -> None: ...
    @property
    def includeScaffoldsWithoutAttachments(self) -> bool:
        """remove attachment points from scaffolds and include the result (default: True)"""
    @includeScaffoldsWithoutAttachments.setter
    def includeScaffoldsWithoutAttachments(self, value: bool) -> None: ...
    @property
    def keepOnlyFirstFragment(self) -> bool:
        """keep only the first fragment from the bond breaking rule (default: True)"""
    @keepOnlyFirstFragment.setter
    def keepOnlyFirstFragment(self, value: bool) -> None: ...
    @property
    def pruneBeforeFragmenting(self) -> bool:
        """Do a pruning/flattening step before starting fragmenting (default: True)"""
    @pruneBeforeFragmenting.setter
    def pruneBeforeFragmenting(self, value: bool) -> None: ...
def BRICSScaffoldParams() -> ScaffoldNetworkParams:
    """
        Returns parameters for generating scaffolds using BRICS fragmentation rules
    
    """
def CreateScaffoldNetwork(mols: typing.Any, params: ScaffoldNetworkParams) -> ScaffoldNetwork:
    """
        create (and return) a new network from a sequence of molecules
    
    """
def UpdateScaffoldNetwork(mols: typing.Any, network: ScaffoldNetwork, params: ScaffoldNetworkParams) -> None:
    """
        update an existing network by adding molecules
    
    """
