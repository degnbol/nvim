# fix_pybind_stubs: rdkit 2026.3.5
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['Deprotect', 'DeprotectData', 'DeprotectDataVect', 'DeprotectInPlace', 'GetDeprotections']
class DeprotectData(Boost.Python.instance):
    """
    DeprotectData class, contains a single deprotection reaction and information
    
     deprotectdata.deprotection_class - functional group being protected
     deprotectdata.reaction_smarts - reaction smarts used for deprotection
     deprotectdata.abbreviation - common abbreviation for the protecting group
     deprotectdata.full_name - full name for the protecting group
    
    
    """
    __instance_size__: typing.ClassVar[int] = 160
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, deprotection_class: str, reaction_smarts: str, abbreviation: str, full_name: str) -> None:
        ...
    def isValid(self) -> bool:
        """
            Returns True if the DeprotectData has a valid reaction
        
        """
    @property
    def abbreviation(*args, **kwargs):
        ...
    @property
    def deprotection_class(*args, **kwargs):
        ...
    @property
    def example(*args, **kwargs):
        ...
    @property
    def full_name(*args, **kwargs):
        ...
    @property
    def reaction_smarts(*args, **kwargs):
        ...
class DeprotectDataVect(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::Deprotect::DeprotectData*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::Deprotect::DeprotectData, std::__1::allocator<RDKit::Deprotect::DeprotectData>>&>)
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
def Deprotect(mol: Mol, deprotections: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Return the deprotected version of the molecule.
    
    """
def DeprotectInPlace(mol: Mol, deprotections: typing.Any = None) -> bool:
    """
        Deprotects the molecule in place.
    
    """
def GetDeprotections() -> DeprotectDataVect:
    """
        Return the default list of deprotections
    
    """
