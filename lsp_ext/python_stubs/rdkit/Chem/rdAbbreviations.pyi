# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing functions for working with molecular abbreviations
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['AbbreviationDefinition', 'CondenseAbbreviationSubstanceGroups', 'CondenseMolAbbreviations', 'GetDefaultAbbreviations', 'GetDefaultLinkers', 'LabelMolAbbreviations', 'ParseAbbreviations', 'ParseLinkers']
class AbbreviationDefinition(Boost.Python.instance):
    """
    Abbreviation Definition
    """
    __instance_size__: typing.ClassVar[int] = 168
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def displayLabel(self) -> str:
        """the label in a drawing when the bond comes from the right (default: '')"""
    @displayLabel.setter
    def displayLabel(self, value: str) -> None: ...
    @property
    def displayLabelW(self) -> str:
        """the label in a drawing when the bond comes from the west (default: '')"""
    @displayLabelW.setter
    def displayLabelW(self, value: str) -> None: ...
    @property
    def includesXBonds(self) -> bool:
        """whether or not the abbreviation definition includes bonds to non-abbreviation atoms (default: True)"""
    @includesXBonds.setter
    def includesXBonds(self, value: bool) -> None: ...
    @property
    def label(self) -> str:
        """the label (default: '')"""
    @label.setter
    def label(self, value: str) -> None: ...
    @property
    def mol(*args, **kwargs):
        """
        the query molecule (should have a dummy as the first atom if includesXBonds is true)
        """
    @mol.setter
    def mol(*args, **kwargs):
        ...
class _vectN5RDKit13Abbreviations22AbbreviationDefinitionE(Boost.Python.instance):
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
                boost::python::objects::iterator_range<boost::python::return_internal_reference<1ul, boost::python::default_call_policies>, std::__1::__wrap_iter<RDKit::Abbreviations::AbbreviationDefinition*>> __iter__(boost::python::back_reference<std::__1::vector<RDKit::Abbreviations::AbbreviationDefinition, std::__1::allocator<RDKit::Abbreviations::AbbreviationDefinition>>&>)
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
def CondenseAbbreviationSubstanceGroups(mol: Mol) -> rdkit.Chem.Mol:
    """
        Finds and replaces abbreviation (i.e. "SUP") substance groups in a molecule. The result is not sanitized.
    
    """
def CondenseMolAbbreviations(mol: Mol, abbrevs: typing.Any, maxCoverage: float = 0.4, sanitize: bool = True) -> rdkit.Chem.Mol:
    """
        Finds and replaces abbreviations in a molecule. The result is not sanitized.
    
    """
def GetDefaultAbbreviations() -> ...:
    """
        returns a list of the default abbreviation definitions
    
    """
def GetDefaultLinkers() -> ...:
    """
        returns a list of the default linker definitions
    
    """
def LabelMolAbbreviations(mol: Mol, abbrevs: typing.Any, maxCoverage: float = 0.4) -> rdkit.Chem.Mol:
    """
        Finds abbreviations and adds to them to a molecule as "SUP" SubstanceGroups
    
    """
def ParseAbbreviations(text: str, removeExtraDummies: bool = False, allowConnectionToDummies: bool = False) -> ...:
    """
        Returns a set of abbreviation definitions from a string.  Format of the text data:  A series of lines, each of which contains: label SMARTS displayLabel displayLabelW  Where label is the label used for the abbreviation, SMARTS is the SMARTS definition of the abbreviation, displayLabel is used in drawings to render the abbreviations and displayLabelW is the display label if a bond comes in from the right.  The 'displayLabel' and 'displayLabelW' fields are optional.  Use dummies in the SMARTS to indicate attachment points. The assumption is that the first atom is a dummy (one will be added if this is not true) and that the second atom is the surrogate for the rest of the group.
    
    """
def ParseLinkers(text: str) -> ...:
    """
        Returns a set of linker definitions from a string.  Equivalent to calling ParseAbbreviations(text, True True).
    
    """
