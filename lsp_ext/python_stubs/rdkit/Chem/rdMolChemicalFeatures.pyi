# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing from chemical feature and functions to generate the
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['BuildFeatureFactory', 'BuildFeatureFactoryFromString', 'GetAtomMatch', 'MolChemicalFeature', 'MolChemicalFeatureFactory']
class MolChemicalFeature(Boost.Python.instance):
    """
    Class to represent a chemical feature.
        These chemical features may or may not have been derived from molecule object;
        i.e. it is possible to have a chemical feature that was created just from its type
        and location.
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
    def ClearCache(self) -> None:
        """
            Clears the cache used to store position information.
        
        """
    def GetActiveConformer(self) -> int:
        """
            Gets the conformer to use.
        
        """
    def GetAtomIds(self) -> typing.Any:
        """
            Get the IDs of the atoms that participate in the feature
        
        """
    def GetFactory(self) -> MolChemicalFeatureFactory:
        """
            Get the factory used to generate this feature
        
        """
    def GetFamily(self) -> str:
        """
            Get the family to which the feature belongs; donor, acceptor, etc.
        
        """
    def GetId(self) -> int:
        """
            Returns the identifier of the feature
            
        
        """
    def GetMol(self) -> rdkit.Chem.Mol:
        """
            Get the molecule used to derive the features
        
        """
    @typing.overload
    def GetPos(self, confId: int) -> Point3D:
        """
            Get the location of the chemical feature
        
        """
    @typing.overload
    def GetPos(self) -> Point3D:
        """
            Get the location of the default chemical feature (first position)
        
        """
    def GetType(self) -> str:
        """
            Get the specific type for the feature
        
        """
    def SetActiveConformer(self, confId: int) -> None:
        """
            Sets the conformer to use (must be associated with a molecule).
        
        """
class MolChemicalFeatureFactory(Boost.Python.instance):
    """
    Class to featurize a molecule
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
    def GetFeatureDefs(self) -> dict:
        """
            Get a dictionary with SMARTS definitions for each feature type
        
        """
    def GetFeatureFamilies(self) -> tuple:
        """
            Get a tuple of feature types
        
        """
    def GetMolFeature(self, mol: Mol, idx: int, includeOnly: str = '', recompute: bool = True, confId: int = -1) -> MolChemicalFeature:
        """
            returns a particular feature (by index)
        
        """
    def GetNumFeatureDefs(self) -> int:
        """
            Get the number of feature definitions
        
        """
    def GetNumMolFeatures(self, mol: Mol, includeOnly: str = '') -> int:
        """
            Get the number of features the molecule has
        
        """
def BuildFeatureFactory(fileName: str) -> MolChemicalFeatureFactory:
    """
        Construct a feature factory given a feature definition in a file
    
    """
def BuildFeatureFactoryFromString(fdefString: str) -> MolChemicalFeatureFactory:
    """
        Construct a feature factory given a feature definition block
    
    """
def GetAtomMatch(featMatch: typing.Any, maxAts: int = 1024) -> typing.Any:
    """
        Returns an empty list if any of the features passed in share an atom.
         Otherwise a list of lists of atom indices is returned.
        
    
    """
