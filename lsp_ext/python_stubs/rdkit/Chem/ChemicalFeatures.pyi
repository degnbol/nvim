# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from rdkit.Chem.rdChemicalFeatures import FreeChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeature
from rdkit.Chem.rdMolChemicalFeatures import MolChemicalFeatureFactory
__all__: list[str] = ['FreeChemicalFeature', 'MCFF_GetFeaturesForMol', 'MolChemicalFeature', 'MolChemicalFeatureFactory']
def MCFF_GetFeaturesForMol(self, mol, includeOnly = '', confId = -1):
    ...

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdMolChemicalFeatures import BuildFeatureFactory as BuildFeatureFactory
from rdkit.Chem.rdMolChemicalFeatures import BuildFeatureFactoryFromString as BuildFeatureFactoryFromString
from rdkit.Chem.rdMolChemicalFeatures import GetAtomMatch as GetAtomMatch
