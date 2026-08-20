# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from rdkit.ForceField import rdForceField
from rdkit.ForceField.rdForceField import ForceField
from rdkit.ForceField.rdForceField import MMFFMolProperties
__all__: list[str] = ['ForceField', 'MMFFMolProperties', 'rdForceField']

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdForceFieldHelpers import CreateEmptyForceFieldForMol as CreateEmptyForceFieldForMol
from rdkit.Chem.rdForceFieldHelpers import GetUFFAngleBendParams as GetUFFAngleBendParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFBondStretchParams as GetUFFBondStretchParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFInversionParams as GetUFFInversionParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFTorsionParams as GetUFFTorsionParams
from rdkit.Chem.rdForceFieldHelpers import GetUFFVdWParams as GetUFFVdWParams
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeForceField as MMFFGetMoleculeForceField
from rdkit.Chem.rdForceFieldHelpers import MMFFGetMoleculeProperties as MMFFGetMoleculeProperties
from rdkit.Chem.rdForceFieldHelpers import MMFFHasAllMoleculeParams as MMFFHasAllMoleculeParams
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule as MMFFOptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs as MMFFOptimizeMoleculeConfs
from rdkit.Chem.rdForceFieldHelpers import MMFFSanitizeMolecule as MMFFSanitizeMolecule
from rdkit.Chem.rdForceFieldHelpers import OptimizeMolecule as OptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import OptimizeMoleculeConfs as OptimizeMoleculeConfs
from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField as UFFGetMoleculeForceField
from rdkit.Chem.rdForceFieldHelpers import UFFHasAllMoleculeParams as UFFHasAllMoleculeParams
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule as UFFOptimizeMolecule
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMoleculeConfs as UFFOptimizeMoleculeConfs
