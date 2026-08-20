# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from rdkit.ForceField.rdForceField import ForceField
from rdkit.ForceField.rdForceField import MMFFMolProperties
from .rdForceField import *
__all__: list[str] = ['ForceField', 'MMFFMolProperties', 'rdForceField']

# present at runtime, absent from the generated stub:
from . import rdForceField as rdForceField
