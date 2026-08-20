# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from rdkit.Chem.rdMolCatalog import MolCatalog
from rdkit.Chem.rdMolCatalog import MolCatalogEntry
__all__: list[str] = ['MolCatalog', 'MolCatalogEntry']

# present at runtime, absent from the generated stub:
from rdkit.Chem.rdMolCatalog import CreateMolCatalog as CreateMolCatalog
