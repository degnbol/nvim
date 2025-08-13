from Bio.File import as_handle as as_handle
from Bio.PDB.MMCIF2Dict import MMCIF2Dict as MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBConstructionException as PDBConstructionException, PDBConstructionWarning as PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder
from _typeshed import Incomplete

class MMCIFParser:
    header: Incomplete
    line_counter: int
    build_structure: Incomplete
    auth_chains: Incomplete
    auth_residues: Incomplete
    QUIET: Incomplete
    def __init__(self, structure_builder: Incomplete | None = None, auth_chains: bool = True, auth_residues: bool = True, QUIET: bool = False) -> None: ...
    def get_structure(self, structure_id, filename): ...

class FastMMCIFParser:
    line_counter: int
    build_structure: Incomplete
    auth_chains: Incomplete
    auth_residues: Incomplete
    QUIET: Incomplete
    def __init__(self, structure_builder: Incomplete | None = None, auth_chains: bool = True, auth_residues: bool = True, QUIET: bool = False) -> None: ...
    def get_structure(self, structure_id, filename): ...
