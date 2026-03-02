from Bio.File import as_handle as as_handle
from Bio.PDB.PDBExceptions import PDBConstructionException as PDBConstructionException, PDBConstructionWarning as PDBConstructionWarning
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder
from _typeshed import Incomplete

class PDBParser:
    structure_builder: Incomplete
    header: Incomplete
    trailer: Incomplete
    line_counter: int
    PERMISSIVE: Incomplete
    QUIET: Incomplete
    is_pqr: Incomplete
    def __init__(self, PERMISSIVE: bool = True, get_header: bool = False, structure_builder: Incomplete | None = None, QUIET: bool = False, is_pqr: bool = False) -> None: ...
    def get_structure(self, id, file): ...
    def get_header(self): ...
    def get_trailer(self): ...
