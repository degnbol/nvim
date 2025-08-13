from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.ccealign import run_cealign as run_cealign
from Bio.PDB.qcprot import QCPSuperimposer as QCPSuperimposer
from _typeshed import Incomplete

class CEAligner:
    window_size: Incomplete
    max_gap: Incomplete
    rms: Incomplete
    refcoord: Incomplete
    def __init__(self, window_size: int = 8, max_gap: int = 30) -> None: ...
    def get_guide_coord_from_structure(self, structure): ...
    def set_reference(self, structure) -> None: ...
    def align(self, structure, transform: bool = True, *, final_optimization: bool = True) -> None: ...
