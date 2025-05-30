from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.SVDSuperimposer import SVDSuperimposer as SVDSuperimposer
from _typeshed import Incomplete

class Superimposer:
    rotran: Incomplete
    rms: Incomplete
    def __init__(self) -> None: ...
    def set_atoms(self, fixed, moving) -> None: ...
    def apply(self, atom_list) -> None: ...
