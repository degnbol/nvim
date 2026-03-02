from Bio.Data import PDBData as PDBData
from Bio.PDB import Selection as Selection
from Bio.PDB.Polypeptide import is_aa as is_aa
from _typeshed import Incomplete
from collections.abc import Generator

class StructureAlignment:
    map12: Incomplete
    map21: Incomplete
    duos: Incomplete
    def __init__(self, fasta_align, m1, m2, si: int = 0, sj: int = 1) -> None: ...
    def get_maps(self): ...
    def get_iterator(self) -> Generator[Incomplete]: ...
