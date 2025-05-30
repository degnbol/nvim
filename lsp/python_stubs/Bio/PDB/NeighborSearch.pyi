from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.Selection import entity_levels as entity_levels, unfold_entities as unfold_entities, uniqueify as uniqueify
from _typeshed import Incomplete

class NeighborSearch:
    atom_list: Incomplete
    coords: Incomplete
    kdt: Incomplete
    def __init__(self, atom_list, bucket_size: int = 10) -> None: ...
    def search(self, center, radius, level: str = 'A'): ...
    def search_all(self, radius, level: str = 'A'): ...
