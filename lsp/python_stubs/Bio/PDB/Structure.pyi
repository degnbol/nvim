from Bio.PDB.Entity import Entity as Entity
from Bio.PDB.Model import Model as Model
from _typeshed import Incomplete
from collections.abc import Generator

class Structure(Entity[None, 'Model']):
    level: str
    def __init__(self, id) -> None: ...
    def get_models(self) -> Generator[Incomplete, Incomplete]: ...
    def get_chains(self) -> Generator[Incomplete, Incomplete]: ...
    def get_residues(self) -> Generator[Incomplete, Incomplete]: ...
    def get_atoms(self) -> Generator[Incomplete, Incomplete]: ...
    def atom_to_internal_coordinates(self, verbose: bool = False) -> None: ...
    def internal_to_atom_coordinates(self, verbose: bool = False) -> None: ...
