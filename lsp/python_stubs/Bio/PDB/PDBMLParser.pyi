from Bio.PDB.Structure import Structure as Structure
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder
from _typeshed import Incomplete
from os import PathLike
from typing import TextIO

class PDBMLParser:
    structure_builder: Incomplete
    def __init__(self) -> None: ...
    def get_structure(self, source: int | str | bytes | PathLike | TextIO) -> Structure: ...
