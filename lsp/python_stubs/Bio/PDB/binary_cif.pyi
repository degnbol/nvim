from Bio import MissingPythonDependencyError as MissingPythonDependencyError
from Bio.PDB.Structure import Structure as Structure
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder

class BinaryCIFParser:
    def __init__(self) -> None: ...
    def get_structure(self, id: str | None, source: str) -> Structure: ...
