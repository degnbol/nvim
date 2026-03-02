from Bio.PDB.PDBIO import Select as Select, StructureIO as StructureIO
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder
from _typeshed import Incomplete

mmcif_order: Incomplete

class MMCIFIO(StructureIO):
    def __init__(self) -> None: ...
    dic: Incomplete
    def set_dict(self, dic) -> None: ...
    def save(self, filepath, select=..., preserve_atom_numbering: bool = False) -> None: ...
