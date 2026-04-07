from Bio.Align import Alignment as Alignment, MultipleSeqAlignment as MultipleSeqAlignment, PairwiseAligner as PairwiseAligner
from Bio.Data import PDBData as PDBData
from Bio.PDB import Selection as Selection
from Bio.PDB.Model import Model as Model
from Bio.PDB.Polypeptide import is_aa as is_aa
from Bio.PDB.Residue import Residue as Residue
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from collections.abc import Generator

class StructureAlignment:
    aligner: Incomplete
    map12: Incomplete
    map21: Incomplete
    duos: Incomplete
    def __init__(self, fasta_align: MultipleSeqAlignment | Alignment | None = None, m1: Model | None = None, m2: Model | None = None, si: int = 0, sj: int = 1, aligner: PairwiseAligner | None = None) -> None: ...
    def get_maps(self): ...
    def get_iterator(self) -> Generator[Incomplete]: ...
