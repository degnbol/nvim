from .Interfaces import SequenceIterator as SequenceIterator, _TextIOSource
from Bio import BiopythonParserWarning as BiopythonParserWarning
from Bio.Data.PDBData import protein_letters_3to1 as protein_letters_3to1, protein_letters_3to1_extended as protein_letters_3to1_extended
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from collections.abc import Generator

def AtomIterator(pdb_id, structure) -> Generator[Incomplete]: ...

class PdbSeqresIterator(SequenceIterator):
    modes: str
    cache: Incomplete
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...

class PdbAtomIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...

PDBX_POLY_SEQ_SCHEME_FIELDS: Incomplete
STRUCT_REF_FIELDS: Incomplete
STRUCT_REF_SEQ_FIELDS: Incomplete

class CifSeqresIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...

class CifAtomIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...
