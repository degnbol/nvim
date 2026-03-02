from .Interfaces import AlignmentIterator as AlignmentIterator, SequentialAlignmentWriter as SequentialAlignmentWriter
from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class PhylipWriter(SequentialAlignmentWriter):
    def write_alignment(self, alignment, id_width=...) -> None: ...

class PhylipIterator(AlignmentIterator):
    id_width: Incomplete
    def __next__(self): ...

class RelaxedPhylipWriter(PhylipWriter):
    def write_alignment(self, alignment) -> None: ...

class RelaxedPhylipIterator(PhylipIterator): ...

class SequentialPhylipWriter(SequentialAlignmentWriter):
    def write_alignment(self, alignment, id_width=...) -> None: ...

class SequentialPhylipIterator(PhylipIterator):
    def __next__(self): ...

def sanitize_name(name, width: Incomplete | None = None): ...
