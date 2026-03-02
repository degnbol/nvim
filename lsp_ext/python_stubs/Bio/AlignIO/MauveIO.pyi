from .Interfaces import AlignmentIterator as AlignmentIterator, SequentialAlignmentWriter as SequentialAlignmentWriter
from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

XMFA_HEADER_REGEX: Incomplete
XMFA_HEADER_REGEX_BIOPYTHON: Incomplete
ID_LINE_FMT: str

class MauveWriter(SequentialAlignmentWriter):
    def __init__(self, *args, **kwargs) -> None: ...
    def write_alignment(self, alignment) -> None: ...

class MauveIterator(AlignmentIterator):
    ids: Incomplete
    sequences: Incomplete
    def __next__(self): ...
