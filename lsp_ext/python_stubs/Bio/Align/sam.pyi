from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq, UndefinedSequenceError as UndefinedSequenceError, reverse_complement as reverse_complement
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    md: Incomplete
    def __init__(self, target, md: bool = False) -> None: ...
    def write_header(self, stream, alignments) -> None: ...
    def format_alignment(self, alignment, md: Incomplete | None = None): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
