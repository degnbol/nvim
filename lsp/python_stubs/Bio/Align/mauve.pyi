from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    def __init__(self, target, metadata: Incomplete | None = None, identifiers: Incomplete | None = None) -> None: ...
    def write_header(self, stream, alignments) -> None: ...
    def write_file(self, stream, alignments): ...
    def format_alignment(self, alignment): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
