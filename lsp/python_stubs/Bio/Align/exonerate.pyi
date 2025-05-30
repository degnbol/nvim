from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    format_alignment: Incomplete
    def __init__(self, target, fmt: str = 'vulgar') -> None: ...
    def write_header(self, stream, alignments) -> None: ...
    def write_footer(self, stream) -> None: ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
