from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    def write_header(self, stream, alignments) -> None: ...
    def format_alignment(self, alignment): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
