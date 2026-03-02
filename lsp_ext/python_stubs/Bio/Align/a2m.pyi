from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    def format_alignment(self, alignment): ...
    write_alignments: Incomplete

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
