from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Nexus import Nexus as Nexus
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    interleave: Incomplete
    def __init__(self, target, interleave: Incomplete | None = None) -> None: ...
    def write_file(self, stream, alignments): ...
    def format_alignment(self, alignment, interleave: Incomplete | None = None): ...
    def write_alignment(self, alignment, stream, interleave: Incomplete | None = None) -> None: ...
    def write_alignments(self, stream, alignments): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
