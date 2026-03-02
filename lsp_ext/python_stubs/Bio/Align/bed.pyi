from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    bedN: Incomplete
    def __init__(self, target, bedN: int = 12) -> None: ...
    def format_alignment(self, alignment): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
