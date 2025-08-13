from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq, UndefinedSequenceError as UndefinedSequenceError, reverse_complement as reverse_complement
from Bio.SeqFeature import CompoundLocation as CompoundLocation, ExactPosition as ExactPosition, SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    fmt: str
    header: Incomplete
    wildcard: Incomplete
    mask: Incomplete
    def __init__(self, target, header: bool = True, mask: Incomplete | None = None, wildcard: str = 'N') -> None: ...
    def write_header(self, stream, alignments) -> None: ...
    def format_alignment(self, alignment): ...

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
