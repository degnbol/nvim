from Bio.Align import Alignment as Alignment, Alignments as Alignments, bigbed as bigbed
from Bio.Align.bigbed import AutoSQLTable as AutoSQLTable, Field as Field
from Bio.Seq import Seq as Seq, UndefinedSequenceError as UndefinedSequenceError, reverse_complement as reverse_complement
from Bio.SeqFeature import Location as Location, SeqFeature as SeqFeature
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

declaration: Incomplete

class AlignmentWriter(bigbed.AlignmentWriter):
    fmt: str
    cds: Incomplete
    fa: Incomplete
    mask: Incomplete
    wildcard: Incomplete
    def __init__(self, target, targets: Incomplete | None = None, compress: bool = True, extraIndex=(), cds: bool = False, fa: bool = False, mask: Incomplete | None = None, wildcard: str = 'N') -> None: ...
    def write_file(self, stream, alignments): ...

class AlignmentIterator(bigbed.AlignmentIterator):
    fmt: str
