from Bio.Align import Alignment as Alignment, Alignments as Alignments, bigbed as bigbed, maf as maf
from Bio.Align.bigbed import AutoSQLTable as AutoSQLTable, Field as Field
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

declaration: Incomplete

class AlignmentWriter(bigbed.AlignmentWriter):
    fmt: str
    def __init__(self, target, targets: Incomplete | None = None, compress: bool = True, blockSize: int = 256, itemsPerSlot: int = 512) -> None: ...
    def write_file(self, stream, alignments): ...

class AlignmentIterator(bigbed.AlignmentIterator, maf.AlignmentIterator):
    fmt: str
    mode: str
    reference: Incomplete
    def __init__(self, source) -> None: ...
