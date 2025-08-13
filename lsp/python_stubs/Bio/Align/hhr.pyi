from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
    def __len__(self) -> int: ...
