from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentWriter as AlignmentWriter
from Bio.Nexus import Nexus as Nexus
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from collections.abc import Iterator
from typing import IO

def NexusIterator(handle: IO[str], seq_count: int | None = None) -> Iterator[MultipleSeqAlignment]: ...

class NexusWriter(AlignmentWriter):
    def write_file(self, alignments): ...
    def write_alignment(self, alignment, interleave: Incomplete | None = None) -> None: ...
