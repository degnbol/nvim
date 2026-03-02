from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator as AlignmentIterator, SequentialAlignmentWriter as SequentialAlignmentWriter
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class ClustalWriter(SequentialAlignmentWriter):
    def write_alignment(self, alignment) -> None: ...

class ClustalIterator(AlignmentIterator):
    def __next__(self): ...
