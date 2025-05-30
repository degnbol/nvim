from .Interfaces import AlignmentIterator as AlignmentIterator
from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class MsfIterator(AlignmentIterator):
    def __next__(self): ...
