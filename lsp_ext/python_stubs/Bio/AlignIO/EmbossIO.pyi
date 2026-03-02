from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.AlignIO.Interfaces import AlignmentIterator as AlignmentIterator
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class EmbossIterator(AlignmentIterator):
    def __next__(self): ...
