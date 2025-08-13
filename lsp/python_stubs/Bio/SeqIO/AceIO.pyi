from .Interfaces import SequenceIterator as SequenceIterator, _TextIOSource
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from Bio.Sequencing import Ace as Ace
from _typeshed import Incomplete

class AceIterator(SequenceIterator):
    modes: str
    ace_contigs: Incomplete
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...
