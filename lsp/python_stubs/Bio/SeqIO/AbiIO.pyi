from .Interfaces import SequenceIterator as SequenceIterator
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AbiIterator(SequenceIterator):
    modes: str
    trim: Incomplete
    def __init__(self, source, trim: bool = False) -> None: ...
    def __next__(self): ...
