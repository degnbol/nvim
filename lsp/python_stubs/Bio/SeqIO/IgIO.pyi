from .Interfaces import SequenceIterator as SequenceIterator
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class IgIterator(SequenceIterator):
    modes: str
    def __init__(self, source) -> None: ...
    def __next__(self): ...
