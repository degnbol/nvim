from .Interfaces import SequenceIterator as SequenceIterator
from Bio.Seq import Seq as Seq
from Bio.SeqFeature import SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from Bio.SeqRecord import SeqRecord as SeqRecord

class GckIterator(SequenceIterator):
    modes: str
    def __init__(self, source) -> None: ...
    def __next__(self): ...
