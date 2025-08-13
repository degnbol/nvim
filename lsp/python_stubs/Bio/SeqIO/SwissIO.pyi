from .Interfaces import SequenceIterator as SequenceIterator, _TextIOSource
from Bio import SeqFeature as SeqFeature, SwissProt as SwissProt
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class SwissIterator(SequenceIterator):
    modes: str
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...
