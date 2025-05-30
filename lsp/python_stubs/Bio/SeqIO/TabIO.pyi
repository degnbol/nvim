from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class TabIterator(SequenceIterator):
    modes: str
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class TabWriter(SequenceWriter):
    modes: str
    def write_record(self, record) -> None: ...

def as_tab(record): ...
