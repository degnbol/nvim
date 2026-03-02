from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class PirIterator(SequenceIterator):
    modes: str
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class PirWriter(SequenceWriter):
    modes: str
    wrap: Incomplete
    record2title: Incomplete
    code: Incomplete
    def __init__(self, handle, wrap: int = 60, record2title: Incomplete | None = None, code: Incomplete | None = None) -> None: ...
    def write_record(self, record) -> None: ...
