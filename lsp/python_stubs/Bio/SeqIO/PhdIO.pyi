from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter, _TextIOSource
from Bio.SeqRecord import SeqRecord as SeqRecord
from Bio.Sequencing import Phd as Phd
from collections.abc import Iterator as Iterator

class PhdIterator(SequenceIterator):
    modes: str
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...

class PhdWriter(SequenceWriter):
    modes: str
    def write_record(self, record) -> None: ...
