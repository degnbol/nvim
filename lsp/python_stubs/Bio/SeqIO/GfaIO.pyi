from .Interfaces import SequenceIterator as SequenceIterator, _TextIOSource
from Bio import BiopythonWarning as BiopythonWarning
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord

class Gfa1Iterator(SequenceIterator):
    modes: str
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...

class Gfa2Iterator(SequenceIterator):
    modes: str
    def __init__(self, source: _TextIOSource) -> None: ...
    def __next__(self): ...
