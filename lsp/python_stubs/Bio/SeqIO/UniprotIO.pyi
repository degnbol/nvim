from .Interfaces import SequenceIterator as SequenceIterator, _BytesIOSource
from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning, SeqFeature as SeqFeature
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

NS: str
REFERENCE_JOURNAL: str

class UniprotIterator(SequenceIterator):
    modes: str
    return_raw_comments: Incomplete
    def __init__(self, source: _BytesIOSource, alphabet: None = None, return_raw_comments: bool = False) -> None: ...
    def __next__(self): ...
