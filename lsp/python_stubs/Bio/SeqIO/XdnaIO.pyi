from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter
from Bio import BiopythonWarning as BiopythonWarning
from Bio.Seq import Seq as Seq
from Bio.SeqFeature import ExactPosition as ExactPosition, SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from Bio.SeqRecord import SeqRecord as SeqRecord

class XdnaIterator(SequenceIterator):
    modes: str
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class XdnaWriter(SequenceWriter):
    modes: str
    def write_records(self, records): ...
