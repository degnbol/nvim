from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter
from Bio import StreamModeError as StreamModeError
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

def ReadRocheXmlManifest(handle): ...

class SffIterator(SequenceIterator):
    modes: str
    read_header_fmt: str
    read_header_size: Incomplete
    trim: Incomplete
    read_flow_fmt: Incomplete
    read_flow_size: Incomplete
    def __init__(self, source, alphabet: Incomplete | None = None, trim: bool = False) -> None: ...
    def __next__(self): ...

class _SffTrimIterator(SffIterator):
    def __init__(self, source) -> None: ...

class SffWriter(SequenceWriter):
    modes: str
    def __init__(self, target, index: bool = True, xml: Incomplete | None = None) -> None: ...
    def write_records(self, records): ...
    def write_header(self) -> None: ...
    def write_record(self, record) -> None: ...
