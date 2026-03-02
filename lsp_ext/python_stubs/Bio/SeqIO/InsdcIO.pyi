from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter
from Bio import BiopythonWarning as BiopythonWarning, SeqFeature as SeqFeature, SeqIO as SeqIO
from Bio.GenBank.Scanner import EmblScanner as EmblScanner, GenBankScanner as GenBankScanner
from Bio.Seq import UndefinedSequenceError as UndefinedSequenceError
from _typeshed import Incomplete

class GenBankIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class EmblIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class ImgtIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class GenBankCdsFeatureIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class EmblCdsFeatureIterator(SequenceIterator):
    modes: str
    records: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self): ...

class _InsdcWriter(SequenceWriter):
    MAX_WIDTH: int
    QUALIFIER_INDENT: int
    QUALIFIER_INDENT_STR: Incomplete
    QUALIFIER_INDENT_TMP: str
    FTQUAL_NO_QUOTE: Incomplete
    modes: str

class GenBankWriter(_InsdcWriter):
    HEADER_WIDTH: int
    QUALIFIER_INDENT: int
    STRUCTURED_COMMENT_START: str
    STRUCTURED_COMMENT_END: str
    STRUCTURED_COMMENT_DELIM: str
    LETTERS_PER_LINE: int
    SEQUENCE_INDENT: int
    def write_record(self, record) -> None: ...

class EmblWriter(_InsdcWriter):
    HEADER_WIDTH: int
    QUALIFIER_INDENT: int
    QUALIFIER_INDENT_STR: Incomplete
    QUALIFIER_INDENT_TMP: str
    FEATURE_HEADER: str
    LETTERS_PER_BLOCK: int
    BLOCKS_PER_LINE: int
    LETTERS_PER_LINE: Incomplete
    POSITION_PADDING: int
    def write_record(self, record) -> None: ...

class ImgtWriter(EmblWriter):
    HEADER_WIDTH: int
    QUALIFIER_INDENT: int
    QUALIFIER_INDENT_STR: Incomplete
    QUALIFIER_INDENT_TMP: str
    FEATURE_HEADER: str
