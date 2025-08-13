from Bio import BiopythonParserWarning as BiopythonParserWarning
from Bio.File import as_handle as as_handle
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from collections.abc import Generator

class InsdcScanner:
    RECORD_START: str
    HEADER_WIDTH: int
    FEATURE_START_MARKERS: Incomplete
    FEATURE_END_MARKERS: Incomplete
    FEATURE_QUALIFIER_INDENT: int
    FEATURE_QUALIFIER_SPACER: str
    SEQUENCE_HEADERS: Incomplete
    debug: Incomplete
    handle: Incomplete
    line: Incomplete
    def __init__(self, debug: int = 0) -> None: ...
    def set_handle(self, handle) -> None: ...
    def find_start(self): ...
    def parse_header(self): ...
    def parse_features(self, skip: bool = False): ...
    def parse_feature(self, feature_key, lines): ...
    def parse_footer(self): ...
    def feed(self, handle, consumer, do_features: bool = True): ...
    def parse(self, handle, do_features: bool = True): ...
    def parse_records(self, handle, do_features: bool = True) -> Generator[Incomplete]: ...
    def parse_cds_features(self, handle, alphabet: Incomplete | None = None, tags2id=('protein_id', 'locus_tag', 'product')) -> Generator[Incomplete]: ...

class EmblScanner(InsdcScanner):
    RECORD_START: str
    HEADER_WIDTH: int
    FEATURE_START_MARKERS: Incomplete
    FEATURE_END_MARKERS: Incomplete
    FEATURE_QUALIFIER_INDENT: int
    FEATURE_QUALIFIER_SPACER: Incomplete
    SEQUENCE_HEADERS: Incomplete
    EMBL_INDENT = HEADER_WIDTH
    EMBL_SPACER: Incomplete
    line: Incomplete
    def parse_footer(self): ...

class _ImgtScanner(EmblScanner):
    FEATURE_START_MARKERS: Incomplete
    line: Incomplete
    def parse_features(self, skip: bool = False): ...

class GenBankScanner(InsdcScanner):
    RECORD_START: str
    HEADER_WIDTH: int
    FEATURE_START_MARKERS: Incomplete
    FEATURE_END_MARKERS: list[str]
    FEATURE_QUALIFIER_INDENT: int
    FEATURE_QUALIFIER_SPACER: Incomplete
    SEQUENCE_HEADERS: Incomplete
    GENBANK_INDENT = HEADER_WIDTH
    GENBANK_SPACER: Incomplete
    STRUCTURED_COMMENT_START: str
    STRUCTURED_COMMENT_END: str
    STRUCTURED_COMMENT_DELIM: str
    line: Incomplete
    def parse_footer(self): ...
