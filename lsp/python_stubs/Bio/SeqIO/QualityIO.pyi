import abc
from .Interfaces import SequenceIterator as SequenceIterator, SequenceWriter as SequenceWriter, _TextIOSource
from Bio import BiopythonParserWarning as BiopythonParserWarning, BiopythonWarning as BiopythonWarning, StreamModeError as StreamModeError
from Bio.File import as_handle as as_handle
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from abc import abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Callable

SANGER_SCORE_OFFSET: int
SOLEXA_SCORE_OFFSET: int
INVALID_CHAR_CODE: int
INVALID_CHAR: Incomplete

def solexa_quality_from_phred(phred_quality: float) -> float: ...
def phred_quality_from_solexa(solexa_quality: float) -> float: ...
def FastqGeneralIterator(source: _TextIOSource) -> Iterator[tuple[str, str, str]]: ...

class FastqIteratorAbstractBaseClass(SequenceIterator[str], metaclass=abc.ABCMeta):
    modes: str
    @property
    @abstractmethod
    def q_mapping(self): ...
    @property
    @abstractmethod
    def q_key(self): ...
    line: Incomplete
    def __init__(self, source) -> None: ...
    def __next__(self) -> SeqRecord: ...

class FastqPhredIterator(FastqIteratorAbstractBaseClass):
    q_mapping: Incomplete
    q_key: str
    def __init__(self, source: _TextIOSource, alphabet: None = None) -> None: ...

class FastqSolexaIterator(FastqIteratorAbstractBaseClass):
    q_mapping: Incomplete
    q_key: str
    def __init__(self, source: _TextIOSource, alphabet: None = None) -> None: ...

class FastqIlluminaIterator(FastqIteratorAbstractBaseClass):
    q_mapping: Incomplete
    q_key: str
    def __init__(self, source: _TextIOSource, alphabet: None = None) -> None: ...

class QualPhredIterator(SequenceIterator):
    modes: str
    def __init__(self, source: _TextIOSource, alphabet: None = None) -> None: ...
    def __next__(self) -> SeqRecord: ...

class FastqPhredWriter(SequenceWriter):
    modes: str
    def write_record(self, record: SeqRecord) -> None: ...

def as_fastq(record: SeqRecord) -> str: ...

class QualPhredWriter(SequenceWriter):
    modes: str
    wrap: int | None
    record2title: Incomplete
    def __init__(self, handle: _TextIOSource, wrap: int = 60, record2title: Callable[[SeqRecord], str] | None = None) -> None: ...
    def write_record(self, record: SeqRecord) -> None: ...

def as_qual(record: SeqRecord) -> str: ...

class FastqSolexaWriter(SequenceWriter):
    modes: str
    def write_record(self, record: SeqRecord) -> None: ...

def as_fastq_solexa(record: SeqRecord) -> str: ...

class FastqIlluminaWriter(SequenceWriter):
    modes: str
    def write_record(self, record: SeqRecord) -> None: ...

def as_fastq_illumina(record: SeqRecord) -> str: ...
def PairedFastaQualIterator(fasta_source: _TextIOSource, qual_source: _TextIOSource, alphabet: None = None) -> Iterator[SeqRecord]: ...

@dataclass
class InvalidCharError(ValueError):
    full_string: str
    index: int
    details: str
    r: int = ...
