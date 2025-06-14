import abc
from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning, MissingPythonDependencyError as MissingPythonDependencyError
from Bio.Align import _codonaligner, _pairwisealigner, substitution_matrices as substitution_matrices
from Bio.Data import CodonTable as CodonTable
from Bio.Seq import MutableSeq as MutableSeq, Seq as Seq, UndefinedSequenceError as UndefinedSequenceError, reverse_complement as reverse_complement, translate as translate
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from abc import ABC, abstractmethod
from typing import NamedTuple

class AlignmentCounts(NamedTuple):
    gaps: Incomplete
    identities: Incomplete
    mismatches: Incomplete

class MultipleSeqAlignment:
    annotations: Incomplete
    column_annotations: Incomplete
    def __init__(self, records, alphabet: Incomplete | None = None, annotations: Incomplete | None = None, column_annotations: Incomplete | None = None) -> None: ...
    def __format__(self, format_spec) -> str: ...
    def __iter__(self): ...
    def __len__(self) -> int: ...
    def get_alignment_length(self): ...
    def extend(self, records) -> None: ...
    def append(self, record) -> None: ...
    def __add__(self, other): ...
    def __getitem__(self, index): ...
    def __delitem__(self, index) -> None: ...
    def sort(self, key: Incomplete | None = None, reverse: bool = False): ...
    @property
    def substitutions(self): ...
    @property
    def alignment(self): ...

class Alignment:
    @classmethod
    def infer_coordinates(cls, lines): ...
    @classmethod
    def parse_printed_alignment(cls, lines): ...
    sequences: Incomplete
    coordinates: Incomplete
    def __init__(self, sequences, coordinates: Incomplete | None = None) -> None: ...
    def __array__(self, dtype: Incomplete | None = None): ...
    def __add__(self, other): ...
    @property
    def frequencies(self): ...
    @property
    def target(self): ...
    @target.setter
    def target(self, value) -> None: ...
    @property
    def query(self): ...
    @query.setter
    def query(self, value) -> None: ...
    def __eq__(self, other): ...
    def __ne__(self, other): ...
    def __lt__(self, other): ...
    def __le__(self, other): ...
    def __gt__(self, other): ...
    def __ge__(self, other): ...
    def __getitem__(self, key): ...
    def __format__(self, format_spec) -> str: ...
    def format(self, fmt: str = '', *args, **kwargs): ...
    def __len__(self) -> int: ...
    @property
    def length(self): ...
    @property
    def shape(self): ...
    @property
    def aligned(self): ...
    @property
    def indices(self): ...
    @property
    def inverse_indices(self): ...
    def sort(self, key: Incomplete | None = None, reverse: bool = False) -> None: ...
    def map(self, alignment): ...
    def mapall(self, alignments): ...
    @property
    def substitutions(self): ...
    def counts(self): ...
    def reverse_complement(self): ...

class AlignmentsAbstractBaseClass(ABC, metaclass=abc.ABCMeta):
    def __iter__(self): ...
    @abstractmethod
    def __next__(self): ...
    @abstractmethod
    def rewind(self): ...
    @abstractmethod
    def __len__(self): ...

class Alignments(AlignmentsAbstractBaseClass, list):
    def __init__(self, alignments=()) -> None: ...
    def __next__(self): ...
    def rewind(self) -> None: ...
    def __len__(self) -> int: ...

class PairwiseAlignments(AlignmentsAbstractBaseClass):
    sequences: Incomplete
    score: Incomplete
    def __init__(self, seqA, seqB, score, paths) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index): ...
    def __next__(self): ...
    def rewind(self) -> None: ...

class PairwiseAligner(_pairwisealigner.PairwiseAligner):
    substitution_matrix: Incomplete
    open_gap_score: float
    extend_gap_score: float
    def __init__(self, scoring: Incomplete | None = None, **kwargs) -> None: ...
    def __setattr__(self, key, value) -> None: ...
    def align(self, seqA, seqB, strand: str = '+'): ...
    def score(self, seqA, seqB, strand: str = '+'): ...

class CodonAligner(_codonaligner.CodonAligner):
    codon_table: Incomplete
    def __init__(self, codon_table: Incomplete | None = None, anchor_len: int = 10) -> None: ...
    def score(self, seqA, seqB): ...
    def align(self, seqA, seqB): ...

formats: Incomplete

def write(alignments, target, fmt, *args, **kwargs): ...
def parse(source, fmt): ...
def read(handle, fmt): ...
