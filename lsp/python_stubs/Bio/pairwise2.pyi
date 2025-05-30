from .cpairwise2 import rint as rint
from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning, BiopythonWarning as BiopythonWarning
from Bio.Align import substitution_matrices as substitution_matrices
from _typeshed import Incomplete
from typing import NamedTuple

MAX_ALIGNMENTS: int

class Alignment(NamedTuple):
    seqA: Incomplete
    seqB: Incomplete
    score: Incomplete
    start: Incomplete
    end: Incomplete

class align:
    class alignment_function:
        match2args: Incomplete
        penalty2args: Incomplete
        function_name: Incomplete
        align_type: Incomplete
        param_names: Incomplete
        __doc__: Incomplete
        def __init__(self, name) -> None: ...
        def decode(self, *args, **keywds): ...
        def __call__(self, *args, **keywds): ...
    def __getattr__(self, attr): ...

def rint(x, precision=...): ...

class identity_match:
    match: Incomplete
    mismatch: Incomplete
    def __init__(self, match: int = 1, mismatch: int = 0) -> None: ...
    def __call__(self, charA, charB): ...

class dictionary_match:
    score_dict: Incomplete
    symmetric: Incomplete
    def __init__(self, score_dict, symmetric: int = 1) -> None: ...
    def __call__(self, charA, charB): ...

class affine_penalty:
    penalize_extend_when_opening: Incomplete
    def __init__(self, open, extend, penalize_extend_when_opening: int = 0) -> None: ...
    def __call__(self, index, length): ...

def calc_affine_penalty(length, open, extend, penalize_extend_when_opening): ...
def print_matrix(matrix) -> None: ...
def format_alignment(align1, align2, score, begin, end, full_sequences: bool = False): ...
