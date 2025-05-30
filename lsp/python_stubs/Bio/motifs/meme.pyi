from Bio import Align as Align, BiopythonDeprecationWarning as BiopythonDeprecationWarning, motifs as motifs
from Bio.Seq import Seq as Seq
from _typeshed import Incomplete

def read(handle): ...

class Motif(motifs.Motif):
    evalue: float
    num_occurrences: int
    name: Incomplete
    id: Incomplete
    alt_id: Incomplete
    def __init__(self, alphabet: Incomplete | None = None, alignment: Incomplete | None = None) -> None: ...

class Instance(Seq):
    sequence_name: str
    sequence_id: str
    start: int
    pvalue: float
    strand: int
    length: int
    motif_name: str
    def __init__(self, *args, **kwds) -> None: ...

class Record(list):
    version: str
    datafile: str
    command: str
    alphabet: str
    sequences: Incomplete
    def __init__(self) -> None: ...
    def __getitem__(self, key): ...
