from Bio import Align as Align, motifs as motifs
from Bio.Seq import Seq as Seq
from _typeshed import Incomplete

def read(handle): ...

class Motif(motifs.Motif):
    evalue: float
    num_occurrences: int
    name: Incomplete
    id: Incomplete
    alt_id: Incomplete
    def __init__(self, alphabet=None, alignment=None) -> None: ...

class Record(list):
    version: str
    datafile: str
    command: str
    alphabet: str
    sequences: Incomplete
    def __init__(self) -> None: ...
    def __getitem__(self, key): ...
