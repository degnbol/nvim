from Bio import motifs as motifs
from _typeshed import Incomplete

def read(handle): ...

class Record(list):
    version: str
    datafile: str
    command: str
    alphabet: Incomplete
    background: Incomplete
    sequences: Incomplete
    def __init__(self) -> None: ...
    def __getitem__(self, key): ...
