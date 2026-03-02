from Bio.motifs import meme as meme
from _typeshed import Incomplete

class Record(list):
    sequences: Incomplete
    version: str
    database: str
    diagrams: Incomplete
    alphabet: Incomplete
    strand_handling: str
    def __init__(self) -> None: ...
    def __getitem__(self, key): ...

def read(handle): ...
