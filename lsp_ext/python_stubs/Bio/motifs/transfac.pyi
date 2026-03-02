from Bio import motifs as motifs
from _typeshed import Incomplete

class Motif(motifs.Motif, dict):
    multiple_value_keys: Incomplete
    reference_keys: Incomplete
    def __getitem__(self, key): ...

class Record(list):
    version: Incomplete
    def __init__(self) -> None: ...

def read(handle, strict: bool = True): ...
def write(motifs): ...
