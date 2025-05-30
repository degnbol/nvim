from Bio.Align import Alignment as Alignment
from Bio.Seq import Seq as Seq
from Bio.motifs import Motif as Motif
from _typeshed import Incomplete

class Record(list):
    parameters: Incomplete
    def __init__(self) -> None: ...

def read(handle): ...
