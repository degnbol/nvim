from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from _typeshed import Incomplete

class DifferentialCutsite:
    start: Incomplete
    enzyme: Incomplete
    cuts_in: Incomplete
    blocked_in: Incomplete
    def __init__(self, **kwds) -> None: ...

class AlignmentHasDifferentLengthsError(Exception): ...

class CAPSMap:
    sequences: Incomplete
    size: Incomplete
    length: Incomplete
    alignment: Incomplete
    enzymes: Incomplete
    def __init__(self, alignment, enzymes: Incomplete | None = None) -> None: ...
