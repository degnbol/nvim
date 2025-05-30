from Bio import MissingPythonDependencyError as MissingPythonDependencyError
from _typeshed import Incomplete

class ParserError(ValueError):
    def __init__(self, *args) -> None: ...

class Record:
    version: Incomplete
    GridCornerUL: Incomplete
    GridCornerUR: Incomplete
    GridCornerLR: Incomplete
    GridCornerLL: Incomplete
    DatHeader: Incomplete
    Algorithm: Incomplete
    AlgorithmParameters: Incomplete
    NumberCells: Incomplete
    intensities: Incomplete
    stdevs: Incomplete
    npix: Incomplete
    nrows: Incomplete
    ncols: Incomplete
    nmask: Incomplete
    mask: Incomplete
    noutliers: Incomplete
    outliers: Incomplete
    modified: Incomplete
    def __init__(self) -> None: ...

def read(handle, version: Incomplete | None = None): ...
