from . import Residues as Residues
from _typeshed import Incomplete
from collections.abc import Generator

class Record:
    sid: str
    residues: Incomplete
    sccs: str
    sunid: str
    hierarchy: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...

class Index(dict):
    filename: Incomplete
    def __init__(self, filename) -> None: ...
    def __getitem__(self, key): ...
