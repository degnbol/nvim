from .Residues import Residues as Residues
from _typeshed import Incomplete
from collections.abc import Generator

class Record:
    sid: str
    residues: Incomplete
    hierarchy: str
    def __init__(self, line: Incomplete | None = None) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
