from _typeshed import Incomplete
from collections.abc import Generator

class Record:
    sunid: str
    nodetype: str
    sccs: str
    name: str
    description: str
    def __init__(self, line: Incomplete | None = None) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
