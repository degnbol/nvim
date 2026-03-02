from _typeshed import Incomplete
from collections.abc import Generator

name_wrap: Incomplete
id_wrap: Incomplete

class Record:
    entry: str
    name: Incomplete
    definition: str
    orthology: Incomplete
    organism: str
    position: str
    motif: Incomplete
    dblinks: Incomplete
    def __init__(self) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
