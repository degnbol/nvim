from _typeshed import Incomplete
from collections.abc import Generator

name_wrap: Incomplete
id_wrap: Incomplete
struct_wrap: Incomplete

class Record:
    entry: str
    name: Incomplete
    formula: str
    mass: str
    pathway: Incomplete
    enzyme: Incomplete
    structures: Incomplete
    dblinks: Incomplete
    def __init__(self) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
