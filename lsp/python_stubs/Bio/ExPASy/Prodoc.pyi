from _typeshed import Incomplete
from collections.abc import Generator

def read(handle): ...
def parse(handle) -> Generator[Incomplete]: ...

class Record:
    accession: str
    prosite_refs: Incomplete
    text: str
    references: Incomplete
    def __init__(self) -> None: ...

class Reference:
    number: str
    authors: str
    citation: str
    def __init__(self) -> None: ...
