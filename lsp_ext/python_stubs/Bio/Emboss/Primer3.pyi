from _typeshed import Incomplete
from collections.abc import Generator

class Record:
    comments: str
    primers: Incomplete
    def __init__(self) -> None: ...

class Primers:
    size: int
    forward_seq: str
    forward_start: int
    forward_length: int
    forward_tm: float
    forward_gc: float
    reverse_seq: str
    reverse_start: int
    reverse_length: int
    reverse_tm: float
    reverse_gc: float
    internal_seq: str
    internal_start: int
    internal_length: int
    internal_tm: float
    internal_gc: float
    def __init__(self) -> None: ...
    def __len__(self) -> int: ...

def parse(handle) -> Generator[Incomplete]: ...
def read(handle): ...
