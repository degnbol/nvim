from _typeshed import Incomplete
from collections.abc import Generator

def get_indiv(line): ...
def read(handle): ...

class Record:
    handle: Incomplete
    marker_len: int
    comment_line: str
    loci_list: Incomplete
    populations: Incomplete
    stack: Incomplete
    def __init__(self, handle) -> None: ...
    def data_generator(self) -> Generator[Incomplete]: ...
