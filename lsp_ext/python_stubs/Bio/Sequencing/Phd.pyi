from Bio import Seq as Seq
from _typeshed import Incomplete
from collections.abc import Generator

CKEYWORDS: Incomplete

class Record:
    file_name: str
    comments: Incomplete
    sites: Incomplete
    seq: str
    seq_trimmed: str
    def __init__(self) -> None: ...

def read(source): ...
def parse(source) -> Generator[Incomplete]: ...
