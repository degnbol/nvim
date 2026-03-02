from _typeshed import Incomplete
from collections.abc import Generator

def read(handle): ...
def parse(handle) -> Generator[Incomplete]: ...

class Record:
    query: str
    hit: str
    gap_threshold: int
    query_length: int
    query_filtered_length: int
    query_nseqs: int
    query_neffseqs: int
    hit_length: int
    hit_filtered_length: int
    hit_nseqs: int
    hit_neffseqs: int
    sw_score: int
    evalue: int
    query_start: int
    hit_start: int
    query_aln: str
    hit_aln: str
    positives: str
    def __init__(self) -> None: ...
    def query_coverage(self): ...
    def hit_coverage(self): ...
