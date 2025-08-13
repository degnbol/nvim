from ._base import _BaseHmmerTextIndexer
from _typeshed import Incomplete

__all__ = ['Hmmer3TextParser', 'Hmmer3TextIndexer']

class Hmmer3TextParser:
    handle: Incomplete
    line: Incomplete
    def __init__(self, handle) -> None: ...
    def __iter__(self): ...

class Hmmer3TextIndexer(_BaseHmmerTextIndexer):
    qresult_start: bytes
    qresult_end: bytes
    def __iter__(self): ...
