from ._base import _BaseInfernalParser
from Bio.SearchIO._index import SearchIndexer
from _typeshed import Incomplete

__all__ = ['InfernalTabParser', 'InfernalTabIndexer']

class InfernalTabParser(_BaseInfernalParser):
    handle: Incomplete
    fmt: Incomplete
    def __init__(self, handle, _fmt=None) -> None: ...
    line: Incomplete
    def __iter__(self): ...

class InfernalTabIndexer(SearchIndexer):
    def __init__(self, filename) -> None: ...
    def __iter__(self): ...
    def get_raw(self, offset): ...
