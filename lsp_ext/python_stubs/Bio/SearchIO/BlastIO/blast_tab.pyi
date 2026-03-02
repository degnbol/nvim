from Bio.SearchIO._index import SearchIndexer
from _typeshed import Incomplete

__all__ = ['BlastTabIndexer', 'BlastTabParser', 'BlastTabWriter']

class BlastTabParser:
    handle: Incomplete
    has_comments: Incomplete
    fields: Incomplete
    line: Incomplete
    def __init__(self, handle, comments: bool = False, fields=...) -> None: ...
    def __iter__(self): ...

class BlastTabIndexer(SearchIndexer):
    def __init__(self, filename, comments: bool = False, fields=...) -> None: ...
    def __iter__(self): ...
    def get_raw(self, offset): ...

class BlastTabWriter:
    handle: Incomplete
    has_comments: Incomplete
    fields: Incomplete
    def __init__(self, handle, comments: bool = False, fields=...) -> None: ...
    def write_file(self, qresults): ...
