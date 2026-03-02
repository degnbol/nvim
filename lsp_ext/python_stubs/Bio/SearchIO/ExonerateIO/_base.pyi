import abc
from Bio.SearchIO._index import SearchIndexer as SearchIndexer
from Bio.SearchIO._model import HSP as HSP, HSPFragment as HSPFragment, Hit as Hit, QueryResult as QueryResult
from Bio.SeqUtils import seq1 as seq1
from _typeshed import Incomplete
from abc import ABC, abstractmethod

class _BaseExonerateParser(ABC, metaclass=abc.ABCMeta):
    handle: Incomplete
    has_c4_alignment: bool
    def __init__(self, handle) -> None: ...
    line: Incomplete
    def __iter__(self): ...
    def read_until(self, bool_func) -> None: ...
    @abstractmethod
    def parse_alignment_block(self, header): ...

class _BaseExonerateIndexer(SearchIndexer):
    def get_qresult_id(self, pos) -> None: ...
    def __iter__(self): ...
