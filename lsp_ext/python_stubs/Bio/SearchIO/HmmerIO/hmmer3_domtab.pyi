from .hmmer3_tab import Hmmer3TabIndexer, Hmmer3TabParser
from _typeshed import Incomplete

__all__ = ['Hmmer3DomtabHmmhitParser', 'Hmmer3DomtabHmmqueryParser', 'Hmmer3DomtabHmmhitIndexer', 'Hmmer3DomtabHmmqueryIndexer', 'Hmmer3DomtabHmmhitWriter', 'Hmmer3DomtabHmmqueryWriter']

class Hmmer3DomtabParser(Hmmer3TabParser): ...

class Hmmer3DomtabHmmhitParser(Hmmer3DomtabParser):
    hmm_as_hit: bool

class Hmmer3DomtabHmmqueryParser(Hmmer3DomtabParser):
    hmm_as_hit: bool

class Hmmer3DomtabHmmhitIndexer(Hmmer3TabIndexer): ...
class Hmmer3DomtabHmmqueryIndexer(Hmmer3TabIndexer): ...

class Hmmer3DomtabHmmhitWriter:
    hmm_as_hit: bool
    handle: Incomplete
    def __init__(self, handle) -> None: ...
    def write_file(self, qresults): ...

class Hmmer3DomtabHmmqueryWriter(Hmmer3DomtabHmmhitWriter):
    hmm_as_hit: bool
