import pyrosetta.distributed
from _typeshed import Incomplete

__all__ = ['movers', 'filters']

class ComponentDoc:
    name: Incomplete
    doc: Incomplete
    __doc__: Incomplete
    def __init__(self, name, doc) -> None: ...

class InlineDocs:
    @pyrosetta.distributed.requires_init
    def __dir__(self): ...
    def __getattr__(self, name): ...
    @pyrosetta.distributed.requires_init
    def get_component_doc(self, name): ...

class MoverDocs(InlineDocs): ...
class FilterDocs(InlineDocs): ...

movers: Incomplete
filters: Incomplete
