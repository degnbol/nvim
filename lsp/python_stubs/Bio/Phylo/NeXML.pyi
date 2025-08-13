from Bio.Phylo import BaseTree as BaseTree
from _typeshed import Incomplete

class Tree(BaseTree.Tree):
    weight: Incomplete
    def __init__(self, root: Incomplete | None = None, rooted: bool = False, id: Incomplete | None = None, name: Incomplete | None = None, weight: float = 1.0) -> None: ...

class Clade(BaseTree.Clade):
    comment: Incomplete
    def __init__(self, branch_length: float = 1.0, name: Incomplete | None = None, clades: Incomplete | None = None, confidence: Incomplete | None = None, comment: Incomplete | None = None, **kwargs) -> None: ...
