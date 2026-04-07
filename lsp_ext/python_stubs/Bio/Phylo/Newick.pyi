from Bio.Phylo import BaseTree as BaseTree
from _typeshed import Incomplete

class Tree(BaseTree.Tree):
    weight: Incomplete
    def __init__(self, root=None, rooted: bool = False, id=None, name=None, weight: float = 1.0) -> None: ...

class Clade(BaseTree.Clade):
    comment: Incomplete
    def __init__(self, branch_length=None, name=None, clades=None, confidence=None, comment=None) -> None: ...
