from Bio.Phylo import BaseTree as BaseTree
from _typeshed import Incomplete

class Tree(BaseTree.Tree):
    weight: Incomplete
    attributes: Incomplete
    def __init__(self, root=None, rooted: bool = False, id=None, name=None, weight: float = 1.0) -> None: ...

class Clade(BaseTree.Clade):
    comment: Incomplete
    attributes: Incomplete
    tu_attributes: Incomplete
    edge_attributes: Incomplete
    def __init__(self, branch_length: float = 1.0, name=None, clades=None, confidence=None, comment=None) -> None: ...
