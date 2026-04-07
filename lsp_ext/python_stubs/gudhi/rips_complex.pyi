import gudhi._rips_complex_ext
import gudhi._rips_complex_ext as t
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class RipsComplex(gudhi._rips_complex_ext.Rips_complex_interface):
    def __init__(self, *, points=..., distance_matrix=..., max_edge_length: float = ..., sparse: Incomplete = ...) -> None: ...
    def create_simplex_tree(self, max_dimension: int = ...): ...
