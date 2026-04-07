import gudhi._euclidean_strong_witness_complex_ext
import gudhi._euclidean_strong_witness_complex_ext as t
import gudhi.simplex_tree
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class EuclideanStrongWitnessComplex(gudhi._euclidean_strong_witness_complex_ext.Euclidean_strong_witness_complex_interface):
    def __init__(self, landmarks: Incomplete | None = ..., witnesses: Incomplete | None = ...) -> None: ...
    def create_simplex_tree(self, max_alpha_square: float = ..., limit_dimension: int = ...) -> gudhi.simplex_tree.SimplexTree: ...
