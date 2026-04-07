import gudhi._euclidean_witness_complex_ext
import gudhi._euclidean_witness_complex_ext as t
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class EuclideanWitnessComplex(gudhi._euclidean_witness_complex_ext.Euclidean_witness_complex_interface):
    def __init__(self, landmarks: Incomplete | None = ..., witnesses: Incomplete | None = ...) -> None: ...
    def create_simplex_tree(self, max_alpha_square, limit_dimension: int = ...): ...
