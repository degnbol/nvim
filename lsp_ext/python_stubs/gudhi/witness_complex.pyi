import gudhi._witness_complex_ext
import gudhi._witness_complex_ext as t
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class WitnessComplex(gudhi._witness_complex_ext.Witness_complex_interface):
    def __init__(self, nearest_landmark_table: Incomplete | None = ...) -> None: ...
    def create_simplex_tree(self, max_alpha_square: float = ..., limit_dimension: int = ...): ...
