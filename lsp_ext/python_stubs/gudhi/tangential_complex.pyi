import gudhi._tangential_complex_ext
import gudhi._tangential_complex_ext as t
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class TangentialComplex(gudhi._tangential_complex_ext._Tangential_complex_interface):
    def __init__(self, intrisic_dim, points: Incomplete | None = ..., off_file: str = ...) -> None: ...
    def create_simplex_tree(self): ...
