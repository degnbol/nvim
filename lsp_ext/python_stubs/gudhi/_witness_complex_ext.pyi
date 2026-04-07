import collections.abc
import gudhi._simplex_tree_ext
from typing import overload

class Witness_complex_interface:
    @overload
    def __init__(self, arg: collections.abc.Sequence[collections.abc.Sequence[tuple[int, float]]]) -> None: ...
    @overload
    def __init__(self) -> None: ...
    def create_simplex_tree(self, simplex_tree: gudhi._simplex_tree_ext._Simplex_tree_python_interface, max_alpha_square: float, limit_dimension: int = ...) -> None: ...
