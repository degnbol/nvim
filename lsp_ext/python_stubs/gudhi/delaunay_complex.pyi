import collections.abc
import gudhi._delaunay_complex_ext
import gudhi._delaunay_complex_ext as t
import gudhi.simplex_tree
import typing
from _typeshed import Incomplete
from gudhi.simplex_tree import SimplexTree as SimplexTree

class DelaunayComplex(gudhi._delaunay_complex_ext.Delaunay_complex_interface):
    def __init__(self, points: collections.abc.Sequence = ..., weights: Incomplete = ..., precision: typing.Literal = ...) -> None: ...
    def create_simplex_tree(self, max_alpha_square: float = ..., filtration: Incomplete = ..., output_squared_values: bool = ...) -> gudhi.simplex_tree.SimplexTree: ...

class AlphaComplex(DelaunayComplex):
    def create_simplex_tree(self, max_alpha_square: float = ..., default_filtration_value: bool = ..., output_squared_values: bool = ...) -> gudhi.simplex_tree.SimplexTree: ...

class DelaunayCechComplex(DelaunayComplex):
    def __init__(self, points: list = ..., precision: str = ...) -> None: ...
    def create_simplex_tree(self, max_alpha_square: float = ..., output_squared_values: bool = ...) -> gudhi.simplex_tree.SimplexTree: ...
