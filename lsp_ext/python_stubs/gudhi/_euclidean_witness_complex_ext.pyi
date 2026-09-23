import collections.abc
import gudhi._simplex_tree_ext
from typing import overload

class Euclidean_witness_complex_interface:
    @overload
    def __init__(self) -> None:
        """
        __init__(self) -> None
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: collections.abc.Sequence[collections.abc.Sequence[float]], /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: ndarray[dtype=float64, shape=(*, *), writable=False], /) -> None
        """
        ...
    @overload
    def __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: collections.abc.Sequence[collections.abc.Sequence[float]]) -> None:
        """
        __init__(self) -> None
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: collections.abc.Sequence[collections.abc.Sequence[float]], /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: ndarray[dtype=float64, shape=(*, *), writable=False], /) -> None
        """
        ...
    def create_simplex_tree(self, simplex_tree: gudhi._simplex_tree_ext._Simplex_tree_python_interface, max_alpha_square: float, limit_dimension: int = ...) -> None:
        """create_simplex_tree(self, simplex_tree: gudhi._simplex_tree_ext._Simplex_tree_python_interface, max_alpha_square: float, limit_dimension: int = 18446744073709551615) -> None"""
        ...
    def get_point(self, vertex: int) -> list[float]:
        """
        get_point(self, vertex: int) -> list[float]

        This function returns the point corresponding to a given vertex.

        :param vertex: The vertex.
        :type vertex: int.
        :returns:  The point.
        :rtype: list of float
           
        """
        ...
