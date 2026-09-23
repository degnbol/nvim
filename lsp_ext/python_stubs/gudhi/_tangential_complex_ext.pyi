import collections.abc
import gudhi._simplex_tree_ext
from typing import overload

class _Tangential_complex_interface:
    @overload
    def __init__(self, intrinsic_dim: int, points: collections.abc.Sequence[collections.abc.Sequence[float]] = ...) -> None:
        """
        __init__(self, intrinsic_dim: int, points: collections.abc.Sequence[collections.abc.Sequence[float]] = []) -> None
        __init__(self, arg0: int, arg1: ndarray[dtype=float64, shape=(*, *), writable=False], /) -> None
        __init__(self, arg0: int, arg1: str, /) -> None
        """
        ...
    @overload
    def __init__(self, arg0: int, arg1) -> None:
        """
        __init__(self, intrinsic_dim: int, points: collections.abc.Sequence[collections.abc.Sequence[float]] = []) -> None
        __init__(self, arg0: int, arg1: ndarray[dtype=float64, shape=(*, *), writable=False], /) -> None
        __init__(self, arg0: int, arg1: str, /) -> None
        """
        ...
    def compute_tangential_complex(self) -> None:
        """
        compute_tangential_complex(self) -> None

        This function computes the tangential complex.

        Raises:
            ValueError: In debug mode, if the computed star dimension is too low.
                    Try to set a bigger maximal edge length value with
                    :meth:`set_max_squared_edge_length` if this happens.
           
        """
        ...
    def create_simplex_tree(self, arg: gudhi._simplex_tree_ext._Simplex_tree_python_interface) -> None:
        """create_simplex_tree(self, arg: gudhi._simplex_tree_ext._Simplex_tree_python_interface, /) -> None"""
        ...
    def fix_inconsistencies_using_perturbation(self, max_perturb: float, time_limit: float) -> None:
        """
        fix_inconsistencies_using_perturbation(self, max_perturb: float, time_limit: float) -> None

        Attempts to fix inconsistencies by perturbing the point positions.

        :param max_perturb: Maximum length of the translations used by the
            perturbation.
        :type max_perturb: double
        :param time_limit: Time limit in seconds. If -1, no time limit is set.
        :type time_limit: double
           
        """
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
    def num_inconsistent_simplices(self) -> int:
        """
        num_inconsistent_simplices(self) -> int

        :returns:  The number of inconsistent simplices.
        :rtype: unsigned
           
        """
        ...
    def num_inconsistent_stars(self) -> int:
        """
        num_inconsistent_stars(self) -> int

        :returns:  The number of stars containing at least one inconsistent simplex.
        :rtype: unsigned
           
        """
        ...
    def num_simplices(self) -> int:
        """
        num_simplices(self) -> int

        :returns:  Total number of simplices in stars (including duplicates that appear in several stars).
        :rtype: unsigned
           
        """
        ...
    def num_vertices(self) -> int:
        """
        num_vertices(self) -> int

        :returns:  The number of vertices.
        :rtype: unsigned
           
        """
        ...
    def set_max_squared_edge_length(self, max_squared_edge_length: float) -> None:
        """
        set_max_squared_edge_length(self, max_squared_edge_length: float) -> None

        Sets the maximal possible squared edge length for the edges in the
        triangulations.

        :param max_squared_edge_length: Maximal possible squared edge length.
        :type max_squared_edge_length: double

        If the maximal edge length value is too low
        :meth:`compute_tangential_complex`
        will throw an exception in debug mode.
           
        """
        ...
