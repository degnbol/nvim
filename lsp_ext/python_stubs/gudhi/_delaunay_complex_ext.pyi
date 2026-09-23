import collections.abc
import enum
import gudhi._simplex_tree_ext
import nanobind
from typing import Callable, ClassVar

class Delaunay_complex_interface:
    get_float_relative_precision: ClassVar[nanobind.nb_func] = ...
    """
    get_float_relative_precision() -> float

    Get the float relative precision of filtration values computation when constructing with :code:`precision = 'safe'`
    (the default).

    Returns:
        The float relative precision.
        
    """
    set_float_relative_precision: ClassVar[nanobind.nb_func] = ...
    """
    set_float_relative_precision(precision: float) -> None

    Set the float relative precision of filtration values computation when constructing with :code:`precision = 'safe'`
    (the default).

    Args:
        precision: When constructing :class:`~gudhi.AlphaComplex`, :class:`~gudhi.DelaunayCechComplex`, or
            :class:`~gudhi.DelaunayComplex` with :code:`precision = 'safe'` (the default), one can set the float relative
            precision of filtration values computed. Default is :code:`1e-5` (cf.
            :func:`~gudhi.DelaunayComplex.get_float_relative_precision`). For more details, please refer to
            `CGAL::Lazy_exact_nt<NT>::set_relative_precision_of_to_double <https://doc.cgal.org/latest/Number_types/classCGAL_1_1Lazy__exact__nt.html>`_

    :raises ValueError: If precision is not in (0, 1).
        
    """
    def __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: collections.abc.Sequence[float], arg2: bool, arg3: bool) -> None:
        """
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: collections.abc.Sequence[float], arg2: bool, arg3: bool, /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: ndarray[dtype=float64, shape=(*), writable=False], arg2: bool, arg3: bool, /) -> None
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: bool, arg2: bool, /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: bool, arg2: bool, /) -> None
        """
        ...
    def create_simplex_tree(self, arg0: gudhi._simplex_tree_ext._Simplex_tree_python_interface, arg1: float, arg2: Filtration, arg3: bool) -> None:
        """create_simplex_tree(self, arg0: gudhi._simplex_tree_ext._Simplex_tree_python_interface, arg1: float, arg2: gudhi._delaunay_complex_ext.Filtration, arg3: bool, /) -> None"""
        ...
    def get_point(self, vertex: int) -> list[float]:
        """
        get_point(self, vertex: int) -> list[float]

        This function returns the point corresponding to a given vertex from the :class:`~gudhi.SimplexTree` (the
        same as the k-th input point, where `k=vertex`)

        Args:
            vertex: The vertex.
        Returns:
            the point.

        :raises IndexError: In case the point has no associated vertex in the diagram (because of weights or because it
            is a duplicate).
        
        """
        ...

class Filtration(enum.Enum):
    __new__: ClassVar[Callable] = ...
    ALPHA: ClassVar[Filtration] = ...
    """Alpha Complex"""
    CECH: ClassVar[Filtration] = ...
    """Delaunay Cech Complex"""
    NONE: ClassVar[Filtration] = ...
    """Default Delaunay Complex"""
    _generate_next_value_: ClassVar[Callable] = ...
    _member_map_: ClassVar[dict] = ...
    _member_names_: ClassVar[list] = ...
    _member_type_: ClassVar[type[object]] = ...
    _unhashable_values_: ClassVar[list] = ...
    _use_args_: ClassVar[bool] = ...
    _value2member_map_: ClassVar[dict] = ...
    _value_repr_: ClassVar[None] = ...
    __nb_enum__: ClassVar[PyCapsule] = ...
