import collections.abc
import gudhi._simplex_tree_ext

class Rips_complex_interface:
    def __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: float, arg2: bool) -> None:
        """
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: float, arg2: bool, /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: float, arg2: bool, /) -> None
        __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[float]], arg1: float, arg2: float, arg3: bool, /) -> None
        __init__(self, arg0: ndarray[dtype=float64, shape=(*, *), writable=False], arg1: float, arg2: float, arg3: bool, /) -> None
        """
        ...
    def create_simplex_tree(self, arg0: gudhi._simplex_tree_ext._Simplex_tree_python_interface, arg1: int) -> None:
        """create_simplex_tree(self, arg0: gudhi._simplex_tree_ext._Simplex_tree_python_interface, arg1: int, /) -> None"""
        ...
