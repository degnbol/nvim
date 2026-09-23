import collections.abc

class _Bitmap_cubical_complex_interface:
    def __init__(self, arg0: collections.abc.Sequence[int], arg1: collections.abc.Sequence[float], arg2: bool) -> None:
        """
        __init__(self, arg0: collections.abc.Sequence[int], arg1: collections.abc.Sequence[float], arg2: bool, /) -> None
        __init__(self, arg: str, /) -> None
        """
        ...
    def dimension(self) -> int:
        """
        dimension(self) -> int

        This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
           
        """
        ...
    def num_simplices(self) -> int:
        """
        num_simplices(self) -> int

        This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
           
        """
        ...

class _Cubical_complex_persistence_interface:
    def __init__(self, arg0: _Bitmap_cubical_complex_interface, arg1: bool) -> None:
        """__init__(self, arg0: gudhi._cubical_complex_ext._Bitmap_cubical_complex_interface, arg1: bool, /) -> None"""
        ...

class _Periodic_cubical_complex_interface:
    def __init__(self, arg0: collections.abc.Sequence[int], arg1: collections.abc.Sequence[float], arg2: collections.abc.Sequence[bool], arg3: bool) -> None:
        """
        __init__(self, arg0: collections.abc.Sequence[int], arg1: collections.abc.Sequence[float], arg2: collections.abc.Sequence[bool], arg3: bool, /) -> None
        __init__(self, arg: str, /) -> None
        """
        ...
    def dimension(self) -> int:
        """
        dimension(self) -> int

        This function returns the dimension of the complex.

        :returns:  int -- the complex dimension.
           
        """
        ...
    def num_simplices(self) -> int:
        """
        num_simplices(self) -> int

        This function returns the number of all cubes in the complex.

        :returns:  int -- the number of all cubes in the complex.
           
        """
        ...

class _Periodic_cubical_complex_persistence_interface:
    def __init__(self, arg0: _Periodic_cubical_complex_interface, arg1: bool) -> None:
        """__init__(self, arg0: gudhi._cubical_complex_ext._Periodic_cubical_complex_interface, arg1: bool, /) -> None"""
        ...
