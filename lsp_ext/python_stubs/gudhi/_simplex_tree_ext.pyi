import collections.abc
from typing import overload

class _Simplex_tree_persistence_interface:
    def __init__(self, arg0: _Simplex_tree_python_interface, arg1: bool) -> None:
        """__init__(self, arg0: gudhi._simplex_tree_ext._Simplex_tree_python_interface, arg1: bool, /) -> None"""
        ...

class _Simplex_tree_python_interface:
    @overload
    def __init__(self) -> None:
        """
        __init__(self) -> None
        __init__(self, arg: gudhi._simplex_tree_ext._Simplex_tree_python_interface) -> None
        """
        ...
    @overload
    def __init__(self, arg: _Simplex_tree_python_interface) -> None:
        """
        __init__(self) -> None
        __init__(self, arg: gudhi._simplex_tree_ext._Simplex_tree_python_interface) -> None
        """
        ...
    def assign_filtration(self, simplex: collections.abc.Sequence[int], filtration: float) -> None:
        """
        assign_filtration(self, simplex: collections.abc.Sequence[int], filtration: float) -> None

        This function assigns a new filtration value to a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :param filtration:  The new filtration value.
        :type filtration:  float

        .. note::

            Beware that after this operation, the structure may not be a valid filtration anymore, a simplex could have a lower
            filtration value than one of its faces. Callers are responsible for fixing this (with more
            :meth:`assign_filtration` or :meth:`make_filtration_non_decreasing` for instance) before calling any function that
            relies on the filtration property, like :meth:`persistence`.
           
        """
        ...
    def clear(self) -> None:
        """
        clear(self) -> None

        Remove all the simplices, leaving an empty complex.
           
        """
        ...
    def dimension(self) -> int:
        """
        dimension(self) -> int

        This function returns the dimension of the simplicial complex.

        :returns:  the simplicial complex dimension.
        :rtype:  int

        .. note::

            This function is not constant time because it can recompute dimension if required (can be triggered by
            :func:`remove_maximal_simplex` or  :func:`prune_above_filtration` methods).
           
        """
        ...
    def euler_characteristic(self) -> int:
        """
        euler_characteristic(self) -> int

        This function computes and returns the Euler characteristic of the non-filtered underlying complex represented
        by the simplex tree.

        :returns:  The Euler characteristic.
        :rtype:  int
           
        """
        ...
    def expansion(self, max_dimension: int) -> None:
        """
        expansion(self, max_dimension: int) -> None

        Expands the simplex tree containing only its one skeleton until dimension max_dim.

        The expanded simplicial complex until dimension :math:`d` attached to a graph :math:`G` is the maximal simplicial
        complex of dimension at most :math:`d` admitting the graph :math:`G` as :math:`1`-skeleton.
        The filtration value assigned to a simplex is the maximal filtration value of one of its edges.

        The simplex tree must contain no simplex of dimension bigger than 1 when calling the method.

        :param max_dimension: The maximal dimension.
        :type max_dimension: int
           
        """
        ...
    def expansion_with_blocker(self, max_dim: int, blocker_func: collections.abc.Callable[[collections.abc.Sequence[int]], bool]) -> None:
        """
        expansion_with_blocker(self, max_dim: int, blocker_func: collections.abc.Callable[[collections.abc.Sequence[int]], bool]) -> None

        Expands the SimplexTree containing only a graph. Simplices corresponding to cliques in the graph are added
        incrementally, faces before cofaces, unless the simplex has dimension larger than `max_dim` or `blocker_func`
        returns `True` for this simplex.

        The function identifies a candidate simplex whose faces are all already in the complex, inserts it with a
        filtration value corresponding to the maximum of the filtration values of the faces, then calls `blocker_func`
        with this new simplex (represented as a list of int). If `blocker_func` returns `True`, the simplex is removed,
        otherwise it is kept. The algorithm then proceeds with the next candidate.

        .. warning::

            Several candidates of the same dimension may be inserted simultaneously before calling `blocker_func`, so
            if you examine the complex in `blocker_func`, you may hit a few simplices of the same dimension that have
            not been vetted by `blocker_func` yet, or have already been rejected but not yet removed.

        :param max_dim: Expansion maximal dimension value.
        :type max_dim: int
        :param blocker_func: Blocker oracle.
        :type blocker_func: Callable[[List[int]], bool]
           
        """
        ...
    def extend_filtration(self) -> None:
        """
        extend_filtration(self) -> None

        Extend filtration for computing extended persistence. This function only uses the filtration values at the
        0-dimensional simplices, and computes the extended persistence diagram induced by the lower-star filtration
        computed with these values.

        .. note::

            Note that after calling this function, the filtration values are actually modified within the simplex tree.
            The function :func:`extended_persistence` retrieves the original values.

        .. note::

            Note that this code creates an extra vertex internally, so you should make sure that the simplex tree does not
            contain a vertex with the largest possible value (i.e., 4294967295).

        This `notebook <https://github.com/GUDHI/TDA-tutorial/blob/master/tutorials/Tuto-GUDHI-extended-persistence.ipynb>`_
        explains how to compute an extension of persistence called extended persistence.
           
        """
        ...
    def filtration(self, simplex: collections.abc.Sequence[int]) -> float:
        """
        filtration(self, simplex: collections.abc.Sequence[int]) -> float

        This function returns the filtration value for a given N-simplex in this simplicial complex, or +infinity if it is not
        in the complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :returns:  The simplicial complex filtration value.
        :rtype:  float
           
        """
        ...
    def find(self, simplex: collections.abc.Sequence[int]) -> bool:
        """
        find(self, simplex: collections.abc.Sequence[int]) -> bool

        This function returns if the N-simplex was found in the simplicial complex or not.

        :param simplex: The N-simplex to find, represented by a list of vertex.
        :type simplex: list of int
        :returns:  `True` if the simplex was found, `False` otherwise.
        :rtype:  bool
           
        """
        ...
    def get_boundaries(self, simplex: collections.abc.Sequence[int]) -> collections.abc.Iterator[tuple[list[int], float]]:
        """
        get_boundaries(self, simplex: collections.abc.Sequence[int]) -> collections.abc.Iterator[tuple[list[int], float]]

        This function returns an iterator over the boundaries of a given N-simplex.
        If you do not need the filtration values, the boundary can also be obtained as
        :code:`itertools.combinations(simplex,len(simplex)-1)`.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int.
        :returns:  The (simplices of the) boundary of a simplex
        :rtype:  Iterator over tuples(simplex, filtration)
          
        """
        ...
    def get_cofaces(self, simplex: collections.abc.Sequence[int], dimension: int) -> list[tuple[list[int], float]]:
        """
        get_cofaces(self, simplex: collections.abc.Sequence[int], dimension: int) -> list[tuple[list[int], float]]

        This function returns the cofaces of a given N-simplex with a given codimension.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :param codimension: The codimension. If codimension = 0, all cofaces are returned (equivalent of get_star function)
        :type codimension: int
        :returns:  The (simplices of the) cofaces of a simplex
        :rtype:  list of tuples(simplex, filtration)
           
        """
        ...
    def get_filtration(self) -> collections.abc.Iterator[tuple[list[int], float]]:
        """
        get_filtration(self) -> collections.abc.Iterator[tuple[list[int], float]]

        This function returns an iterator over simplices and their given
        filtration values sorted by increasing filtration values.

        :returns:  The simplices sorted by increasing filtration values.
        :rtype:  Iterator over tuples(simplex, filtration)
          
        """
        ...
    def get_simplices(self) -> collections.abc.Iterator[tuple[list[int], float]]:
        """
        get_simplices(self) -> collections.abc.Iterator[tuple[list[int], float]]

        This function returns an iterator over simplices and their given filtration values.

        :returns:  The simplices.
        :rtype:  Iterator over tuples(simplex, filtration)
          
        """
        ...
    def get_skeleton(self, dimension: int) -> collections.abc.Iterator[tuple[list[int], float]]:
        """
        get_skeleton(self, dimension: int) -> collections.abc.Iterator[tuple[list[int], float]]

        This function returns an iterator over the (simplices of the) skeleton of a maximum given dimension.

        :param dimension: The skeleton dimension value.
        :type dimension: int
        :returns:  The (simplices of the) skeleton of a maximum dimension.
        :rtype:  Iterator over tuples(simplex, filtration)
          
        """
        ...
    def get_star(self, simplex: collections.abc.Sequence[int]) -> list[tuple[list[int], float]]:
        """
        get_star(self, simplex: collections.abc.Sequence[int]) -> list[tuple[list[int], float]]

        This function returns the star of a given N-simplex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int
        :returns:  The (simplices of the) star of a simplex.
        :rtype:  list of tuples(simplex, filtration)
           
        """
        ...
    def insert(self, simplex: collections.abc.Sequence[int], filtration: float = ...) -> bool:
        """
        insert(self, simplex: collections.abc.Sequence[int], filtration: float = 0.0) -> bool

        This function inserts the given N-simplex and its subfaces with the given filtration value (default value is '0.0'). If
        some of those simplices are already present with a higher filtration value, their filtration value is lowered.

        :param simplex: The N-simplex to insert, represented by a list of vertex.
        :type simplex: list of int
        :param filtration: The filtration value of the simplex.
        :type filtration: float
        :returns:  `True` if the simplex was not yet in the complex, `False` otherwise (whatever its original filtration value)
        :rtype:  bool
           
        """
        ...
    def is_empty(self) -> bool:
        """
        is_empty(self) -> bool

        This function returns whether the simplicial complex is empty.

        :returns:  `True` if the simplicial complex is empty.
        :rtype:  bool
           
        """
        ...
    def make_filtration_non_decreasing(self) -> bool:
        """
        make_filtration_non_decreasing(self) -> bool

        This function ensures that each simplex has a higher filtration value than its faces by increasing the filtration
        values.

        :returns: `True` if any filtration value was modified, `False` if the filtration was already non-decreasing.
        :rtype: bool
           
        """
        ...
    def num_simplices(self) -> int:
        """
        num_simplices(self) -> int

        This function returns the number of simplices of the simplicial complex.

        :returns:  the simplicial complex number of simplices.
        :rtype:  int
           
        """
        ...
    def num_simplices_by_dimension(self, *args, **kwargs):
        """
        num_simplices_by_dimension(self) -> numpy.ndarray[dtype=uint64]

        Computes and returns the number of simplices of each dimension in the complex.

        :returns:  Array containing at index `d` the number of simplices of dimension `d` in the complex.
        :rtype:  1-dimensional numpy.array of length `maximal dimension + 1`.
           
        """
        ...
    def num_vertices(self) -> int:
        """
        num_vertices(self) -> int

        This function returns the number of vertices of the simplicial complex.

        :returns:  The simplicial complex number of vertices.
        :rtype:  int
           
        """
        ...
    def prune_above_dimension(self, dimension: int) -> bool:
        """
        prune_above_dimension(self, dimension: int) -> bool

        Remove all simplices of dimension greater than a given value.

        :param dimension: Maximum dimension value.
        :type dimension: int
        :returns: The modification information.
        :rtype: bool
           
        """
        ...
    def prune_above_filtration(self, filtration: float) -> bool:
        """
        prune_above_filtration(self, filtration: float) -> bool

        Prune above filtration value given as parameter.

        :param filtration: Maximum threshold value.
        :type filtration: float
        :returns: The filtration modification information.
        :rtype: bool

        .. note::

            Note that the dimension of the simplicial complex may be lower after calling :func:`prune_above_filtration` than it
            was before. However, :func:`upper_bound_dimension` will return the old value, which remains a  valid upper bound.
            If you care, you can call :func:`dimension` method to recompute the exact dimension.
           
        """
        ...
    def remove_maximal_simplex(self, simplex: collections.abc.Sequence[int]) -> None:
        """
        remove_maximal_simplex(self, simplex: collections.abc.Sequence[int]) -> None

        This function removes a given maximal N-simplex from the simplicial complex.

        :param simplex: The N-simplex, represented by a list of vertex.
        :type simplex: list of int

        .. note::

            The dimension of the simplicial complex may be lower after calling remove_maximal_simplex than it was before.
            However, :func:`upper_bound_dimension` method will return the old value, which remains a valid upper bound. If you
            care, you can call :func:`dimension` to recompute the exact dimension.
           
        """
        ...
    def reset_filtration(self, filtration: float, min_dim: int = ...) -> None:
        """
        reset_filtration(self, filtration: float, min_dim: int = 0) -> None

        This function resets the filtration value of all the simplices of dimension at least min_dim. Resets all the simplex
        tree when `min_dim = 0`.
        `reset_filtration` may break the filtration property with `min_dim > 0`, and it is the user's responsibility to make it
        a valid filtration (using a large enough `filt_value`, or calling `make_filtration_non_decreasing` afterwards for
        instance).

        :param filtration: New threshold value.
        :type filtration: float.
        :param min_dim: The minimal dimension. Default value is 0.
        :type min_dim: int.
           
        """
        ...
    def set_dimension(self, dimension: int, exact: bool = ...) -> None:
        """
        set_dimension(self, dimension: int, exact: bool = True) -> None

        This function sets the dimension of the simplicial complex.

        :param dimension: The new dimension value.
        :type dimension: int

        .. note::

            This function must be used with caution because it disables dimension recomputation when required (this
            recomputation can be triggered by  :func:`remove_maximal_simplex` or :func:`prune_above_filtration` ).
           
        """
        ...
    def upper_bound_dimension(self) -> int:
        """
        upper_bound_dimension(self) -> int

        This function returns a valid dimension upper bound of the simplicial complex.

        :returns:  an upper bound on the dimension of the simplicial complex.
        :rtype:  int
           
        """
        ...
    def __eq__(self, other: _Simplex_tree_python_interface) -> bool:
        """
        __eq__(self, other: gudhi._simplex_tree_ext._Simplex_tree_python_interface) -> bool

        Equality operator in order to compare 2 SimplexTree data structures.

        :returns: `True` if the 2 complexes have the same simplices with the same filtration values, `False` otherwise.
        :rtype: bool
           
        """
        ...
