import collections.abc

class _Nerve_gic_interface:
    def __init__(self) -> None:
        """__init__(self) -> None"""
        ...
    def compute_PD(self) -> list[tuple[float, float]]:
        """
        compute_PD(self) -> list[tuple[float, float]]

        Computes the extended persistence diagram of the complex.
        
        """
        ...
    def compute_confidence_level_from_distance(self, distance: float) -> float:
        """
        compute_confidence_level_from_distance(self, distance: float) -> float

        Computes the confidence level of a specific bottleneck distance threshold.

        :param distance: Bottleneck distance.
        :type distance: double

        :rtype: double
        :returns: Confidence level.
        
        """
        ...
    def compute_distance_from_confidence_level(self, alpha: float) -> float:
        """
        compute_distance_from_confidence_level(self, alpha: float) -> float

        Computes the bottleneck distance threshold corresponding to a
        specific confidence level.

        :param alpha: Confidence level.
        :type alpha: double

        :rtype: double
        :returns: Bottleneck distance.
        
        """
        ...
    def compute_distribution(self, N: int = ...) -> None:
        """
        compute_distribution(self, N: int = 100) -> None

        Computes bootstrapped distances distribution.

        :param N: Loop number (default value is 100).
        :type alpha: int
        
        """
        ...
    def compute_p_value(self) -> float:
        """
        compute_p_value(self) -> float

        Computes the p-value, i.e. the opposite of the confidence level of
        the largest bottleneck distance preserving the points
        persistence diagram of the output simplicial complex.

        :rtype: double
        :returns: p-value.
        
        """
        ...
    def find_simplices(self) -> None:
        """
        find_simplices(self) -> None

        Computes the simplices of the simplicial complex.
        
        """
        ...
    def plot_dot(self) -> None:
        """
        plot_dot(self) -> None

        Creates a .dot file called SC.dot for neato (part of the graphviz
        package) once the simplicial complex is computed to get a visualization of
        its 1-skeleton in a .pdf file.
        
        """
        ...
    def plot_off(self) -> None:
        """
        plot_off(self) -> None

        Creates a .off file called SC.off for 3D visualization, which contains
        the 2-skeleton of the GIC. This function assumes that the cover has been
        computed with Voronoi. If data points are in 1D or 2D, the remaining
        coordinates of the points embedded in 3D are set to 0.
        
        """
        ...
    def set_automatic_resolution(self) -> float:
        """
        set_automatic_resolution(self) -> float

        Computes the optimal length of intervals (i.e. the smallest interval
        length avoiding discretization artifacts - see :cite:`Carriere17c`) for a
        functional cover.

        :rtype: double
        :returns: resolution interval length used to compute the cover.
        
        """
        ...
    def set_color_from_coordinate(self, k: int = ...) -> None:
        """
        set_color_from_coordinate(self, k: int = 0) -> None

        Computes the function used to color the nodes of the simplicial
        complex from the k-th coordinate.

        :param k: Coordinate to use (start at 0). Default value is 0.
        :type k: int
        
        """
        ...
    def set_color_from_range(self, color: collections.abc.Sequence[float]) -> None:
        """
        set_color_from_range(self, color: collections.abc.Sequence[float]) -> None

        Computes the function used to color the nodes of the simplicial
        complex from a vector stored in memory.

        :param color: Input vector of values.
        :type color: vector[double]
        
        """
        ...
    def set_cover_from_Voronoi(self, m: int = ...) -> None:
        """
        set_cover_from_Voronoi(self, m: int = 100) -> None

        Creates the cover C from the Voronoï cells of a subsampling of the point cloud.

        :param m: Number of points in the subsample. Default value is 100.
        :type m: int
        
        """
        ...
    def set_cover_from_function(self) -> None:
        """
        set_cover_from_function(self) -> None

        Creates a cover C from the preimages of the function f.
        
        """
        ...
    def set_cover_from_range(self, assignments: collections.abc.Sequence[collections.abc.Sequence[int]]) -> None:
        """
        set_cover_from_range(self, assignments: collections.abc.Sequence[collections.abc.Sequence[int]]) -> None

        Creates a cover C from a vector stored in memory.

        :param assignments: Vector containing the assignments of the points to their corresponding cover elements. For instance, if the i-th point belongs to the 1st and 3rd cover elements, then assignments[i] = [1,3].
        :type assignments: List[List[int]]
        
        """
        ...
    def set_distances_from_range(self, distance_matrix: collections.abc.Sequence[collections.abc.Sequence[float]]) -> None:
        """
        set_distances_from_range(self, distance_matrix: collections.abc.Sequence[collections.abc.Sequence[float]]) -> None

        Reads and stores the input distance matrix from a vector stored in memory.

        :param distance_matrix: Input vector containing the distance matrix.
        :type distance_matrix: vector[vector[double]]
        
        """
        ...
    def set_function_from_coordinate(self, k: int) -> None:
        """
        set_function_from_coordinate(self, k: int) -> None

        Creates the function f from the k-th coordinate of the point cloud.

        :param k: Coordinate to use (start at 0).
        :type k: int
        
        """
        ...
    def set_function_from_range(self, function: collections.abc.Sequence[float]) -> None:
        """
        set_function_from_range(self, function: collections.abc.Sequence[float]) -> None

        Creates the function f from a vector stored in memory.

        :param function: Input vector of values.
        :type function: vector[double]
        
        """
        ...
    def set_gain(self, g: float = ...) -> None:
        """
        set_gain(self, g: float = 0.3) -> None

        Sets a gain from a value stored in memory.

        :param g: Gain (default value is 0.3).
        :type g: double
        
        """
        ...
    def set_graph_from_OFF(self) -> None:
        """
        set_graph_from_OFF(self) -> None

        Creates a graph G from the triangulation given by the input OFF file.
        
        """
        ...
    def set_graph_from_automatic_rips(self, N: int = ...) -> float:
        """
        set_graph_from_automatic_rips(self, N: int = 100) -> float

        Creates a graph G from a Rips complex whose threshold value is
        automatically tuned with subsampling - see :cite:`Carriere17c`.

        :param N: Number of subsampling iteration (the default reasonable value is 100, but there is no guarantee on how to choose it).
        :type N: int
        :rtype: double
        :returns: Delta threshold used for computing the Rips complex.
        
        """
        ...
    def set_graph_from_rips(self, threshold: float) -> None:
        """
        set_graph_from_rips(self, threshold: float) -> None

        Creates a graph G from a Rips complex.

        :param threshold: Threshold value for the Rips complex.
        :type threshold: double
        
        """
        ...
    def set_mask(self, nodemask: int) -> None:
        """
        set_mask(self, nodemask: int) -> None

        Sets the mask, which is a threshold integer such that nodes in the
        complex that contain a number of data points which is less than or
        equal to this threshold are not displayed.

        :param nodemask: Threshold.
        :type nodemask: int
        
        """
        ...
    def set_point_cloud_from_range(self, cloud: collections.abc.Sequence[collections.abc.Sequence[float]]) -> None:
        """
        set_point_cloud_from_range(self, cloud: collections.abc.Sequence[collections.abc.Sequence[float]]) -> None

        Reads and stores the input point cloud from a vector stored in memory.

        :param cloud: Input vector containing the point cloud.
        :type cloud: vector[vector[double]]
        
        """
        ...
    def set_resolution_with_interval_length(self, resolution: float) -> None:
        """
        set_resolution_with_interval_length(self, resolution: float) -> None

        Sets a length of intervals from a value stored in memory.

        :param resolution: Length of intervals.
        :type resolution: double
        
        """
        ...
    def set_resolution_with_interval_number(self, resolution: int) -> None:
        """
        set_resolution_with_interval_number(self, resolution: int) -> None

        "Sets a number of intervals from a value stored in memory.

        :param resolution: Number of intervals.
        :type resolution: int
        
        """
        ...
    def set_subsampling(self, constant: float, power: float) -> None:
        """
        set_subsampling(self, constant: float, power: float) -> None

        Sets the constants used to subsample the data set. These constants
        are explained in :cite:`Carriere17c`.

        :param constant: Constant.
        :type constant: double

        :param power: Power.
        :type resolution: double
        
        """
        ...
    def set_type(self, type: str) -> None:
        """
        set_type(self, type: str) -> None

        Specifies whether the type of the output simplicial complex.

        :param type: either "GIC" or "Nerve".
        :type type: string
        
        """
        ...
    def set_verbose(self, verbose: bool = ...) -> None:
        """
        set_verbose(self, verbose: bool = False) -> None

        Specifies whether the program should display information or not.

        :param verbose: true = display info, false = do not display info.
        :type verbose: boolean
        
        """
        ...
    def subcolor(self, c: int) -> float:
        """
        subcolor(self, c: int) -> float

        Returns the mean color value corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int

        :rtype: float
        :returns: Mean color value of data points.
        
        """
        ...
    def subpopulation(self, c: int) -> list[int]:
        """
        subpopulation(self, c: int) -> list[int]

        Returns the data subset corresponding to a specific node of the
        created complex.

        :param c: ID of the node.
        :type c: int

        :rtype: vector[int]
        :returns: Vector of IDs of data points.
        
        """
        ...
    def write_info(self) -> None:
        """
        write_info(self) -> None

        Creates a .txt file called SC.txt describing the 1-skeleton, which can
        then be plotted with e.g. KeplerMapper.
        
        """
        ...
