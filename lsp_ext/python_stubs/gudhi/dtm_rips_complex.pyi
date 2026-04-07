import gudhi.weighted_rips_complex
from _typeshed import Incomplete
from gudhi.point_cloud.dtm import DistanceToMeasure as DistanceToMeasure
from gudhi.weighted_rips_complex import WeightedRipsComplex as WeightedRipsComplex

class DTMRipsComplex(gudhi.weighted_rips_complex.WeightedRipsComplex):
    def __init__(self, points: Incomplete | None = ..., distance_matrix: Incomplete | None = ..., k: int = ..., q: int = ..., max_filtration: float = ...) -> None: ...
