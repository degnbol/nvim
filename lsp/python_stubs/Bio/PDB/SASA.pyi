from _typeshed import Incomplete

__all__ = ['ShrakeRupley']

class ShrakeRupley:
    probe_radius: Incomplete
    n_points: Incomplete
    radii_dict: Incomplete
    def __init__(self, probe_radius: float = 1.4, n_points: int = 100, radii_dict: Incomplete | None = None) -> None: ...
    def compute(self, entity, level: str = 'A') -> None: ...
