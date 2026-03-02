from _typeshed import Incomplete
from collections.abc import Generator

class ColorSpiral:
    a: Incomplete
    b: Incomplete
    v_init: Incomplete
    v_final: Incomplete
    jitter: Incomplete
    def __init__(self, a: int = 1, b: float = 0.33, v_init: float = 0.85, v_final: float = 0.5, jitter: float = 0.05) -> None: ...
    def get_colors(self, k, offset: float = 0.1) -> Generator[Incomplete]: ...

def get_colors(k, **kwargs): ...
def get_color_dict(l, **kwargs): ...
