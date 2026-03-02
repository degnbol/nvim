from _typeshed import Incomplete

class CrossLink:
    featureA: Incomplete
    featureB: Incomplete
    color: Incomplete
    border: Incomplete
    flip: Incomplete
    def __init__(self, featureA, featureB, color=..., border: Incomplete | None = None, flip: bool = False) -> None: ...
    @property
    def startA(self): ...
    @property
    def endA(self): ...
    @property
    def startB(self): ...
    @property
    def endB(self): ...
