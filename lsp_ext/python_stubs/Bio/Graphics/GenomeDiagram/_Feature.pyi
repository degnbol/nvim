from ._Colors import ColorTranslator as ColorTranslator
from _typeshed import Incomplete

class Feature:
    parent: Incomplete
    id: Incomplete
    color: Incomplete
    border: Incomplete
    hide: int
    sigil: str
    arrowhead_length: float
    arrowshaft_height: float
    name_qualifiers: Incomplete
    label: Incomplete
    label_font: str
    label_size: int
    label_color: Incomplete
    label_angle: int
    label_position: Incomplete
    label_strand: Incomplete
    def __init__(self, parent: Incomplete | None = None, feature_id: Incomplete | None = None, feature: Incomplete | None = None, color=..., label: int = 0, border: Incomplete | None = None, colour: Incomplete | None = None) -> None: ...
    def set_feature(self, feature) -> None: ...
    def get_feature(self): ...
    def set_colour(self, colour) -> None: ...
    def set_color(self, color) -> None: ...
    def __getattr__(self, name): ...
