from ._FeatureSet import FeatureSet as FeatureSet
from ._GraphSet import GraphSet as GraphSet
from _typeshed import Incomplete

class Track:
    height: Incomplete
    name: Incomplete
    hide: Incomplete
    start: Incomplete
    end: Incomplete
    greytrack: Incomplete
    greytrack_labels: Incomplete
    greytrack_fontsize: Incomplete
    greytrack_font: Incomplete
    greytrack_font_rotation: Incomplete
    greytrack_fontcolor: Incomplete
    scale: Incomplete
    scale_format: Incomplete
    scale_color: Incomplete
    scale_font: Incomplete
    scale_fontsize: Incomplete
    scale_fontangle: Incomplete
    scale_ticks: Incomplete
    scale_largeticks: Incomplete
    scale_smallticks: Incomplete
    scale_largetick_interval: Incomplete
    scale_smalltick_interval: Incomplete
    scale_largetick_labels: Incomplete
    scale_smalltick_labels: Incomplete
    axis_labels: Incomplete
    def __init__(self, name: Incomplete | None = None, height: int = 1, hide: int = 0, greytrack: int = 0, greytrack_labels: int = 5, greytrack_fontsize: int = 8, greytrack_font: str = 'Helvetica', greytrack_font_rotation: int = 0, greytrack_font_color=..., scale: int = 1, scale_format: Incomplete | None = None, scale_color=..., scale_font: str = 'Helvetica', scale_fontsize: int = 6, scale_fontangle: int = 45, scale_largeticks: float = 0.5, scale_ticks: int = 1, scale_smallticks: float = 0.3, scale_largetick_interval: float = 1000000.0, scale_smalltick_interval: float = 10000.0, scale_largetick_labels: int = 1, scale_smalltick_labels: int = 0, axis_labels: int = 1, start: Incomplete | None = None, end: Incomplete | None = None, greytrack_font_colour: Incomplete | None = None, scale_colour: Incomplete | None = None) -> None: ...
    def add_set(self, set) -> None: ...
    def new_set(self, type: str = 'feature', **args): ...
    def del_set(self, set_id) -> None: ...
    def get_sets(self): ...
    def get_ids(self): ...
    def range(self): ...
    def to_string(self, verbose: int = 0): ...
    def __getitem__(self, key): ...
