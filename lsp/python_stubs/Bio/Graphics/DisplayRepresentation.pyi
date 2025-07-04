from Bio.Graphics.BasicChromosome import ChromosomeSegment as ChromosomeSegment, TelomereSegment as TelomereSegment
from _typeshed import Incomplete

RAINBOW_COLORS: Incomplete

class ChromosomeCounts:
    def __init__(self, segment_names, color_scheme=...) -> None: ...
    def add_count(self, segment_name, count: int = 1) -> None: ...
    def scale_segment_value(self, segment_name, scale_value: Incomplete | None = None) -> None: ...
    def add_label(self, segment_name, label) -> None: ...
    def set_scale(self, segment_name, scale) -> None: ...
    def get_segment_info(self): ...
    def fill_chromosome(self, chromosome): ...
