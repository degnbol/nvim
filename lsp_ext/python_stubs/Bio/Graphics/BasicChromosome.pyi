from _typeshed import Incomplete
from reportlab.graphics.widgetbase import Widget

class _ChromosomeComponent(Widget):
    def __init__(self) -> None: ...
    def add(self, component) -> None: ...
    def remove(self, component) -> None: ...
    def draw(self) -> None: ...

class Organism(_ChromosomeComponent):
    page_size: Incomplete
    title_size: int
    output_format: Incomplete
    def __init__(self, output_format: str = 'pdf') -> None: ...
    def draw(self, output_file, title): ...

class Chromosome(_ChromosomeComponent):
    start_x_position: int
    end_x_position: int
    start_y_position: int
    end_y_position: int
    title_size: int
    scale_num: Incomplete
    label_size: int
    chr_percent: float
    label_sep_percent: Incomplete
    def __init__(self, chromosome_name) -> None: ...
    def subcomponent_size(self): ...
    def draw(self, cur_drawing) -> None: ...

class ChromosomeSegment(_ChromosomeComponent):
    start_x_position: int
    end_x_position: int
    start_y_position: int
    end_y_position: int
    scale: int
    fill_color: Incomplete
    label: Incomplete
    label_size: int
    chr_percent: float
    def __init__(self) -> None: ...
    def draw(self, cur_drawing) -> None: ...

class AnnotatedChromosomeSegment(ChromosomeSegment):
    bp_length: Incomplete
    features: Incomplete
    default_feature_color: Incomplete
    name_qualifiers: Incomplete
    label_sep_percent: Incomplete
    def __init__(self, bp_length, features, default_feature_color=..., name_qualifiers=('gene', 'label', 'name', 'locus_tag', 'product')) -> None: ...

class TelomereSegment(ChromosomeSegment):
    def __init__(self, inverted: int = 0) -> None: ...

class SpacerSegment(ChromosomeSegment):
    def draw(self, cur_diagram) -> None: ...
