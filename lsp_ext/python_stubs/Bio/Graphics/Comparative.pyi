from _typeshed import Incomplete

class ComparativeScatterPlot:
    number_of_columns: int
    page_size: Incomplete
    title_size: int
    output_format: Incomplete
    display_info: Incomplete
    color_choices: Incomplete
    shape_choices: Incomplete
    def __init__(self, output_format: str = 'pdf') -> None: ...
    def draw_to_file(self, output_file, title): ...
