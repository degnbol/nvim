from _typeshed import Incomplete

class DistributionPage:
    distributions: Incomplete
    number_of_columns: int
    page_size: Incomplete
    title_size: int
    output_format: Incomplete
    def __init__(self, output_format: str = 'pdf') -> None: ...
    def draw(self, output_file, title): ...

class BarChartDistribution:
    display_info: Incomplete
    x_axis_title: str
    y_axis_title: str
    chart_title: str
    chart_title_size: int
    padding_percent: float
    def __init__(self, display_info: Incomplete | None = None) -> None: ...
    def draw(self, cur_drawing, start_x, start_y, end_x, end_y) -> None: ...

class LineDistribution:
    def __init__(self) -> None: ...
    def draw(self, cur_drawing, start_x, start_y, end_x, end_y) -> None: ...
