from Bio import BiopythonParserWarning as BiopythonParserWarning
from _typeshed import Incomplete
from collections.abc import Generator

class PlateRecord:
    id: Incomplete
    qualifiers: Incomplete
    def __init__(self, plateid, wells: Incomplete | None = None) -> None: ...
    def __getitem__(self, index): ...
    def __setitem__(self, key, value) -> None: ...
    def __delitem__(self, key) -> None: ...
    def __iter__(self): ...
    def __contains__(self, wellid) -> bool: ...
    def __len__(self) -> int: ...
    def __eq__(self, other): ...
    def __add__(self, plate): ...
    def __sub__(self, plate): ...
    def get_row(self, row) -> Generator[Incomplete, None, Incomplete]: ...
    def get_column(self, column) -> Generator[Incomplete, None, Incomplete]: ...
    def subtract_control(self, control: str = 'A01', wells: Incomplete | None = None): ...

class WellRecord:
    plate: Incomplete
    id: Incomplete
    max: Incomplete
    min: Incomplete
    average_height: Incomplete
    area: Incomplete
    plateau: Incomplete
    slope: Incomplete
    lag: Incomplete
    v: Incomplete
    y0: Incomplete
    model: Incomplete
    def __init__(self, wellid, plate: Incomplete | None = None, signals: Incomplete | None = None) -> None: ...
    def __setitem__(self, time, signal) -> None: ...
    def __getitem__(self, time): ...
    def __iter__(self): ...
    def __eq__(self, other): ...
    def __add__(self, well): ...
    def __sub__(self, well): ...
    def __len__(self) -> int: ...
    def get_raw(self): ...
    def get_times(self): ...
    def get_signals(self): ...
    def fit(self, function=('gompertz', 'logistic', 'richards')): ...

def JsonIterator(handle) -> Generator[Incomplete]: ...
def CsvIterator(handle) -> Generator[Incomplete]: ...

class JsonWriter:
    plates: Incomplete
    def __init__(self, plates) -> None: ...
    def write(self, handle): ...
