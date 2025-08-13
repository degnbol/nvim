from _typeshed import Incomplete

class InputRecord:
    primer_info: Incomplete
    def __init__(self) -> None: ...
    def add_primer_set(self, primer_name, first_primer_seq, second_primer_seq) -> None: ...

class OutputRecord:
    amplifiers: Incomplete
    def __init__(self) -> None: ...

class Amplifier:
    hit_info: str
    length: int
    def __init__(self) -> None: ...

def read(handle): ...
