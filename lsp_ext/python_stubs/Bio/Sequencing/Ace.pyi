from _typeshed import Incomplete
from collections.abc import Generator

class rd:
    name: str
    padded_bases: Incomplete
    info_items: Incomplete
    read_tags: Incomplete
    sequence: str
    def __init__(self) -> None: ...

class qa:
    qual_clipping_start: Incomplete
    qual_clipping_end: Incomplete
    align_clipping_start: Incomplete
    align_clipping_end: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class ds:
    chromat_file: str
    phd_file: str
    time: str
    chem: str
    dye: str
    template: str
    direction: str
    def __init__(self, line: Incomplete | None = None) -> None: ...

class af:
    name: str
    coru: Incomplete
    padded_start: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class bs:
    name: str
    padded_start: Incomplete
    padded_end: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class rt:
    name: str
    tag_type: str
    program: str
    padded_start: Incomplete
    padded_end: Incomplete
    date: str
    comment: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class ct:
    name: str
    tag_type: str
    program: str
    padded_start: Incomplete
    padded_end: Incomplete
    date: str
    notrans: str
    info: Incomplete
    comment: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class wa:
    tag_type: str
    program: str
    date: str
    info: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class wr:
    name: str
    aligned: str
    program: str
    date: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class Reads:
    rd: Incomplete
    qa: Incomplete
    ds: Incomplete
    rt: Incomplete
    wr: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

class Contig:
    name: str
    nbases: Incomplete
    nreads: Incomplete
    nsegments: Incomplete
    uorc: Incomplete
    sequence: str
    quality: Incomplete
    af: Incomplete
    bs: Incomplete
    reads: Incomplete
    ct: Incomplete
    wa: Incomplete
    def __init__(self, line: Incomplete | None = None) -> None: ...

def parse(source) -> Generator[Incomplete, Incomplete]: ...

class ACEFileRecord:
    ncontigs: Incomplete
    nreads: Incomplete
    contigs: Incomplete
    wa: Incomplete
    def __init__(self) -> None: ...
    def sort(self) -> None: ...

def read(handle): ...
