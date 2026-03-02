from _typeshed import Incomplete
from collections.abc import Generator

class SequenceLine:
    acc: str
    nid: str
    lid: str
    pid: str
    clone: str
    image: str
    is_image: bool
    end: str
    mgc: str
    seqtype: str
    trace: str
    text: Incomplete
    def __init__(self, text: Incomplete | None = None) -> None: ...

class ProtsimLine:
    org: str
    protgi: str
    protid: str
    pct: str
    aln: str
    text: Incomplete
    def __init__(self, text: Incomplete | None = None) -> None: ...

class STSLine:
    acc: str
    unists: str
    text: Incomplete
    def __init__(self, text: Incomplete | None = None) -> None: ...

class Record:
    ID: str
    species: str
    title: str
    symbol: str
    cytoband: str
    express: Incomplete
    restr_expr: str
    gnm_terminus: str
    gene_id: str
    locuslink: str
    homol: str
    chromosome: str
    protsim: Incomplete
    sequence: Incomplete
    sts: Incomplete
    txmap: Incomplete
    def __init__(self) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
def read(handle): ...
