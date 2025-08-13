from _typeshed import Incomplete
from collections.abc import Generator

def parse(handle) -> Generator[Incomplete]: ...
def read(handle): ...

class Record:
    name: str
    type: str
    accession: str
    created: str
    data_update: str
    info_update: str
    pdoc: str
    description: str
    pattern: str
    matrix: Incomplete
    rules: Incomplete
    prorules: Incomplete
    postprocessing: Incomplete
    nr_sp_release: str
    nr_sp_seqs: str
    nr_total: Incomplete
    nr_positive: Incomplete
    nr_unknown: Incomplete
    nr_false_pos: Incomplete
    nr_false_neg: Incomplete
    nr_partial: Incomplete
    cc_taxo_range: str
    cc_max_repeat: str
    cc_site: Incomplete
    cc_skip_flag: str
    dr_positive: Incomplete
    dr_false_neg: Incomplete
    dr_false_pos: Incomplete
    dr_potential: Incomplete
    dr_unknown: Incomplete
    pdb_structs: Incomplete
    def __init__(self) -> None: ...
