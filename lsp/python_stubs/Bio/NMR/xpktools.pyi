from _typeshed import Incomplete

HEADERLEN: int

class XpkEntry:
    fields: Incomplete
    def __init__(self, entry, headline) -> None: ...

class Peaklist:
    firstline: Incomplete
    axislabels: Incomplete
    dataset: Incomplete
    sw: Incomplete
    sf: Incomplete
    datalabels: Incomplete
    data: Incomplete
    def __init__(self, infn) -> None: ...
    dict: Incomplete
    def residue_dict(self, index): ...
    def write_header(self, outfn) -> None: ...

def replace_entry(line, fieldn, newentry): ...
def data_table(fn_list, datalabel, keyatom): ...
