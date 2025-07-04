from Bio.PopGen.GenePop import get_indiv as get_indiv
from _typeshed import Incomplete

def read(fname): ...

class FileRecord:
    comment_line: str
    loci_list: Incomplete
    fname: Incomplete
    def __init__(self, fname) -> None: ...
    current_pop: int
    current_ind: int
    def start_read(self) -> None: ...
    def skip_header(self) -> None: ...
    def seek_position(self, pop, indiv) -> None: ...
    def skip_population(self): ...
    def get_individual(self): ...
    def remove_population(self, pos, fname) -> None: ...
    def remove_locus_by_position(self, pos, fname) -> None: ...
    def remove_loci_by_position(self, positions, fname) -> None: ...
    def remove_locus_by_name(self, name, fname) -> None: ...
    def remove_loci_by_name(self, names, fname) -> None: ...
