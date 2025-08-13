from Bio.Align import Alignment as Alignment, interfaces as interfaces
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class AlignmentIterator(interfaces.AlignmentIterator):
    fmt: str
    gf_mapping: Incomplete
    gr_mapping: Incomplete
    gc_mapping: Incomplete
    gs_mapping: Incomplete

class AlignmentWriter(interfaces.AlignmentWriter):
    gf_mapping: Incomplete
    gs_mapping: Incomplete
    gr_mapping: Incomplete
    gc_mapping: Incomplete
    fmt: str
    def format_alignment(self, alignment): ...
