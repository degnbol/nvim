from .Interfaces import AlignmentIterator as AlignmentIterator, SequentialAlignmentWriter as SequentialAlignmentWriter
from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class StockholmWriter(SequentialAlignmentWriter):
    pfam_gr_mapping: Incomplete
    pfam_gc_mapping: Incomplete
    pfam_gs_mapping: Incomplete
    def write_alignment(self, alignment) -> None: ...

class StockholmIterator(AlignmentIterator):
    pfam_gr_mapping: Incomplete
    pfam_gc_mapping: Incomplete
    pfam_gs_mapping: Incomplete
    ids: Incomplete
    sequences: Incomplete
    seq_annotation: Incomplete
    seq_col_annotation: Incomplete
    def __next__(self): ...
