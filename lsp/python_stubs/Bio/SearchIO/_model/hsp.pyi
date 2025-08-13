from ._base import _BaseHSP
from Bio import BiopythonWarning as BiopythonWarning
from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.SearchIO._utils import allitems as allitems, fragcascade as fragcascade, fullcascade as fullcascade, getattr_str as getattr_str, singleitem as singleitem
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class HSP(_BaseHSP):
    output_index: Incomplete
    def __init__(self, fragments=(), output_index: int = -1) -> None: ...
    def __iter__(self): ...
    def __contains__(self, fragment) -> bool: ...
    def __len__(self) -> int: ...
    def __bool__(self) -> bool: ...
    def __getitem__(self, idx): ...
    def __setitem__(self, idx, fragments) -> None: ...
    def __delitem__(self, idx) -> None: ...
    aln_span: Incomplete
    hit_start: Incomplete
    query_start: Incomplete
    hit_end: Incomplete
    query_end: Incomplete
    hit_span: Incomplete
    query_span: Incomplete
    hit_range: Incomplete
    query_range: Incomplete
    hit_inter_ranges: Incomplete
    query_inter_ranges: Incomplete
    hit_inter_spans: Incomplete
    query_inter_spans: Incomplete
    is_fragmented: Incomplete
    hit_description: Incomplete
    query_description: Incomplete
    hit_id: Incomplete
    query_id: Incomplete
    molecule_type: Incomplete
    fragment: Incomplete
    hit: Incomplete
    query: Incomplete
    aln: Incomplete
    aln_annotation: Incomplete
    hit_features: Incomplete
    query_features: Incomplete
    hit_strand: Incomplete
    query_strand: Incomplete
    hit_frame: Incomplete
    query_frame: Incomplete
    fragments: Incomplete
    hit_all: Incomplete
    query_all: Incomplete
    aln_all: Incomplete
    aln_annotation_all: Incomplete
    hit_features_all: Incomplete
    query_features_all: Incomplete
    hit_strand_all: Incomplete
    query_strand_all: Incomplete
    hit_frame_all: Incomplete
    query_frame_all: Incomplete
    hit_start_all: Incomplete
    query_start_all: Incomplete
    hit_end_all: Incomplete
    query_end_all: Incomplete
    hit_span_all: Incomplete
    query_span_all: Incomplete
    hit_range_all: Incomplete
    query_range_all: Incomplete

class HSPFragment(_BaseHSP):
    aln_annotation: Incomplete
    def __init__(self, hit_id: str = '<unknown id>', query_id: str = '<unknown id>', hit: Incomplete | None = None, query: Incomplete | None = None, molecule_type: Incomplete | None = None) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...
    hit: Incomplete
    query: Incomplete
    aln: Incomplete
    molecule_type: Incomplete
    aln_span: Incomplete
    hit_description: Incomplete
    query_description: Incomplete
    hit_id: Incomplete
    query_id: Incomplete
    hit_features: Incomplete
    query_features: Incomplete
    hit_strand: Incomplete
    query_strand: Incomplete
    hit_frame: Incomplete
    query_frame: Incomplete
    hit_start: Incomplete
    query_start: Incomplete
    hit_end: Incomplete
    query_end: Incomplete
    hit_span: Incomplete
    query_span: Incomplete
    hit_range: Incomplete
    query_range: Incomplete
