from Bio.Align import MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Seq import Seq as Seq
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete
from collections.abc import Generator
from xml.sax.handler import ContentHandler

def fmt_(value, format_spec: str = '%s', default_str: str = '<unknown>'): ...

class Header:
    application: str
    version: str
    date: str
    reference: str
    query: str
    query_letters: Incomplete
    database: str
    database_sequences: Incomplete
    database_letters: Incomplete
    def __init__(self) -> None: ...

class Description:
    title: str
    score: Incomplete
    bits: Incomplete
    e: Incomplete
    num_alignments: Incomplete
    def __init__(self) -> None: ...

class DescriptionExt(Description):
    items: Incomplete
    def __init__(self) -> None: ...
    title: Incomplete
    def append_item(self, item) -> None: ...

class DescriptionExtItem:
    id: Incomplete
    title: Incomplete
    accession: Incomplete
    taxid: Incomplete
    sciname: Incomplete
    def __init__(self) -> None: ...

class Alignment:
    title: str
    hit_id: str
    hit_def: str
    length: Incomplete
    hsps: Incomplete
    def __init__(self) -> None: ...

class HSP:
    score: Incomplete
    bits: Incomplete
    expect: Incomplete
    num_alignments: Incomplete
    identities: Incomplete
    positives: Incomplete
    gaps: Incomplete
    align_length: Incomplete
    strand: Incomplete
    frame: Incomplete
    query: str
    query_start: Incomplete
    query_end: Incomplete
    match: str
    sbjct: str
    sbjct_start: Incomplete
    sbjct_end: Incomplete
    def __init__(self) -> None: ...

class MultipleAlignment:
    alignment: Incomplete
    def __init__(self) -> None: ...
    def to_generic(self): ...

class Round:
    number: Incomplete
    reused_seqs: Incomplete
    new_seqs: Incomplete
    alignments: Incomplete
    multiple_alignment: Incomplete
    def __init__(self) -> None: ...

class DatabaseReport:
    database_name: Incomplete
    posted_date: Incomplete
    num_letters_in_database: Incomplete
    num_sequences_in_database: Incomplete
    ka_params: Incomplete
    gapped: int
    ka_params_gap: Incomplete
    def __init__(self) -> None: ...

class Parameters:
    matrix: str
    gap_penalties: Incomplete
    sc_match: Incomplete
    sc_mismatch: Incomplete
    num_hits: Incomplete
    num_sequences: Incomplete
    num_good_extends: Incomplete
    num_seqs_better_e: Incomplete
    hsps_no_gap: Incomplete
    hsps_prelim_gapped: Incomplete
    hsps_prelim_gapped_attemped: Incomplete
    hsps_gapped: Incomplete
    query_id: Incomplete
    query_length: Incomplete
    database_length: Incomplete
    effective_hsp_length: Incomplete
    effective_query_length: Incomplete
    effective_database_length: Incomplete
    effective_search_space: Incomplete
    effective_search_space_used: Incomplete
    frameshift: Incomplete
    threshold: Incomplete
    window_size: Incomplete
    dropoff_1st_pass: Incomplete
    gap_x_dropoff: Incomplete
    gap_x_dropoff_final: Incomplete
    gap_trigger: Incomplete
    blast_cutoff: Incomplete
    def __init__(self) -> None: ...

class Blast(Header, DatabaseReport, Parameters):
    descriptions: Incomplete
    alignments: Incomplete
    multiple_alignment: Incomplete
    def __init__(self) -> None: ...

class PSIBlast(Header, DatabaseReport, Parameters):
    rounds: Incomplete
    converged: int
    def __init__(self) -> None: ...

class _XMLparser(ContentHandler):
    def __init__(self, debug: int = 0) -> None: ...
    def startElement(self, name, attr) -> None: ...
    def characters(self, ch) -> None: ...
    def endElement(self, name) -> None: ...

class BlastParser(_XMLparser):
    def __init__(self, debug: int = 0) -> None: ...
    def reset(self) -> None: ...
    def set_hit_id(self) -> None: ...
    def set_hit_def(self) -> None: ...
    def set_hit_accession(self) -> None: ...
    def set_hit_len(self) -> None: ...

def read(handle, debug: int = 0): ...
def parse(handle, debug: int = 0) -> Generator[Incomplete]: ...
