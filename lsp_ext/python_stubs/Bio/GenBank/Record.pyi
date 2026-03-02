from _typeshed import Incomplete

class Record:
    GB_LINE_LENGTH: int
    GB_BASE_INDENT: int
    GB_FEATURE_INDENT: int
    GB_INTERNAL_INDENT: int
    GB_OTHER_INTERNAL_INDENT: int
    GB_FEATURE_INTERNAL_INDENT: int
    GB_SEQUENCE_INDENT: int
    BASE_FORMAT: Incomplete
    INTERNAL_FORMAT: Incomplete
    OTHER_INTERNAL_FORMAT: Incomplete
    BASE_FEATURE_FORMAT: Incomplete
    INTERNAL_FEATURE_FORMAT: Incomplete
    SEQUENCE_FORMAT: Incomplete
    accession: Incomplete
    base_counts: str
    comment: str
    contig: str
    data_file_division: str
    date: str
    db_source: str
    dblinks: Incomplete
    definition: str
    features: Incomplete
    gi: str
    keywords: Incomplete
    locus: str
    molecule_type: str
    nid: str
    organism: str
    origin: str
    pid: str
    primary: Incomplete
    projects: Incomplete
    references: Incomplete
    residue_type: str
    segment: str
    sequence: str
    size: str
    source: str
    taxonomy: Incomplete
    topology: str
    version: str
    wgs: str
    wgs_scafld: Incomplete
    def __init__(self) -> None: ...

class Reference:
    number: str
    bases: str
    authors: str
    consrtm: str
    title: str
    journal: str
    medline_id: str
    pubmed_id: str
    remark: str
    def __init__(self) -> None: ...

class Feature:
    key: Incomplete
    location: Incomplete
    qualifiers: Incomplete
    def __init__(self, key: str = '', location: str = '') -> None: ...

class Qualifier:
    key: Incomplete
    value: Incomplete
    def __init__(self, key: str = '', value: str = '') -> None: ...
