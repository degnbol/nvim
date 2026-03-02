from Bio.SeqFeature import Position as Position, SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from _typeshed import Incomplete
from collections.abc import Generator

class SwissProtParserError(ValueError):
    line: Incomplete
    def __init__(self, *args, line: Incomplete | None = None) -> None: ...

class Record:
    entry_name: Incomplete
    data_class: Incomplete
    molecule_type: Incomplete
    sequence_length: Incomplete
    accessions: Incomplete
    created: Incomplete
    sequence_update: Incomplete
    annotation_update: Incomplete
    description: Incomplete
    gene_name: Incomplete
    organism: Incomplete
    organelle: str
    organism_classification: Incomplete
    taxonomy_id: Incomplete
    host_organism: Incomplete
    host_taxonomy_id: Incomplete
    references: Incomplete
    comments: Incomplete
    cross_references: Incomplete
    keywords: Incomplete
    features: Incomplete
    protein_existence: str
    seqinfo: Incomplete
    sequence: str
    def __init__(self) -> None: ...

class Reference:
    number: Incomplete
    positions: Incomplete
    comments: Incomplete
    references: Incomplete
    authors: Incomplete
    title: Incomplete
    location: Incomplete
    def __init__(self) -> None: ...

class FeatureTable(SeqFeature): ...

def parse(source) -> Generator[Incomplete]: ...
def read(source): ...
