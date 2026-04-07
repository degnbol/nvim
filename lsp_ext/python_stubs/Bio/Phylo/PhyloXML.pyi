from Bio import BiopythonWarning as BiopythonWarning
from Bio.Align import Alignment as Alignment, MultipleSeqAlignment as MultipleSeqAlignment
from Bio.Phylo import BaseTree as BaseTree
from Bio.Seq import Seq as Seq
from Bio.SeqFeature import SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class PhyloXMLWarning(BiopythonWarning): ...
class PhyloElement(BaseTree.TreeElement): ...

class Phyloxml(PhyloElement):
    attributes: Incomplete
    phylogenies: Incomplete
    other: Incomplete
    def __init__(self, attributes, phylogenies=None, other=None) -> None: ...
    def __getitem__(self, index): ...
    def __iter__(self): ...
    def __len__(self) -> int: ...

class Other(PhyloElement):
    tag: Incomplete
    namespace: Incomplete
    attributes: Incomplete
    value: Incomplete
    children: Incomplete
    def __init__(self, tag, namespace=None, attributes=None, value=None, children=None) -> None: ...
    def __iter__(self): ...

class Phylogeny(PhyloElement, BaseTree.Tree):
    root: Incomplete
    rooted: Incomplete
    rerootable: Incomplete
    branch_length_unit: Incomplete
    type: Incomplete
    name: Incomplete
    id: Incomplete
    description: Incomplete
    date: Incomplete
    confidences: Incomplete
    clade_relations: Incomplete
    sequence_relations: Incomplete
    properties: Incomplete
    other: Incomplete
    def __init__(self, root=None, rooted: bool = True, rerootable=None, branch_length_unit=None, type=None, name=None, id=None, description=None, date=None, confidences=None, clade_relations=None, sequence_relations=None, properties=None, other=None) -> None: ...
    @classmethod
    def from_tree(cls, tree, **kwargs): ...
    @classmethod
    def from_clade(cls, clade, **kwargs): ...
    def as_phyloxml(self): ...
    def to_phyloxml_container(self, **kwargs): ...
    def to_alignment(self): ...
    @property
    def alignment(self): ...
    confidence: Incomplete

class Clade(PhyloElement, BaseTree.Clade):
    branch_length: Incomplete
    id_source: Incomplete
    name: Incomplete
    width: Incomplete
    color: Incomplete
    node_id: Incomplete
    events: Incomplete
    binary_characters: Incomplete
    date: Incomplete
    confidences: Incomplete
    taxonomies: Incomplete
    sequences: Incomplete
    distributions: Incomplete
    references: Incomplete
    properties: Incomplete
    clades: Incomplete
    other: Incomplete
    def __init__(self, branch_length=None, id_source=None, name=None, width=None, color=None, node_id=None, events=None, binary_characters=None, date=None, confidences=None, taxonomies=None, sequences=None, distributions=None, references=None, properties=None, clades=None, other=None) -> None: ...
    @classmethod
    def from_clade(cls, clade, **kwargs): ...
    def to_phylogeny(self, **kwargs): ...
    confidence: Incomplete
    taxonomy: Incomplete

class BranchColor(PhyloElement, BaseTree.BranchColor):
    def __init__(self, *args, **kwargs) -> None: ...

class Accession(PhyloElement):
    value: Incomplete
    source: Incomplete
    def __init__(self, value, source) -> None: ...

class Annotation(PhyloElement):
    re_ref: Incomplete
    ref: Incomplete
    source: Incomplete
    evidence: Incomplete
    type: Incomplete
    desc: Incomplete
    confidence: Incomplete
    uri: Incomplete
    properties: Incomplete
    def __init__(self, ref=None, source=None, evidence=None, type=None, desc=None, confidence=None, uri=None, properties=None) -> None: ...

class BinaryCharacters(PhyloElement):
    type: Incomplete
    gained_count: Incomplete
    lost_count: Incomplete
    present_count: Incomplete
    absent_count: Incomplete
    gained: Incomplete
    lost: Incomplete
    present: Incomplete
    absent: Incomplete
    def __init__(self, type=None, gained_count=None, lost_count=None, present_count=None, absent_count=None, gained=None, lost=None, present=None, absent=None) -> None: ...

class CladeRelation(PhyloElement):
    distance: Incomplete
    type: Incomplete
    id_ref_0: Incomplete
    id_ref_1: Incomplete
    confidence: Incomplete
    def __init__(self, type, id_ref_0, id_ref_1, distance=None, confidence=None) -> None: ...

class Confidence(float, PhyloElement):
    def __new__(cls, value, type: str = 'unknown'): ...
    @property
    def value(self): ...

class Date(PhyloElement):
    value: Incomplete
    unit: Incomplete
    desc: Incomplete
    minimum: Incomplete
    maximum: Incomplete
    def __init__(self, value=None, unit=None, desc=None, minimum=None, maximum=None) -> None: ...

class Distribution(PhyloElement):
    desc: Incomplete
    points: Incomplete
    polygons: Incomplete
    def __init__(self, desc=None, points=None, polygons=None) -> None: ...

class DomainArchitecture(PhyloElement):
    length: Incomplete
    domains: Incomplete
    def __init__(self, length=None, domains=None) -> None: ...

class Events(PhyloElement):
    ok_type: Incomplete
    type: Incomplete
    duplications: Incomplete
    speciations: Incomplete
    losses: Incomplete
    confidence: Incomplete
    def __init__(self, type=None, duplications=None, speciations=None, losses=None, confidence=None) -> None: ...
    def items(self): ...
    def keys(self): ...
    def values(self): ...
    def __len__(self) -> int: ...
    def __getitem__(self, key): ...
    def __setitem__(self, key, val) -> None: ...
    def __delitem__(self, key) -> None: ...
    def __iter__(self): ...
    def __contains__(self, key) -> bool: ...

class Id(PhyloElement):
    value: Incomplete
    provider: Incomplete
    def __init__(self, value, provider=None) -> None: ...

class MolSeq(PhyloElement):
    re_value: Incomplete
    value: Incomplete
    is_aligned: Incomplete
    def __init__(self, value, is_aligned=None) -> None: ...

class Point(PhyloElement):
    geodetic_datum: Incomplete
    lat: Incomplete
    long: Incomplete
    alt: Incomplete
    alt_unit: Incomplete
    def __init__(self, geodetic_datum, lat, long, alt=None, alt_unit=None) -> None: ...

class Polygon(PhyloElement):
    points: Incomplete
    def __init__(self, points=None) -> None: ...

class Property(PhyloElement):
    re_ref: Incomplete
    ok_applies_to: Incomplete
    ok_datatype: Incomplete
    unit: Incomplete
    id_ref: Incomplete
    value: Incomplete
    ref: Incomplete
    applies_to: Incomplete
    datatype: Incomplete
    def __init__(self, value, ref, applies_to, datatype, unit=None, id_ref=None) -> None: ...

class ProteinDomain(PhyloElement):
    value: Incomplete
    start: Incomplete
    end: Incomplete
    confidence: Incomplete
    id: Incomplete
    def __init__(self, value, start, end, confidence=None, id=None) -> None: ...
    @classmethod
    def from_seqfeature(cls, feat): ...
    def to_seqfeature(self): ...

class Reference(PhyloElement):
    re_doi: Incomplete
    doi: Incomplete
    desc: Incomplete
    def __init__(self, doi=None, desc=None) -> None: ...

class Sequence(PhyloElement):
    types: Incomplete
    re_symbol: Incomplete
    type: Incomplete
    id_ref: Incomplete
    id_source: Incomplete
    symbol: Incomplete
    accession: Incomplete
    name: Incomplete
    location: Incomplete
    mol_seq: Incomplete
    uri: Incomplete
    domain_architecture: Incomplete
    annotations: Incomplete
    other: Incomplete
    def __init__(self, type=None, id_ref=None, id_source=None, symbol=None, accession=None, name=None, location=None, mol_seq=None, uri=None, domain_architecture=None, annotations=None, other=None) -> None: ...
    @classmethod
    def from_seqrecord(cls, record, is_aligned=None): ...
    def to_seqrecord(self): ...

class SequenceRelation(PhyloElement):
    ok_type: Incomplete
    distance: Incomplete
    type: Incomplete
    id_ref_0: Incomplete
    id_ref_1: Incomplete
    confidence: Incomplete
    def __init__(self, type, id_ref_0, id_ref_1, distance=None, confidence=None) -> None: ...

class Taxonomy(PhyloElement):
    re_code: Incomplete
    ok_rank: Incomplete
    id_source: Incomplete
    id: Incomplete
    code: Incomplete
    scientific_name: Incomplete
    authority: Incomplete
    rank: Incomplete
    uri: Incomplete
    common_names: Incomplete
    synonyms: Incomplete
    other: Incomplete
    def __init__(self, id_source=None, id=None, code=None, scientific_name=None, authority=None, rank=None, uri=None, common_names=None, synonyms=None, other=None) -> None: ...

class Uri(PhyloElement):
    value: Incomplete
    desc: Incomplete
    type: Incomplete
    def __init__(self, value, desc=None, type=None) -> None: ...
