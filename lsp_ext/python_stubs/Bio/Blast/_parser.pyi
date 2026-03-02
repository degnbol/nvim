from Bio import Entrez as Entrez
from Bio.Align import Alignment as Alignment
from Bio.Blast import HSP as HSP, Hit as Hit, Record as Record
from Bio.Seq import Seq as Seq, reverse_complement as reverse_complement
from Bio.SeqFeature import SeqFeature as SeqFeature, SimpleLocation as SimpleLocation
from Bio.SeqRecord import SeqRecord as SeqRecord
from _typeshed import Incomplete

class DTDHandler:
    parser: Incomplete
    start_methods: Incomplete
    end_methods: Incomplete
    def __init__(self) -> None: ...
    def parseFile(self, filename) -> None: ...

class SchemaHandler:
    parser: Incomplete
    start_methods: Incomplete
    end_methods: Incomplete
    def __init__(self, parser) -> None: ...

class _HSP_cache: ...

class XMLHandler:
    schema_namespace: str
    def __init__(self, parser) -> None: ...
    def __iter__(self): ...
