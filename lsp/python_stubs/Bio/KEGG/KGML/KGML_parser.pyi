from Bio.KEGG.KGML.KGML_pathway import Component as Component, Entry as Entry, Graphics as Graphics, Pathway as Pathway, Reaction as Reaction, Relation as Relation
from _typeshed import Incomplete
from collections.abc import Generator

def read(handle): ...
def parse(handle) -> Generator[Incomplete]: ...

class KGMLParser:
    entry: Incomplete
    def __init__(self, elem) -> None: ...
    pathway: Incomplete
    def parse(self): ...
