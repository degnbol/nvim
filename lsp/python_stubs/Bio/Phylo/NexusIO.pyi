from Bio.Nexus import Nexus as Nexus
from Bio.Phylo import Newick as Newick, NewickIO as NewickIO
from _typeshed import Incomplete
from collections.abc import Generator

NEX_TEMPLATE: str
TREE_TEMPLATE: str

def parse(handle) -> Generator[Incomplete, None, Incomplete]: ...
def write(obj, handle, **kwargs): ...
