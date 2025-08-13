from _typeshed import Incomplete
from collections.abc import Generator

rxn_wrap: Incomplete
name_wrap: Incomplete
id_wrap: Incomplete
struct_wrap: Incomplete

class Record:
    entry: str
    name: Incomplete
    classname: Incomplete
    sysname: Incomplete
    reaction: Incomplete
    substrate: Incomplete
    product: Incomplete
    inhibitor: Incomplete
    cofactor: Incomplete
    effector: Incomplete
    comment: Incomplete
    pathway: Incomplete
    genes: Incomplete
    disease: Incomplete
    structures: Incomplete
    dblinks: Incomplete
    def __init__(self) -> None: ...

def parse(handle) -> Generator[Incomplete]: ...
def read(handle): ...
