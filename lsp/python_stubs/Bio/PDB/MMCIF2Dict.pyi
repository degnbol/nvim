from Bio.File import as_handle as as_handle
from _typeshed import Incomplete

class MMCIF2Dict(dict):
    quote_chars: Incomplete
    whitespace_chars: Incomplete
    def __init__(self, filename) -> None: ...
