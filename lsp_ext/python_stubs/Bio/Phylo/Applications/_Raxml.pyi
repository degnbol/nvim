from Bio.Application import AbstractCommandline as AbstractCommandline
from _typeshed import Incomplete

class RaxmlCommandline(AbstractCommandline):
    parameters: Incomplete
    parsimony_seed: int
    def __init__(self, cmd: str = 'raxmlHPC', **kwargs) -> None: ...
