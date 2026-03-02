from .search import PairQueryExecutor as PairQueryExecutor, SingleQueryExecutor as SingleQueryExecutor, StructureDatabase as StructureDatabase, StructurePairQuery as StructurePairQuery, StructureSingleQuery as StructureSingleQuery
from _typeshed import Incomplete
from pyrosetta.utility.array import structured_array_to_basic as structured_array_to_basic

class StructureSearchManager:
    source_residues: Incomplete
    target_atoms: Incomplete
    db: Incomplete
    def __init__(self, source_residues) -> None: ...
    def pair_query(self, query_components, query_tolerance, primary_range=None): ...
    def single_query(self, query_component, query_tolerance): ...
