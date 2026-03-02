from .h5_importer import assign_by_field_names as assign_by_field_names, h5py as h5py, requires_h5py as requires_h5py
from _typeshed import Incomplete
from pyrosetta.rosetta.protocols.indexed_structure_store import StructureStoreManager as StructureStoreManager, StructureStoreProvider as StructureStoreProvider

class H5StructureStoreProvider(StructureStoreProvider):
    instance: Incomplete
    @requires_h5py
    def can_load(self, store_path): ...
    def can_write(self, store_path): ...
    last_residues: Incomplete
    @requires_h5py
    def load_store(self, store_path): ...

def init_H5StructureStoreProvider() -> None: ...
