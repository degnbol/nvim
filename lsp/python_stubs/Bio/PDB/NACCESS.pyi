from Bio.PDB.AbstractPropertyMap import AbstractAtomPropertyMap as AbstractAtomPropertyMap, AbstractResiduePropertyMap as AbstractResiduePropertyMap
from Bio.PDB.PDBIO import PDBIO as PDBIO
from _typeshed import Incomplete

def run_naccess(model, pdb_file, probe_size: Incomplete | None = None, z_slice: Incomplete | None = None, naccess: str = 'naccess', temp_path: str = '/tmp/'): ...
def process_rsa_data(rsa_data): ...
def process_asa_data(rsa_data): ...

class NACCESS(AbstractResiduePropertyMap):
    def __init__(self, model, pdb_file: Incomplete | None = None, naccess_binary: str = 'naccess', tmp_directory: str = '/tmp') -> None: ...

class NACCESS_atomic(AbstractAtomPropertyMap):
    naccess_atom_dict: Incomplete
    def __init__(self, model, pdb_file: Incomplete | None = None, naccess_binary: str = 'naccess', tmp_directory: str = '/tmp') -> None: ...
