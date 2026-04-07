from Bio import BiopythonWarning as BiopythonWarning, MissingPythonDependencyError as MissingPythonDependencyError
from Bio.motifs import jaspar as jaspar, matrix as matrix
from _typeshed import Incomplete

JASPAR_DFLT_COLLECTION: str

class JASPAR5:
    name: Incomplete
    host: Incomplete
    user: Incomplete
    password: Incomplete
    dbh: Incomplete
    def __init__(self, host=None, name=None, user=None, password=None) -> None: ...
    def fetch_motif_by_id(self, id): ...
    def fetch_motifs_by_name(self, name): ...
    def fetch_motifs(self, collection=..., tf_name=None, tf_class=None, tf_family=None, matrix_id=None, tax_group=None, species=None, pazar_id=None, data_type=None, medline=None, min_ic: int = 0, min_length: int = 0, min_sites: int = 0, all: bool = False, all_versions: bool = False): ...
