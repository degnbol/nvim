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
    def __init__(self, host: Incomplete | None = None, name: Incomplete | None = None, user: Incomplete | None = None, password: Incomplete | None = None) -> None: ...
    def fetch_motif_by_id(self, id): ...
    def fetch_motifs_by_name(self, name): ...
    def fetch_motifs(self, collection=..., tf_name: Incomplete | None = None, tf_class: Incomplete | None = None, tf_family: Incomplete | None = None, matrix_id: Incomplete | None = None, tax_group: Incomplete | None = None, species: Incomplete | None = None, pazar_id: Incomplete | None = None, data_type: Incomplete | None = None, medline: Incomplete | None = None, min_ic: int = 0, min_length: int = 0, min_sites: int = 0, all: bool = False, all_versions: bool = False): ...
