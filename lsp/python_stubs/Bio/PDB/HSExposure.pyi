from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap as AbstractPropertyMap
from Bio.PDB.Polypeptide import CaPPBuilder as CaPPBuilder, is_aa as is_aa
from Bio.PDB.vectors import rotaxis as rotaxis
from _typeshed import Incomplete

class _AbstractHSExposure(AbstractPropertyMap):
    ca_cb_list: Incomplete
    def __init__(self, model, radius, offset, hse_up_key, hse_down_key, angle_key: Incomplete | None = None) -> None: ...

class HSExposureCA(_AbstractHSExposure):
    def __init__(self, model, radius: int = 12, offset: int = 0) -> None: ...
    def pcb_vectors_pymol(self, filename: str = 'hs_exp.py') -> None: ...

class HSExposureCB(_AbstractHSExposure):
    def __init__(self, model, radius: int = 12, offset: int = 0) -> None: ...

class ExposureCN(AbstractPropertyMap):
    def __init__(self, model, radius: float = 12.0, offset: int = 0) -> None: ...
