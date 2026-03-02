from Bio.File import as_handle as as_handle
from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.internal_coords import IC_Chain as IC_Chain, IC_Residue as IC_Residue
from Bio.PDB.vectors import homog_scale_mtx as homog_scale_mtx
from _typeshed import Incomplete

def write_SCAD(entity, file, scale: Incomplete | None = None, pdbid: Incomplete | None = None, backboneOnly: bool = False, includeCode: bool = True, maxPeptideBond: Incomplete | None = None, start: Incomplete | None = None, fin: Incomplete | None = None, handle: str = 'protein') -> None: ...

peptide_scad: str
