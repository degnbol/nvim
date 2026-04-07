from Bio.File import as_handle as as_handle
from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.internal_coords import IC_Chain as IC_Chain, IC_Residue as IC_Residue
from Bio.PDB.vectors import homog_scale_mtx as homog_scale_mtx

def write_SCAD(entity, file, scale=None, pdbid=None, backboneOnly: bool = False, includeCode: bool = True, maxPeptideBond=None, start=None, fin=None, handle: str = 'protein') -> None: ...

peptide_scad: str
