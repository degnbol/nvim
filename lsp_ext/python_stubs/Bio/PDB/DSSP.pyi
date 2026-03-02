from Bio.Data.PDBData import protein_letters_3to1 as protein_letters_3to1, residue_sasa_scales as residue_sasa_scales
from Bio.PDB.AbstractPropertyMap import AbstractResiduePropertyMap as AbstractResiduePropertyMap
from Bio.PDB.MMCIF2Dict import MMCIF2Dict as MMCIF2Dict
from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.PDBParser import PDBParser as PDBParser
from _typeshed import Incomplete

residue_max_acc = residue_sasa_scales

def version(version_string): ...
def ss_to_index(ss): ...
def dssp_dict_from_pdb_file(in_file, DSSP: str = 'dssp', dssp_version: str = '3.9.9'): ...
def make_dssp_dict(filename): ...

class DSSP(AbstractResiduePropertyMap):
    residue_max_acc: Incomplete
    def __init__(self, model, in_file, dssp: str = 'dssp', acc_array: str = 'Sander', file_type: str = '') -> None: ...
