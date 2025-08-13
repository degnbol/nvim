from Bio import SeqIO as SeqIO
from Bio.Data.PDBData import protein_letters_1to3 as protein_letters_1to3
from Bio.File import as_handle as as_handle
from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.Residue import Residue as Residue
from Bio.PDB.Structure import Structure as Structure
from Bio.PDB.StructureBuilder import StructureBuilder as StructureBuilder
from Bio.PDB.ic_data import dihedra_primary_defaults as dihedra_primary_defaults, dihedra_secondary_defaults as dihedra_secondary_defaults, dihedra_secondary_xoxt_defaults as dihedra_secondary_xoxt_defaults, hedra_defaults as hedra_defaults, ic_data_backbone as ic_data_backbone, ic_data_sidechains as ic_data_sidechains
from Bio.PDB.internal_coords import AtomKey as AtomKey, Dihedron as Dihedron, Edron as Edron, Hedron as Hedron, IC_Chain as IC_Chain, IC_Residue as IC_Residue
from _typeshed import Incomplete
from typing import TextIO

def read_PIC(file: TextIO, verbose: bool = False, quick: bool = False, defaults: bool = False) -> Structure: ...
def read_PIC_seq(seqRec: SeqIO.SeqRecord, pdbid: str = None, title: str = None, chain: str = None) -> Structure: ...
def enumerate_atoms(entity) -> None: ...
def pdb_date(datestr: str) -> str: ...
def write_PIC(entity, file, pdbid: Incomplete | None = None, chainid: Incomplete | None = None, picFlags: int = ..., hCut: float | None | None = None, pCut: float | None | None = None): ...
