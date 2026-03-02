from Bio.File import as_handle as as_handle
from Bio.PDB.Atom import Atom as Atom
from Bio.PDB.Chain import Chain as Chain
from Bio.PDB.Model import Model as Model
from Bio.PDB.PDBExceptions import PDBException as PDBException
from Bio.PDB.PDBIO import PDBIO as PDBIO
from Bio.PDB.PICIO import enumerate_atoms as enumerate_atoms, pdb_date as pdb_date, read_PIC as read_PIC, write_PIC as write_PIC
from Bio.PDB.Residue import DisorderedResidue as DisorderedResidue, Residue as Residue
from Bio.PDB.Structure import Structure as Structure
from Bio.PDB.internal_coords import IC_Residue as IC_Residue
from typing import Any

def structure_rebuild_test(entity, verbose: bool = False, quick: bool = False) -> dict: ...
def report_IC(entity: Structure | Model | Chain | Residue, reportDict: dict[str, Any] = None, verbose: bool = False) -> dict[str, Any]: ...
def IC_duplicate(entity) -> Structure: ...
def compare_residues(e0: Structure | Model | Chain, e1: Structure | Model | Chain, verbose: bool = False, quick: bool = False, rtol: float = None, atol: float = None) -> dict[str, Any]: ...
def write_PDB(entity: Structure, file: str, pdbid: str = None, chainid: str = None) -> None: ...
