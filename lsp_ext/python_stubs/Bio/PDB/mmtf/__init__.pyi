from .mmtfio import MMTFIO as MMTFIO
from Bio import MissingPythonDependencyError as MissingPythonDependencyError
from Bio.PDB.mmtf.DefaultParser import StructureDecoder as StructureDecoder

def get_from_decoded(decoder): ...

class MMTFParser:
    @staticmethod
    def get_structure_from_url(pdb_id): ...
    @staticmethod
    def get_structure(file_path): ...
