from Bio.PDB.Atom import Atom as Atom
from Bio.PDB.Entity import Entity as Entity
from Bio.PDB.PDBExceptions import PDBException as PDBException
from _typeshed import Incomplete

entity_levels: Incomplete

def uniqueify(items): ...
def get_unique_parents(entity_list): ...
def unfold_entities(entity_list, target_level): ...
