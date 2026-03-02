from Bio import BiopythonWarning as BiopythonWarning, Seq as Seq, SeqUtils as SeqUtils
from _typeshed import Incomplete

DNA_NN1: Incomplete
DNA_NN2: Incomplete
DNA_NN3: Incomplete
DNA_NN4: Incomplete
RNA_NN1: Incomplete
RNA_NN2: Incomplete
RNA_NN3: Incomplete
R_DNA_NN1: Incomplete
DNA_IMM1: Incomplete
DNA_TMM1: Incomplete
DNA_DE1: Incomplete
RNA_DE1: Incomplete

def make_table(oldtable: Incomplete | None = None, values: Incomplete | None = None): ...
def salt_correction(Na: int = 0, K: int = 0, Tris: int = 0, Mg: int = 0, dNTPs: int = 0, method: int = 1, seq: Incomplete | None = None): ...
def chem_correction(melting_temp, DMSO: int = 0, fmd: int = 0, DMSOfactor: float = 0.75, fmdfactor: float = 0.65, fmdmethod: int = 1, GC: Incomplete | None = None): ...
def Tm_Wallace(seq, check: bool = True, strict: bool = True): ...
def Tm_GC(seq, check: bool = True, strict: bool = True, valueset: int = 7, userset: Incomplete | None = None, Na: int = 50, K: int = 0, Tris: int = 0, Mg: int = 0, dNTPs: int = 0, saltcorr: int = 0, mismatch: bool = True): ...
def Tm_NN(seq, check: bool = True, strict: bool = True, c_seq: Incomplete | None = None, shift: int = 0, nn_table: Incomplete | None = None, tmm_table: Incomplete | None = None, imm_table: Incomplete | None = None, de_table: Incomplete | None = None, dnac1: int = 25, dnac2: int = 25, selfcomp: bool = False, Na: int = 50, K: int = 0, Tris: int = 0, Mg: int = 0, dNTPs: int = 0, saltcorr: int = 5): ...
