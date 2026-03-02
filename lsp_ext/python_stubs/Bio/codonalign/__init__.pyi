from Bio import BiopythonExperimentalWarning as BiopythonExperimentalWarning, BiopythonWarning as BiopythonWarning
from Bio.Data import CodonTable as CodonTable
from Bio.SeqRecord import SeqRecord as SeqRecord
from Bio.codonalign.codonalignment import CodonAlignment as CodonAlignment, mktest as mktest
from Bio.codonalign.codonseq import CodonSeq as CodonSeq
from _typeshed import Incomplete

def build(pro_align, nucl_seqs, corr_dict: Incomplete | None = None, gap_char: str = '-', unknown: str = 'X', codon_table: Incomplete | None = None, complete_protein: bool = False, anchor_len: int = 10, max_score: int = 10): ...
