from Bio.Data import IUPACData as IUPACData
from _typeshed import Incomplete

unambiguous_dna_by_name: Incomplete
unambiguous_dna_by_id: Incomplete
unambiguous_rna_by_name: Incomplete
unambiguous_rna_by_id: Incomplete
generic_by_name: Incomplete
generic_by_id: Incomplete
ambiguous_dna_by_name: Incomplete
ambiguous_dna_by_id: Incomplete
ambiguous_rna_by_name: Incomplete
ambiguous_rna_by_id: Incomplete
ambiguous_generic_by_name: Incomplete
ambiguous_generic_by_id: Incomplete
standard_dna_table: Incomplete
standard_rna_table: Incomplete

class TranslationError(Exception): ...

class CodonTable:
    forward_table: dict[str, str]
    back_table: dict[str, str]
    start_codons: list[str]
    stop_codons: list[str]
    nucleotide_alphabet: Incomplete
    protein_alphabet: Incomplete
    def __init__(self, nucleotide_alphabet: str | None = None, protein_alphabet: str | None = None, forward_table: dict[str, str] = ..., back_table: dict[str, str] = ..., start_codons: list[str] = ..., stop_codons: list[str] = ...) -> None: ...

def make_back_table(table, default_stop_codon): ...

class NCBICodonTable(CodonTable):
    nucleotide_alphabet: str | None
    protein_alphabet: Incomplete
    id: Incomplete
    names: Incomplete
    forward_table: Incomplete
    back_table: Incomplete
    start_codons: Incomplete
    stop_codons: Incomplete
    def __init__(self, id, names, table, start_codons, stop_codons) -> None: ...

class NCBICodonTableDNA(NCBICodonTable):
    nucleotide_alphabet: Incomplete

class NCBICodonTableRNA(NCBICodonTable):
    nucleotide_alphabet: Incomplete

class AmbiguousCodonTable(CodonTable):
    def __init__(self, codon_table, ambiguous_nucleotide_alphabet, ambiguous_nucleotide_values, ambiguous_protein_alphabet, ambiguous_protein_values) -> None: ...
    def __getattr__(self, name): ...

def list_possible_proteins(codon, forward_table, ambiguous_nucleotide_values): ...
def list_ambiguous_codons(codons, ambiguous_nucleotide_values): ...

class AmbiguousForwardTable:
    forward_table: Incomplete
    ambiguous_nucleotide: Incomplete
    ambiguous_protein: Incomplete
    def __init__(self, forward_table, ambiguous_nucleotide, ambiguous_protein) -> None: ...
    def __contains__(self, codon) -> bool: ...
    def get(self, codon, failobj: Incomplete | None = None): ...
    def __getitem__(self, codon): ...

def register_ncbi_table(name, alt_name, id, table, start_codons, stop_codons) -> None: ...
