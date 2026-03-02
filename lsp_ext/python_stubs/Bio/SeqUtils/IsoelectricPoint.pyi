from _typeshed import Incomplete

positive_pKs: Incomplete
negative_pKs: Incomplete
pKcterminal: Incomplete
pKnterminal: Incomplete
charged_aas: Incomplete

class IsoelectricPoint:
    sequence: Incomplete
    charged_aas_content: Incomplete
    def __init__(self, protein_sequence, aa_content: Incomplete | None = None) -> None: ...
    def charge_at_pH(self, pH): ...
    def pi(self, pH: float = 7.775, min_: float = 4.05, max_: int = 12): ...
