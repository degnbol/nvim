from Bio import MissingPythonDependencyError as MissingPythonDependencyError
from Bio.KEGG.KGML.KGML_pathway import Pathway as Pathway
from _typeshed import Incomplete

def darken(color, factor: float = 0.7): ...
def color_to_reportlab(color): ...
def get_temp_imagefilename(url): ...

class KGMLCanvas:
    pathway: Incomplete
    show_maps: Incomplete
    show_orthologs: Incomplete
    show_compounds: Incomplete
    show_genes: Incomplete
    show_reaction_entries: Incomplete
    label_compounds: Incomplete
    label_orthologs: Incomplete
    label_reaction_entries: Incomplete
    label_maps: Incomplete
    fontname: Incomplete
    fontsize: Incomplete
    draw_relations: Incomplete
    non_reactant_transparency: float
    import_imagemap: Incomplete
    margins: Incomplete
    def __init__(self, pathway, import_imagemap: bool = False, label_compounds: bool = True, label_orthologs: bool = True, label_reaction_entries: bool = True, label_maps: bool = True, show_maps: bool = False, fontname: str = 'Helvetica', fontsize: int = 6, draw_relations: bool = True, show_orthologs: bool = True, show_compounds: bool = True, show_genes: bool = True, show_reaction_entries: bool = True, margins=(0.02, 0.02)) -> None: ...
    drawing: Incomplete
    def draw(self, filename) -> None: ...
