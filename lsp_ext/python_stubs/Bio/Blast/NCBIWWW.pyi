from Bio import BiopythonWarning as BiopythonWarning
from Bio._utils import function_with_previous as function_with_previous
from _typeshed import Incomplete

email: Incomplete
tool: str
NCBI_BLAST_URL: str

@function_with_previous
def qblast(program, database, sequence, url_base=..., auto_format=None, composition_based_statistics=None, db_genetic_code=None, endpoints=None, entrez_query: str = '(none)', expect: float = 10.0, filter=None, gapcosts=None, genetic_code=None, hitlist_size: int = 50, i_thresh=None, layout=None, lcase_mask=None, matrix_name=None, nucl_penalty=None, nucl_reward=None, other_advanced=None, perc_ident=None, phi_pattern=None, query_file=None, query_believe_defline=None, query_from=None, query_to=None, searchsp_eff=None, service=None, threshold=None, ungapped_alignment=None, word_size=None, short_query=None, alignments: int = 500, alignment_view=None, descriptions: int = 500, entrez_links_new_window=None, expect_low=None, expect_high=None, format_entrez_query=None, format_object=None, format_type: str = 'XML', ncbi_gi=None, results_file=None, show_overview=None, megablast=None, template_type=None, template_length=None, username: str = 'blast', password=None): ...
