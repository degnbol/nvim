from ._base import _BaseExonerateParser
from .exonerate_vulgar import ExonerateVulgarIndexer

__all__ = ['ExonerateCigarParser', 'ExonerateCigarIndexer']

class ExonerateCigarParser(_BaseExonerateParser):
    def parse_alignment_block(self, header): ...

class ExonerateCigarIndexer(ExonerateVulgarIndexer):
    def get_qresult_id(self, pos): ...
