from ._base import _BaseExonerateIndexer, _BaseExonerateParser

__all__ = ['ExonerateVulgarParser', 'ExonerateVulgarIndexer']

class ExonerateVulgarParser(_BaseExonerateParser):
    def parse_alignment_block(self, header): ...

class ExonerateVulgarIndexer(_BaseExonerateIndexer):
    def get_qresult_id(self, pos): ...
    def get_raw(self, offset): ...
