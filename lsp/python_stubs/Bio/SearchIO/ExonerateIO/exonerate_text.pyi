from ._base import _BaseExonerateIndexer, _BaseExonerateParser

__all__ = ['ExonerateTextParser', 'ExonerateTextIndexer']

class ExonerateTextParser(_BaseExonerateParser):
    def parse_alignment_block(self, header): ...

class ExonerateTextIndexer(_BaseExonerateIndexer):
    def get_qresult_id(self, pos): ...
    def get_raw(self, offset): ...
