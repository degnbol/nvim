from _typeshed import Incomplete

class PrintFormat:
    ConsoleWidth: int
    NameWidth: int
    MaxSize: int
    Cmodulo: Incomplete
    PrefWidth: Incomplete
    Indent: int
    linesize: Incomplete
    def print_as(self, what: str = 'list') -> None: ...
    def format_output(self, dct, title: str = '', s1: str = ''): ...
    def print_that(self, dct, title: str = '', s1: str = '') -> None: ...
    def make_format(self, cut=(), title: str = '', nc=(), s1: str = ''): ...
