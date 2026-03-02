from Bio.Application import AbstractCommandline as AbstractCommandline
from _typeshed import Incomplete

class _NcbibaseblastCommandline(AbstractCommandline):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class _NcbiblastCommandline(_NcbibaseblastCommandline):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class _Ncbiblast2SeqCommandline(_NcbiblastCommandline):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class _NcbiblastMain2SeqCommandline(_Ncbiblast2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class NcbiblastpCommandline(_NcbiblastMain2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'blastp', **kwargs) -> None: ...

class NcbiblastnCommandline(_NcbiblastMain2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'blastn', **kwargs) -> None: ...

class NcbiblastxCommandline(_NcbiblastMain2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'blastx', **kwargs) -> None: ...

class NcbitblastnCommandline(_NcbiblastMain2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'tblastn', **kwargs) -> None: ...

class NcbitblastxCommandline(_NcbiblastMain2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'tblastx', **kwargs) -> None: ...

class NcbipsiblastCommandline(_Ncbiblast2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'psiblast', **kwargs) -> None: ...

class NcbirpsblastCommandline(_NcbiblastCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'rpsblast', **kwargs) -> None: ...

class NcbirpstblastnCommandline(_NcbiblastCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'rpstblastn', **kwargs) -> None: ...

class NcbiblastformatterCommandline(_NcbibaseblastCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'blast_formatter', **kwargs) -> None: ...

class NcbideltablastCommandline(_Ncbiblast2SeqCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'deltablast', **kwargs) -> None: ...

class NcbimakeblastdbCommandline(AbstractCommandline):
    parameters: Incomplete
    def __init__(self, cmd: str = 'makeblastdb', **kwargs) -> None: ...
