from Bio.Application import AbstractCommandline as AbstractCommandline
from _typeshed import Incomplete

class _EmbossMinimalCommandLine(AbstractCommandline):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class _EmbossCommandLine(_EmbossMinimalCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: Incomplete | None = None, **kwargs) -> None: ...

class Primer3Commandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'eprimer3', **kwargs) -> None: ...

class PrimerSearchCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'primersearch', **kwargs) -> None: ...

class FDNADistCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fdnadist', **kwargs) -> None: ...

class FTreeDistCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'ftreedist', **kwargs) -> None: ...

class FNeighborCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fneighbor', **kwargs) -> None: ...

class FSeqBootCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fseqboot', **kwargs) -> None: ...

class FDNAParsCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fdnapars', **kwargs) -> None: ...

class FProtParsCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fprotpars', **kwargs) -> None: ...

class FProtDistCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fprotdist', **kwargs) -> None: ...

class FConsenseCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fconsense', **kwargs) -> None: ...

class WaterCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'water', **kwargs) -> None: ...

class NeedleCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'needle', **kwargs) -> None: ...

class NeedleallCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'needleall', **kwargs) -> None: ...

class StretcherCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'stretcher', **kwargs) -> None: ...

class FuzznucCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fuzznuc', **kwargs) -> None: ...

class FuzzproCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'fuzzpro', **kwargs) -> None: ...

class Est2GenomeCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'est2genome', **kwargs) -> None: ...

class ETandemCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'etandem', **kwargs) -> None: ...

class EInvertedCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'einverted', **kwargs) -> None: ...

class PalindromeCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'palindrome', **kwargs) -> None: ...

class TranalignCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'tranalign', **kwargs) -> None: ...

class DiffseqCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'diffseq', **kwargs) -> None: ...

class IepCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'iep', **kwargs) -> None: ...

class SeqretCommandline(_EmbossMinimalCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'seqret', **kwargs) -> None: ...

class SeqmatchallCommandline(_EmbossCommandLine):
    parameters: Incomplete
    def __init__(self, cmd: str = 'seqmatchall', **kwargs) -> None: ...
