from ._paml import Paml as Paml
from _typeshed import Incomplete

class Yn00Error(EnvironmentError): ...

class Yn00(Paml):
    ctl_file: str
    def __init__(self, alignment: Incomplete | None = None, working_dir: Incomplete | None = None, out_file: Incomplete | None = None) -> None: ...
    def write_ctl_file(self) -> None: ...
    alignment: Incomplete
    out_file: Incomplete
    def read_ctl_file(self, ctl_file) -> None: ...
    def run(self, ctl_file: Incomplete | None = None, verbose: bool = False, command: str = 'yn00', parse: bool = True): ...

def read(results_file): ...
