import subprocess
from Bio import BiopythonDeprecationWarning as BiopythonDeprecationWarning
from _typeshed import Incomplete

class ApplicationError(subprocess.CalledProcessError):
    returncode: Incomplete
    cmd: Incomplete
    stdout: Incomplete
    stderr: Incomplete
    def __init__(self, returncode, cmd, stdout: str = '', stderr: str = '') -> None: ...

class AbstractCommandline:
    parameters: Incomplete
    program_name: Incomplete
    def __init__(self, cmd, **kwargs) -> None: ...
    def set_parameter(self, name, value: Incomplete | None = None) -> None: ...
    def __setattr__(self, name, value) -> None: ...
    def __call__(self, stdin: Incomplete | None = None, stdout: bool = True, stderr: bool = True, cwd: Incomplete | None = None, env: Incomplete | None = None): ...

class _AbstractParameter:
    def __init__(self) -> None: ...

class _Option(_AbstractParameter):
    names: Incomplete
    is_filename: Incomplete
    checker_function: Incomplete
    description: Incomplete
    equate: Incomplete
    is_required: Incomplete
    is_set: bool
    value: Incomplete
    def __init__(self, names, description, filename: bool = False, checker_function: Incomplete | None = None, is_required: bool = False, equate: bool = True) -> None: ...

class _Switch(_AbstractParameter):
    names: Incomplete
    description: Incomplete
    is_set: bool
    is_required: bool
    def __init__(self, names, description) -> None: ...

class _Argument(_AbstractParameter):
    names: Incomplete
    is_filename: Incomplete
    checker_function: Incomplete
    description: Incomplete
    is_required: Incomplete
    is_set: bool
    value: Incomplete
    def __init__(self, names, description, filename: bool = False, checker_function: Incomplete | None = None, is_required: bool = False) -> None: ...

class _ArgumentList(_Argument): ...

class _StaticArgument(_AbstractParameter):
    names: Incomplete
    is_required: bool
    is_set: bool
    value: Incomplete
    def __init__(self, value) -> None: ...
