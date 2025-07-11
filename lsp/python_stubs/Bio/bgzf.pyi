from _typeshed import Incomplete
from collections.abc import Generator

def open(filename, mode: str = 'rb'): ...
def make_virtual_offset(block_start_offset, within_block_offset): ...
def split_virtual_offset(virtual_offset): ...
def BgzfBlocks(handle) -> Generator[Incomplete]: ...

class BgzfReader:
    max_cache: Incomplete
    def __init__(self, filename: Incomplete | None = None, mode: str = 'r', fileobj: Incomplete | None = None, max_cache: int = 100) -> None: ...
    def tell(self): ...
    def seek(self, virtual_offset): ...
    def read(self, size: int = -1): ...
    def readline(self): ...
    def __next__(self): ...
    def __iter__(self): ...
    def close(self) -> None: ...
    def seekable(self): ...
    def isatty(self): ...
    def fileno(self): ...
    def __enter__(self): ...
    def __exit__(self, type: type[BaseException] | None, value: BaseException | None, traceback: types.TracebackType | None) -> None: ...

class BgzfWriter:
    compresslevel: Incomplete
    def __init__(self, filename: Incomplete | None = None, mode: str = 'w', fileobj: Incomplete | None = None, compresslevel: int = 6) -> None: ...
    def write(self, data) -> None: ...
    def flush(self) -> None: ...
    def close(self) -> None: ...
    def tell(self): ...
    def seekable(self): ...
    def isatty(self): ...
    def fileno(self): ...
    def __enter__(self): ...
    def __exit__(self, type: type[BaseException] | None, value: BaseException | None, traceback: types.TracebackType | None) -> None: ...
