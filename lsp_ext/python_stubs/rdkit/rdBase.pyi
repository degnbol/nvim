# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing basic definitions for wrapped C++ code

"""
from __future__ import annotations
import typing
__all__: list[str] = ['AttachFileToLog', 'BlockLogs', 'CaptureErrorLog', 'DisableLog', 'EnableLog', 'LogDebugMsg', 'LogErrorMsg', 'LogInfoMsg', 'LogMessage', 'LogStatus', 'LogToCppStreams', 'LogToPythonLogger', 'LogToPythonStderr', 'LogWarningMsg', 'MatchTypeVect', 'SeedRandomNumberGenerator', 'UnsignedLong_Vect', 'VectSizeT', 'VectorOfStringVectors', 'WrapLogs', 'boostVersion', 'ostream', 'rdkitBuild', 'rdkitVersion', 'std_ostream', 'streambuf']
class BlockLogs(Boost.Python.instance):
    """
    Temporarily block logs from outputting while this instance is in scope.
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __enter__(self) -> BlockLogs:
        """
        """
    def __exit__(self, exc_type: typing.Any, exc_value: typing.Any, traceback: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
class CaptureErrorLog(Boost.Python.instance):
    """
    Captures messages from rdErrorLog while this instance is in scope.
    Can be used as a context manager. The ``messages`` property is
    accessible both inside the context and after it exits.
    Nesting is supported: inner captures shadow outer ones.
    
    Example::
    
      with rdBase.CaptureErrorLog() as capture:
          rdkit_function_that_may_fail()
      print(capture.messages)
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __enter__(self) -> CaptureErrorLog:
        """
        """
    def __exit__(self, exc_type: typing.Any, exc_value: typing.Any, traceback: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def messages(self) -> str:
        """Messages captured from rdErrorLog. (default: '')"""
class MatchTypeVect(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class UnsignedLong_Vect(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class VectSizeT(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class VectorOfStringVectors(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _listNSt3__16vectorIiNS_9allocatorIiEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _listNSt3__16vectorIjNS_9allocatorIjEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _listi(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
class _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vectNSt3__16vectorIdNS_9allocatorIdEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vectNSt3__16vectorIiNS_9allocatorIiEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vectNSt3__16vectorIjNS_9allocatorIjEEEE(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vectd(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vecti(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class _vectj(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __iter__(vect):
        ...
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class ostream(std_ostream):
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, python_file_obj: typing.Any, buffer_size: int = 0) -> None:
        ...
class std_ostream(Boost.Python.instance):
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
class streambuf(Boost.Python.instance):
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, python_file_obj: typing.Any, buffer_size: int = 0) -> None:
        ...
def AttachFileToLog(spec: str, filename: str, delay: int = 100) -> None:
    """
        Causes the log to write to a file
    
    """
def DisableLog(spec: str) -> None:
    """
    """
def EnableLog(spec: str) -> None:
    """
    """
def LogDebugMsg(msg: str) -> None:
    """
        Log a message to the RDKit debug logs
    
    """
def LogErrorMsg(msg: str) -> None:
    """
        Log a message to the RDKit error logs
    
    """
def LogInfoMsg(msg: str) -> None:
    """
        Log a message to the RDKit info logs
    
    """
def LogMessage(spec: str, msg: str) -> None:
    """
        Log a message to any rdApp.* log
    
    """
def LogStatus() -> str:
    """
    """
def LogToCppStreams() -> None:
    """
        Initialize RDKit logs with C++ streams
    
    """
def LogToPythonLogger() -> None:
    """
        Initialize RDKit logs with Python's logging module
    
    """
def LogToPythonStderr() -> None:
    """
        Initialize RDKit logs with Python's stderr stream
    
    """
def LogWarningMsg(msg: str) -> None:
    """
        Log a message to the RDKit warning logs
    
    """
def SeedRandomNumberGenerator(seed: int) -> None:
    """
        Provides a seed to the standard C random number generator
        This does not affect pure Python code, but is relevant to some of the RDKit C++ components.
    
    """
def WrapLogs() -> None:
    """
        Tee the RDKit logs to Python's stderr stream
    
    """
def _version() -> str:
    """
        Deprecated, use the constant rdkitVersion instead
    
    """
_iostreamsEnabled: bool = True
_multithreadedEnabled: bool = True
_serializationEnabled: bool = True
boostVersion: str = '1_85'
rdkitBuild: str = 'Darwin|24.6.0|UNIX|AppleClang|64-bit'
rdkitVersion: str = '2026.03.5'
