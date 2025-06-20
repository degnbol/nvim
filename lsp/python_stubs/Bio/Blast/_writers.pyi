import abc
from Bio.Seq import UndefinedSequenceError as UndefinedSequenceError
from _typeshed import Incomplete
from abc import ABC

class BaseXMLWriter(ABC, metaclass=abc.ABCMeta):
    stream: Incomplete
    def __init__(self, stream) -> None: ...
    def write(self, records): ...

class XMLWriter(BaseXMLWriter): ...
class XML2Writer(BaseXMLWriter): ...
