from _typeshed import Incomplete
from typing import NamedTuple

class PoseScoreSerializerBase:
    @staticmethod
    def to_pickle(value): ...
    @staticmethod
    def from_pickle(value): ...
    @staticmethod
    def to_base64(value): ...
    @staticmethod
    def from_base64(value): ...
    @staticmethod
    def to_base64_pickle(value): ...
    @staticmethod
    def from_base64_pickle(value): ...
    @staticmethod
    def bool_from_str(value): ...

class PoseScoreSerializer(PoseScoreSerializerBase):

    class _CustomTypeMetric(NamedTuple):
        type: Incomplete
        prefix: Incomplete
        encode_func: Incomplete
        decode_func: Incomplete
    @staticmethod
    def maybe_encode(value): ...
    @staticmethod
    def maybe_decode(value): ...
