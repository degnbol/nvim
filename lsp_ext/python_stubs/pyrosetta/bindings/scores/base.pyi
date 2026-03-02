from _typeshed import Incomplete
from pyrosetta.bindings.scores.serialization import PoseScoreSerializer as PoseScoreSerializer
from pyrosetta.rosetta.core.scoring import ScoreType as ScoreType
from pyrosetta.rosetta.core.simple_metrics import clear_sm_data as clear_sm_data, get_sm_data as get_sm_data
from pyrosetta.rosetta.core.simple_metrics.metrics import CustomRealValueMetric as CustomRealValueMetric, CustomStringValueMetric as CustomStringValueMetric

class ClobberWarning(UserWarning):
    def __init__(self, msg) -> None: ...

class PoseCacheAccessorBase(PoseScoreSerializer):
    pose: Incomplete
    custom_real_value_metric: Incomplete
    custom_string_value_metric: Incomplete
    def __init__(self, pose) -> None: ...
    def apply(self, metric, key, value) -> None: ...
    def __len__(self) -> int: ...
    def __iter__(self): ...
