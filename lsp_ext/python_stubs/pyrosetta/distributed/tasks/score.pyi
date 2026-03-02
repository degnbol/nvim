import pyrosetta.distributed
import pyrosetta.distributed.tasks.taskbase as taskbase
from _typeshed import Incomplete

class ScorePoseTask(taskbase.TaskBase):
    weights: Incomplete
    patch: Incomplete
    def __init__(self, weights=None, patch=None) -> None: ...
    protocol_lock: Incomplete
    def setup(self) -> None: ...
    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def score_function(self): ...
    @pyrosetta.distributed.requires_init
    def execute(self, pack_or_pose): ...
