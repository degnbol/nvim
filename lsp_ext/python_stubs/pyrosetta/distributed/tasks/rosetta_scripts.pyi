import pyrosetta.distributed
import pyrosetta.distributed.tasks.taskbase as taskbase
from _typeshed import Incomplete
from collections.abc import Generator

def validate(protocol_xml) -> None: ...

class BaseRosettaScriptsTask(taskbase.TaskBase):
    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def parser(self): ...
    protocol_xml: Incomplete
    def __init__(self, protocol_xml) -> None: ...
    default_options: Incomplete
    tag: Incomplete
    protocol_lock: Incomplete
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def setup(self) -> None: ...
    @property
    @pyrosetta.distributed.requires_init
    @pyrosetta.distributed.with_lock
    def parsed_protocol(self): ...
    def execute(self, pack_or_pose): ...

class MultioutputRosettaScriptsTask(BaseRosettaScriptsTask):
    @pyrosetta.distributed.requires_init
    def apply(self, pack_or_pose) -> Generator[Incomplete]: ...

class SingleoutputRosettaScriptsTask(BaseRosettaScriptsTask):
    @pyrosetta.distributed.requires_init
    def apply(self, pack_or_pose): ...
