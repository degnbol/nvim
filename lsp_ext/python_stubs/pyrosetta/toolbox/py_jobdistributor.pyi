from _typeshed import Incomplete
from pyrosetta.bindings.scores import PoseScoreSerializer as PoseScoreSerializer
from pyrosetta.rosetta.core.io.raw_data import ScoreMap as ScoreMap

def output_scorefile(pose, pdb_name, current_name, scorefilepath, scorefxn, nstruct, native_pose=None, additional_decoy_info=None, json_format: bool = True) -> None: ...

class PyJobDistributor:
    pdb_name: Incomplete
    nstruct: Incomplete
    compress: Incomplete
    current_id: Incomplete
    current_name: Incomplete
    scorefxn: Incomplete
    native_pose: Incomplete
    additional_decoy_info: Incomplete
    json_format: bool
    sequence: Incomplete
    def __init__(self, pdb_name, nstruct, scorefxn, compress: bool = False) -> None: ...
    @property
    def job_complete(self): ...
    current_in_progress_name: Incomplete
    def start_decoy(self) -> None: ...
    def output_decoy(self, pose) -> None: ...
