from pyrosetta.bindings.utility import bind_method as bind_method
from pyrosetta.rosetta.core.pose import getPoseExtraFloatScores as getPoseExtraFloatScores, getPoseExtraStringScores as getPoseExtraStringScores

def setPoseExtraScore(pose=None, name=None, value=None) -> None: ...
