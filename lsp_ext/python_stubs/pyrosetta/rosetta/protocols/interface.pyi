import pyrosetta.rosetta.core.pose
import pyrosetta.rosetta.utility

def select_interface_residues(pose: pyrosetta.rosetta.core.pose.Pose, interface: str, interface_distance: int) -> pyrosetta.rosetta.utility.vector1_bool: ...
