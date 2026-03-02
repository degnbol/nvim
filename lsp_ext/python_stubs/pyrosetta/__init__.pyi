import pyrosetta.rosetta as rosetta
from _typeshed import Incomplete
from pyrosetta.distributed.utility.log import LoggingContext as LoggingContext
from pyrosetta.io import Pose as Pose, create_score_function as create_score_function, dump_cif as dump_cif, dump_file as dump_file, dump_mmtf as dump_mmtf, dump_multimodel_pdb as dump_multimodel_pdb, dump_pdb as dump_pdb, dump_scored_pdb as dump_scored_pdb, get_fa_scorefxn as get_fa_scorefxn, get_score_function as get_score_function, pose_from_file as pose_from_file, pose_from_pdb as pose_from_pdb, pose_from_sequence as pose_from_sequence, poses_from_files as poses_from_files, poses_from_multimodel_pdb as poses_from_multimodel_pdb, poses_from_sequences as poses_from_sequences, poses_from_silent as poses_from_silent, poses_to_silent as poses_to_silent
from pyrosetta.rosetta.core.id import AtomID as AtomID
from pyrosetta.rosetta.core.kinematics import FoldTree as FoldTree, MoveMap as MoveMap
from pyrosetta.rosetta.core.scoring import ScoreFunction as ScoreFunction
from pyrosetta.rosetta.protocols.moves import MonteCarlo as MonteCarlo, PyMOLMover as PyMOLMover, RepeatMover as RepeatMover, SequenceMover as SequenceMover, TrialMover as TrialMover
from pyrosetta.rosetta.protocols.simple_moves import SwitchResidueTypeSetMover as SwitchResidueTypeSetMover
from pyrosetta.toolbox import PyJobDistributor as PyJobDistributor, etable_atom_pair_energies as etable_atom_pair_energies
from typing import NamedTuple

logger: Incomplete

class PyRosettaException(Exception): ...

class PythonPyExitCallback(rosetta.utility.py.PyExitCallback):
    def __init__(self) -> None: ...
    def exit_callback(self) -> None: ...

def init(options: str = '-ex1 -ex2aro', extra_options: str = '', set_logging_handler=None, notebook=None, silent: bool = False) -> None: ...
def version(): ...
def Vector1(list_in): ...
def Set(list_in): ...
def generate_nonstandard_residue_set(pose, params_list): ...
def standard_task_factory(): ...
def standard_packer_task(pose): ...

class CD(NamedTuple):
    base: Incomplete
    first: Incomplete
    last: Incomplete
    methods: Incomplete

ScoreTypesRegistry: Incomplete

def defineEnergyMethodCreator(class_, scoreType): ...

class EnergyMethod:
    scoreName: Incomplete
    scoreType: Incomplete
    version: Incomplete
    def __init__(self, scoreName=None, scoreType=None, version: int = 1) -> None: ...
    def __call__(self, original_class): ...
