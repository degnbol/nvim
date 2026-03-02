import pyrosetta.rosetta.core.id
import pyrosetta.rosetta.core.pose
import pyrosetta.rosetta.core.scoring
import pyrosetta.rosetta.numeric
from typing import overload

class CartesianAtom:
    atom_id: pyrosetta.rosetta.core.id.AtomID
    force: pyrosetta.rosetta.numeric.xyzVector_double_t
    index: int
    mass: float
    old_force: pyrosetta.rosetta.numeric.xyzVector_double_t
    old_position: pyrosetta.rosetta.numeric.xyzVector_double_t
    old_velocity: pyrosetta.rosetta.numeric.xyzVector_double_t
    position: pyrosetta.rosetta.numeric.xyzVector_double_t
    res: int
    velocity: pyrosetta.rosetta.numeric.xyzVector_double_t
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: CartesianAtom) -> None: ...
    def assign(self) -> CartesianAtom: ...

class MD_Angle:
    angle: float
    atom_id_1: pyrosetta.rosetta.core.id.AtomID
    atom_id_2: pyrosetta.rosetta.core.id.AtomID
    atom_id_3: pyrosetta.rosetta.core.id.AtomID
    index1: int
    index2: int
    index3: int
    length: float
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: MD_Angle) -> None: ...
    def assign(self) -> MD_Angle: ...

class MD_Bond:
    atom_id_1: pyrosetta.rosetta.core.id.AtomID
    atom_id_2: pyrosetta.rosetta.core.id.AtomID
    index1: int
    index2: int
    length: float
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: MD_Bond) -> None: ...
    def assign(self) -> MD_Bond: ...

class MD_HarmonicDihedral:
    angle: float
    atom_id_1: pyrosetta.rosetta.core.id.AtomID
    atom_id_2: pyrosetta.rosetta.core.id.AtomID
    atom_id_3: pyrosetta.rosetta.core.id.AtomID
    atom_id_4: pyrosetta.rosetta.core.id.AtomID
    index1: int
    index2: int
    index3: int
    index4: int
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: MD_HarmonicDihedral) -> None: ...
    def assign(self) -> MD_HarmonicDihedral: ...

class MolecularDynamics:
    @overload
    def __init__(self, inputpose: pyrosetta.rosetta.core.pose.Pose, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction) -> None: ...
    @overload
    def __init__(self, arg0: MolecularDynamics) -> None: ...
    def assign(self) -> MolecularDynamics: ...
    def doMD(self, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction, Steps: int, startTemp: float, endTemp: float) -> None: ...
    @overload
    def doMinimising(self, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction) -> None: ...
    @overload
    def doMinimising(self, constclasscore) -> void: ...
    @overload
    def testCartesianDerivatives(self, scorefxn: pyrosetta.rosetta.core.scoring.ScoreFunction) -> None: ...
    @overload
    def testCartesianDerivatives(self, constclasscore) -> void: ...
