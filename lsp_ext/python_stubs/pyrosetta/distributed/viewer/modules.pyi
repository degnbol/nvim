import pyrosetta.distributed
from _typeshed import Incomplete
from pyrosetta.rosetta.core.pose import Pose as Pose
from pyrosetta.rosetta.core.pose.full_model_info import get_chains_from_pdb_info as get_chains_from_pdb_info, get_res_num_from_pdb_info as get_res_num_from_pdb_info
from pyrosetta.rosetta.core.select import get_residues_from_subset as get_residues_from_subset
from pyrosetta.rosetta.core.select.residue_selector import ResidueSelector as ResidueSelector, TrueResidueSelector as TrueResidueSelector

class setBackgroundColor:
    color: Incomplete
    def __init__(self, color: int = 4294967295) -> None: ...
    def apply(self, viewer, pose, pdbstring): ...

class setDisulfides:
    color: Incomplete
    radius: Incomplete
    def __init__(self, color: str = 'gold', radius: float = 0.5) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class setHydrogenBonds:
    color: Incomplete
    dashed: Incomplete
    radius: Incomplete
    def __init__(self, color: str = 'black', dashed: bool = True, radius=None) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class setHydrogens:
    color: Incomplete
    radius: Incomplete
    polar_only: Incomplete
    def __init__(self, color: str = 'white', radius: float = 0.05, polar_only: bool = False) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class setStyle:
    residue_selector: Incomplete
    cartoon: Incomplete
    cartoon_color: Incomplete
    style: Incomplete
    colorscheme: Incomplete
    radius: Incomplete
    label: Incomplete
    label_fontsize: Incomplete
    label_background: Incomplete
    label_fontcolor: Incomplete
    command: Incomplete
    def __init__(self, residue_selector=None, cartoon: bool = True, cartoon_color: str = 'spectrum', style: str = 'stick', colorscheme: str = 'blackCarbon', radius: str = '0.1', label: bool = True, label_fontsize: int = 12, label_background: bool = False, label_fontcolor: str = 'black', command=None) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class setSurface:
    residue_selector: Incomplete
    surface_type: Incomplete
    opacity: Incomplete
    color: Incomplete
    colorscheme: Incomplete
    def __init__(self, residue_selector=None, surface_type: str = 'VDW', opacity: float = 0.5, color=None, colorscheme=None) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class setZoom:
    factor: Incomplete
    def __init__(self, factor: int = 2) -> None: ...
    def apply(self, viewer, pose, pdbstring): ...

class setZoomTo:
    residue_selector: Incomplete
    def __init__(self, residue_selector=None) -> None: ...
    @pyrosetta.distributed.requires_init
    def apply(self, viewer, pose, pdbstring): ...

class ViewerInputError(Exception):
    def __init__(self, obj) -> None: ...
