import pyrosetta.rosetta.ObjexxFCL
import pyrosetta.rosetta.core.pack.annealer
from typing import overload

class FlexbbSimAnnealer(pyrosetta.rosetta.core.pack.annealer.SimAnnealerBase):
    def __init__(self, bestrotamer_at_seqpos: pyrosetta.rosetta.ObjexxFCL.FArray1D_int_t, bestenergy: float, start_with_current: bool, ig, rotsets, current_rot_index: pyrosetta.rosetta.ObjexxFCL.FArray1D_int_t, calc_rot_freq: bool, rot_freq: pyrosetta.rosetta.ObjexxFCL.FArray1D_float_t) -> None: ...
    @overload
    def run(self) -> None: ...
    @overload
    def run(self) -> void: ...
