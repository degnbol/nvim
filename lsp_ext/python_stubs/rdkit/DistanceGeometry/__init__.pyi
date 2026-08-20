# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from .DistGeom import *
__all__: list[str] = ['DistGeom']

# present at runtime, absent from the generated stub:
from . import DistGeom as DistGeom
from rdkit.DistanceGeometry.DistGeom import DoTriangleSmoothing as DoTriangleSmoothing
from rdkit.DistanceGeometry.DistGeom import EmbedBoundsMatrix as EmbedBoundsMatrix
