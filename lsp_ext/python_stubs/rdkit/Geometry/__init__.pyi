# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 A module for Geometry stuff	
 
"""
from __future__ import annotations
from rdkit import DataStructs
from rdkit.Geometry.rdGeometry import Point2D
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Geometry.rdGeometry import PointND
from rdkit.Geometry.rdGeometry import UniformGrid3D_
from rdkit.Geometry.rdGeometry import UniformRealValueGrid3D
from .rdGeometry import *
__all__: list[str] = ['DataStructs', 'Point2D', 'Point3D', 'PointND', 'UniformGrid3D_', 'UniformRealValueGrid3D', 'rdGeometry']

# present at runtime, absent from the generated stub:
from rdkit.Geometry.rdGeometry import ComputeDihedralAngle as ComputeDihedralAngle
from rdkit.Geometry.rdGeometry import ComputeGridCentroid as ComputeGridCentroid
from rdkit.Geometry.rdGeometry import ComputeSignedDihedralAngle as ComputeSignedDihedralAngle
from rdkit.Geometry.rdGeometry import FindGridTerminalPoints as FindGridTerminalPoints
from rdkit.Geometry.rdGeometry import ProtrudeDistance as ProtrudeDistance
from rdkit.Geometry.rdGeometry import TanimotoDistance as TanimotoDistance
from rdkit.Geometry.rdGeometry import TverskyIndex as TverskyIndex
from rdkit.Geometry.rdGeometry import UniformGrid3D as UniformGrid3D
from rdkit.Geometry.rdGeometry import WriteGridToFile as WriteGridToFile
from . import rdGeometry as rdGeometry
