# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing geometry objects like points, grids, etc
"""
from __future__ import annotations
import typing
__all__: list[str] = ['ComputeDihedralAngle', 'ComputeGridCentroid', 'ComputeSignedDihedralAngle', 'FindGridTerminalPoints', 'Point2D', 'Point3D', 'PointND', 'ProtrudeDistance', 'TanimotoDistance', 'TverskyIndex', 'UniformGrid3D', 'UniformGrid3D_', 'UniformRealValueGrid3D', 'WriteGridToFile']
class Point2D(Boost.Python.instance):
    """
    A class to represent a two-dimensional point
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 48
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AngleTo(self, other: Point2D) -> float:
        """
            determines the angle between a vector to this point (between 0 and PI)
        
        """
    def DirectionVector(self, other: Point2D) -> Point2D:
        """
            return a normalized direction vector from this point to another
        
        """
    def DotProduct(self, other: Point2D) -> float:
        """
            Dot product with another point
        
        """
    def Length(self) -> float:
        """
            Length of the vector
        
        """
    def LengthSq(self) -> float:
        """
            Square of the length
        
        """
    def Normalize(self) -> None:
        """
            Normalize the vector (using L2 norm)
        
        """
    def SignedAngleTo(self, other: Point2D) -> float:
        """
            determines the signed angle between a vector to this point (between 0 and 2*PI)
        
        """
    def __add__(self, other: Point2D) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, idx: int) -> float:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: Point2D) -> typing.Any:
        """
        """
    def __idiv__(self, scale: float) -> Point2D:
        """
            Scalar division
        
        """
    def __imul__(self, scale: float) -> Point2D:
        """
            Scalar multiplication
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, xv: float, yv: float) -> None:
        ...
    @typing.overload
    def __init__(self, other: Point3D) -> None:
        ...
    def __isub__(self, other: Point2D) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __mul__(self, other: float) -> typing.Any:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: Point2D) -> typing.Any:
        """
        """
    def __truediv__(self, other: float) -> typing.Any:
        """
        """
    @property
    def x(self) -> float:
        """default: 0.0"""
    @x.setter
    def x(self, value: float) -> None: ...
    @property
    def y(self) -> float:
        """default: 0.0"""
    @y.setter
    def y(self, value: float) -> None: ...
class Point3D(Boost.Python.instance):
    """
    A class to represent a three-dimensional point
    The x, y, and z coordinates can be read and written using either attributes
    (i.e. pt.x = 4) or indexing (i.e. pt[0] = 4).
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 56
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AngleTo(self, other: Point3D) -> float:
        """
            determines the angle between a vector to this point (between 0 and PI)
        
        """
    def CrossProduct(self, other: Point3D) -> Point3D:
        """
            Get the cross product between two points
        
        """
    def DirectionVector(self, other: Point3D) -> Point3D:
        """
            return a normalized direction vector from this point to another
        
        """
    def Distance(self, pt2: Point3D) -> float:
        """
            Distance from this point to another point
        
        """
    def DotProduct(self, other: Point3D) -> float:
        """
            Dot product with another point
        
        """
    def Length(self) -> float:
        """
            Length of the vector
        
        """
    def LengthSq(self) -> float:
        """
            Square of the length
        
        """
    def Normalize(self) -> None:
        """
            Normalize the vector (using L2 norm)
        
        """
    def SignedAngleTo(self, other: Point3D) -> float:
        """
            determines the signed angle between a vector to this point (between 0 and 2*PI)
        
        """
    def __add__(self, other: Point3D) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, idx: int) -> float:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: Point3D) -> Point3D:
        """
            Addition to another point
        
        """
    @typing.overload
    def __iadd__(self, other: Point3D) -> typing.Any:
        """
        """
    def __idiv__(self, scale: float) -> Point3D:
        """
            Scalar division
        
        """
    def __imul__(self, scale: float) -> Point3D:
        """
            Scalar multiplication
        
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, xv: float, yv: float, zv: float) -> None:
        ...
    @typing.overload
    def __isub__(self, other: Point3D) -> Point3D:
        """
            Vector difference
        
        """
    @typing.overload
    def __isub__(self, other: Point3D) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __mul__(self, other: float) -> typing.Any:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: Point3D) -> typing.Any:
        """
        """
    def __truediv__(self, other: float) -> typing.Any:
        """
        """
    @property
    def x(self) -> float:
        """default: 0.0"""
    @x.setter
    def x(self, value: float) -> None: ...
    @property
    def y(self) -> float:
        """default: 0.0"""
    @y.setter
    def y(self, value: float) -> None: ...
    @property
    def z(self) -> float:
        """default: 0.0"""
    @z.setter
    def z(self, value: float) -> None: ...
class PointND(Boost.Python.instance):
    """
    A class to represent an N-dimensional point
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 48
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AngleTo(self, other: PointND) -> float:
        """
            determines the angle between a vector to this point (between 0 and PI)
        
        """
    def DirectionVector(self, other: PointND) -> PointND:
        """
            return a normalized direction vector from this point to another
        
        """
    def Distance(self: Point3D, pt2: Point3D) -> float:
        """
            Distance from this point to another point
        
        """
    def DotProduct(self, other: PointND) -> float:
        """
            Dot product with another point
        
        """
    def Length(self) -> float:
        """
            Length of the vector
        
        """
    def LengthSq(self) -> float:
        """
            Square of the length
        
        """
    def Normalize(self) -> None:
        """
            Normalize the vector (using L2 norm)
        
        """
    def __add__(self, other: PointND) -> typing.Any:
        """
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getitem__(self, idx: int) -> float:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    @typing.overload
    def __iadd__(self, other: PointND) -> PointND:
        """
            Addition to another point
        
        """
    @typing.overload
    def __iadd__(self, other: PointND) -> typing.Any:
        """
        """
    def __idiv__(self, scale: float) -> PointND:
        """
            Scalar division
        
        """
    def __imul__(self, scale: float) -> PointND:
        """
            Scalar multiplication
        
        """
    def __init__(self, dim: int) -> None:
        ...
    @typing.overload
    def __isub__(self, other: PointND) -> PointND:
        """
            Vector difference
        
        """
    @typing.overload
    def __isub__(self, other: PointND) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __mul__(self, other: float) -> typing.Any:
        """
        """
    def __setitem__(self, idx: int, val: float) -> float:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
    def __sub__(self, other: PointND) -> typing.Any:
        """
        """
    def __truediv__(self, other: float) -> typing.Any:
        """
        """
class UniformGrid3D_(Boost.Python.instance):
    """
    Class to represent a uniform three-dimensional
        cubic grid. Each grid point can store a poisitive integer value. For the sake
        of efficiency these value can either be binary, fit in 2, 4, 8 or 16 bits
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 96
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def CompareParams(self, other: UniformGrid3D_) -> bool:
        """
            Compare the parameters between two grid object
        
        """
    def GetGridIndex(self, xi: int, yi: int, zi: int) -> int:
        """
            Get the index to the grid point with the three integer indices provided
        
        """
    def GetGridIndices(self, idx: int) -> tuple:
        """
            Returns the integer indices of the grid index provided.
        
        """
    def GetGridPointIndex(self, point: Point3D) -> int:
        """
            Get the index to the grid point closest to the specified point
        
        """
    def GetGridPointLoc(self, pointId: int) -> Point3D:
        """
            Get the location of the specified grid point
        
        """
    def GetNumX(self) -> int:
        """
            Get the number of grid points along x-axis
        
        """
    def GetNumY(self) -> int:
        """
            Get the number of grid points along y-axis
        
        """
    def GetNumZ(self) -> int:
        """
            Get the number of grid points along z-axis
        
        """
    def GetOccupancyVect(self) -> DiscreteValueVect:
        """
            Get the occupancy vector for the grid
        
        """
    def GetOffset(self) -> Point3D:
        """
            Get the location of the center of the grid
        
        """
    def GetSize(self) -> int:
        """
            Get the size of the grid (number of grid points)
        
        """
    def GetSpacing(self) -> float:
        """
            Get the grid spacing
        
        """
    def GetVal(self, id: int) -> int:
        """
            Get the value at the specified grid point
        
        """
    def GetValPoint(self, pt: Point3D) -> int:
        """
            Get the value at the closest grid point
        
        """
    def SetSphereOccupancy(self, center: Point3D, radius: float, stepSize: float, maxLayers: int = -1, ignoreOutOfBound: bool = True) -> None:
        """
            Set the occupancy on the grid for a sphere or specified radius
             and multiple layers around this sphere, with decreasing values of 
            occupancy
            
        
        """
    def SetVal(self, id: int, val: int) -> None:
        """
            Set the value at the specified grid point
        
        """
    def SetValPoint(self, pt: Point3D, val: int) -> None:
        """
            Set the value at grid point closest to the specified point
        
        """
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: UniformGrid3D_) -> typing.Any:
        """
        """
    def __iand__(self, other: UniformGrid3D_) -> typing.Any:
        """
        """
    def __init__(self, pkl: str) -> None:
        ...
    def __ior__(self, other: UniformGrid3D_) -> typing.Any:
        """
        """
    def __isub__(self, other: UniformGrid3D_) -> typing.Any:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
class UniformRealValueGrid3D(Boost.Python.instance):
    """
    Class to represent a uniform three-dimensional
        cubic grid. Each grid point can store a floating point value. 
    """
    __getstate_manages_dict__: typing.ClassVar[bool] = True
    __instance_size__: typing.ClassVar[int] = 120
    __safe_for_unpickling__: typing.ClassVar[bool] = True
    @staticmethod
    def CompareGrids(arg1: UniformRealValueGrid3D, arg2: UniformRealValueGrid3D) -> bool:
        """
            Compare the parameters and values between two grid objects.
        
        """
    @staticmethod
    def CompareParams(arg1: UniformRealValueGrid3D, arg2: UniformRealValueGrid3D) -> bool:
        """
            Compare the parameters between two grid object.
        
        """
    @staticmethod
    def CompareVectors(arg1: UniformRealValueGrid3D, arg2: UniformRealValueGrid3D) -> bool:
        """
            Compare the vector values between two grid objects.
        
        """
    @staticmethod
    def GetGridIndex(arg1: UniformRealValueGrid3D, arg2: int, arg3: int, arg4: int) -> int:
        """
            Get the index to the grid point with the three integer indices provided
        
        """
    @staticmethod
    def GetGridIndices(arg1: UniformRealValueGrid3D, arg2: int) -> tuple:
        """
            Returns the integer indices of the grid index provided.
        
        """
    @staticmethod
    def GetGridPointIndex(arg1: UniformRealValueGrid3D, arg2: Point3D) -> int:
        """
            Get the index to the grid point closest to the specified point
        
        """
    @staticmethod
    def GetGridPointLoc(arg1: UniformRealValueGrid3D, arg2: int) -> Point3D:
        """
            Get the location of the specified grid point
        
        """
    @staticmethod
    def GetNumX(arg1: UniformRealValueGrid3D) -> int:
        """
            Get the number of grid points along x-axis
        
        """
    @staticmethod
    def GetNumY(arg1: UniformRealValueGrid3D) -> int:
        """
            Get the number of grid points along y-axis
        
        """
    @staticmethod
    def GetNumZ(arg1: UniformRealValueGrid3D) -> int:
        """
            Get the number of grid points along z-axis
        
        """
    @staticmethod
    def GetOccupancyVect(arg1: UniformRealValueGrid3D) -> RealValueVect:
        """
            Get the occupancy vector for the grid
        
        """
    @staticmethod
    def GetOffset(arg1: UniformRealValueGrid3D) -> Point3D:
        """
            Get the location of the center of the grid
        
        """
    @staticmethod
    def GetSize(arg1: UniformRealValueGrid3D) -> int:
        """
            Get the size of the grid (number of grid points)
        
        """
    @staticmethod
    def GetSpacing(arg1: UniformRealValueGrid3D) -> float:
        """
            Get the grid spacing
        
        """
    @staticmethod
    def GetVal(arg1: UniformRealValueGrid3D, arg2: int) -> float:
        """
            Get the value at the specified grid point
        
        """
    @staticmethod
    def GetValPoint(arg1: UniformRealValueGrid3D, arg2: Point3D) -> float:
        """
            Get the value at the closest grid point
        
        """
    @staticmethod
    def SetVal(arg1: UniformRealValueGrid3D, arg2: int, arg3: float) -> None:
        """
            Set the value at the specified grid point
        
        """
    @staticmethod
    def SetValPoint(arg1: UniformRealValueGrid3D, arg2: Point3D, arg3: float) -> None:
        """
            Set the value at grid point closest to the specified point
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __getinitargs__(self) -> tuple:
        """
        """
    def __getstate__(self) -> tuple:
        """
        """
    def __iadd__(self, other: UniformRealValueGrid3D) -> typing.Any:
        """
        """
    def __iand__(self, other: UniformRealValueGrid3D) -> typing.Any:
        """
        """
    @typing.overload
    def __init__(self, arg1: str) -> None:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg1: UniformRealValueGrid3D) -> None:
        ...
    @typing.overload
    def __init__(self, dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, offSet: Point3D = None) -> typing.Any:
        ...
    def __ior__(self, other: UniformRealValueGrid3D) -> typing.Any:
        """
        """
    def __isub__(self, other: UniformRealValueGrid3D) -> typing.Any:
        """
        """
    def __setstate__(self, data: tuple) -> None:
        """
        """
def ComputeDihedralAngle(pt1: Point3D, pt2: Point3D, pt3: Point3D, pt4: Point3D) -> float:
    """
        calculates the dihedral angle determined by four Point3D objects
    
    """
def ComputeGridCentroid(grid: UniformGrid3D_, pt: Point3D, windowRadius: float) -> tuple:
    """
        Compute the grid point at the center of sphere around a Point3D
    
    """
def ComputeSignedDihedralAngle(pt1: Point3D, pt2: Point3D, pt3: Point3D, pt4: Point3D) -> float:
    """
        calculates the signed dihedral angle determined by four Point3D objects
    
    """
def FindGridTerminalPoints(grid: UniformGrid3D_, windowRadius: float, inclusionFraction: float) -> tuple:
    """
        Find a grid's terminal points (defined in the subshape algorithm).
    
    """
def ProtrudeDistance(grid1: UniformGrid3D_, grid2: UniformGrid3D_) -> float:
    """
        Compute the protrude distance between two grid objects
    
    """
def TanimotoDistance(grid1: UniformGrid3D_, grid2: UniformGrid3D_) -> float:
    """
        Compute the tanimoto distance between two grid objects
    
    """
def TverskyIndex(grid1: UniformGrid3D_, grid2: UniformGrid3D_, alpha: float, beta: float) -> float:
    """
        Compute the tversky index between two grid objects
    
    """
def UniformGrid3D(dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, valType: DiscreteValueType = ..., offSet: Point3D = None) -> UniformGrid3D_:
    """
        Faking the constructor
    
    """
def WriteGridToFile(grid: UniformGrid3D_, filename: str) -> None:
    """
        Write the grid to a grid file
    
    """
