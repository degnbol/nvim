# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing Trajectory and Snapshot objects
"""
from __future__ import annotations
import typing
__all__: list[str] = ['ReadAmberTrajectory', 'ReadGromosTrajectory', 'Snapshot', 'Trajectory']
class Snapshot(Boost.Python.instance):
    """
    A class which allows storing coordinates from a trajectory
    """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetEnergy(self) -> float:
        """
            returns the energy for this Snapshot
        
        """
    def GetPoint2D(self, pointNum: int) -> Point2D:
        """
            return the coordinates at pointNum as a Point2D object; requires the Trajectory dimension to be == 2
        
        """
    def GetPoint3D(self, pointNum: int) -> Point3D:
        """
            return the coordinates at pointNum as a Point3D object; requires the Trajectory dimension to be >= 2
        
        """
    def SetEnergy(self, energy: float) -> None:
        """
            sets the energy for this Snapshot
        
        """
    @typing.overload
    def __init__(self, coordList: list, energy: float = 0.0) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, other: Snapshot) -> typing.Any:
        ...
class Trajectory(Boost.Python.instance):
    """
    A class which allows storing Snapshots from a trajectory
    """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddConformersToMol(self, mol: Mol, fromCid: int = -1, toCid: int = -1) -> int:
        """
            adds conformations from the Trajectory to mol
            fromCid is the first Snapshot that will be added as a Conformer; defaults to -1 (first available)
            toCid is the last Snapshot that will be added as a Conformer; defaults to -1 (all)
            
        
        """
    def AddSnapshot(self, s: Snapshot) -> int:
        """
            appends Snapshot s to this Trajectory; returns the zero-based index position of the added snapshot
            
        
        """
    def Clear(self) -> None:
        """
            removes all Snapshots from the Trajectory
            
        
        """
    def Dimension(self) -> int:
        """
            returns the dimensionality of this Trajectory's coordinate tuples
        
        """
    def GetSnapshot(self, snapshotNum: int) -> Snapshot:
        """
            returns the Snapshot snapshotNum, where the latter is the zero-based index of the retrieved Snapshot
            
        
        """
    def InsertSnapshot(self, snapshotNum: int, s: Snapshot) -> int:
        """
            inserts Snapshot s into the Trajectory at the position snapshotNum, where the latter is the zero-based index of the Trajectory's Snapshot before which the Snapshot s will be inserted; returns the zero-based index position of the inserted snapshot
            
        
        """
    def NumPoints(self) -> int:
        """
            returns the number of coordinate tuples associated to each Snapshot
        
        """
    def RemoveSnapshot(self, snapshotNum: int) -> int:
        """
            removes Snapshot snapshotNum from the Trajectory, where snapshotNum is the zero-based index of Snapshot to be removed
            
        
        """
    @typing.overload
    def __init__(self, dimension: int, numPoints: int, snapshotList: list = []) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, other: Trajectory) -> typing.Any:
        ...
    def __len__(self) -> int:
        """
        """
def ReadAmberTrajectory(fName: str, traj: Trajectory) -> int:
    """
        reads coordinates from an AMBER trajectory file into the Trajectory object; returns the number of Snapshot objects read in
        
    
    """
def ReadGromosTrajectory(fName: str, traj: Trajectory) -> int:
    """
        reads coordinates from a GROMOS trajectory file into the Trajectory object; returns the number of Snapshot objects read in
        
    
    """
