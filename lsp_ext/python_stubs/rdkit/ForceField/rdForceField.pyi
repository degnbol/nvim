"""
Exposes the ForceField class
"""
from __future__ import annotations
import typing
__all__: list[str] = ['ForceField', 'MMFFMolProperties']
class ForceField(Boost.Python.instance):
    """
    A force field
    """
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def AddDistanceConstraint(self, idx1: int, idx2: int, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
            Adds a distance constraint to the UFF force field (deprecated, use UFFAddDistanceConstraint instead).
        
        """
    def AddExtraPoint(self, x: float, y: float, z: float, fixed: bool = True) -> int:
        """
            Adds an extra point, this can be useful for adding constraints.
        
        """
    def AddFixedPoint(self, idx: int) -> None:
        """
            Adds a fixed point to the force field.
        
        """
    def CalcEnergy(self, pos: typing.Any = None) -> float:
        """
            Returns the energy (in kcal/mol) of the current arrangement
            or of the supplied coordinate list (if non-empty)
        
        """
    def CalcGrad(self, pos: typing.Any = None) -> typing.Any:
        """
            Returns a tuple filled with the per-coordinate gradients
            of the current arrangement or of the supplied coordinate list (if non-empty)
        
        """
    def Dimension(self) -> int:
        """
            Returns the dimension of the ForceField
        
        """
    def GetExtraPointPos(self, idx: int) -> typing.Any:
        """
            returns the location of an extra point as a tuple
        
        """
    def Initialize(self) -> None:
        """
            initializes the force field (call this before minimizing)
        
        """
    def MMFFAddAngleConstraint(self, idx1: int, idx2: int, idx3: int, relative: bool, minAngleDeg: float, maxAngleDeg: float, forceConstant: float) -> None:
        """
            Adds an angle constraint to the MMFF force field; if relative == True, then minAngleDeg and maxAngleDeg are intended as relative to the current angle.
        
        """
    def MMFFAddDistanceConstraint(self, idx1: int, idx2: int, relative: bool, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
            Adds a distance constraint to the MMFF force field; if relative == True, then minLen and maxLen are intended as relative to the current distance.
        
        """
    def MMFFAddPositionConstraint(self, idx: int, maxDispl: float, forceConstant: float) -> None:
        """
            Adds a position constraint to the MMFF force field.
        
        """
    def MMFFAddTorsionConstraint(self, idx1: int, idx2: int, idx3: int, idx4: int, relative: bool, minDihedralDeg: float, maxDihedralDeg: float, forceConstant: float) -> None:
        """
            Adds a dihedral angle constraint to the MMFF force field; if relative == True, then minDihedralDeg and maxDihedralDeg are intended as relative to the current dihedral angle.
        
        """
    def Minimize(self, maxIts: int = 200, forceTol: float = 0.0001, energyTol: float = 1e-06) -> int:
        """
            Runs some minimization iterations.
            
              Returns 0 if the minimization succeeded.
        
        """
    def MinimizeTrajectory(self, snapshotFreq: int, maxIts: int = 200, forceTol: float = 0.0001, energyTol: float = 1e-06) -> tuple:
        """
            Runs some minimization iterations, recording the minimization trajectory every snapshotFreq steps.
            
            Returns a (int, []) tuple; the int is 0 if the minimization succeeded, while the list contains Snapshot objects.
        
        """
    def NumPoints(self) -> int:
        """
            Returns the number of points the ForceField is handling
        
        """
    def Positions(self) -> typing.Any:
        """
            Returns a tuple filled with the coordinates of the
            points the ForceField is handling
        
        """
    def UFFAddAngleConstraint(self, idx1: int, idx2: int, idx3: int, relative: bool, minAngleDeg: float, maxAngleDeg: float, forceConstant: float) -> None:
        """
            Adds an angle constraint to the UFF force field; if relative == True, then minAngleDeg and maxAngleDeg are intended as relative to the current angle.
        
        """
    def UFFAddDistanceConstraint(self, idx1: int, idx2: int, relative: bool, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
            Adds a distance constraint to the UFF force field; if relative == True, then minLen and maxLen are intended as relative to the current distance.
        
        """
    def UFFAddPositionConstraint(self, idx: int, maxDispl: float, forceConstant: float) -> None:
        """
            Adds a position constraint to the UFF force field.
        
        """
    def UFFAddTorsionConstraint(self, idx1: int, idx2: int, idx3: int, idx4: int, relative: bool, minDihedralDeg: float, maxDihedralDeg: float, forceConstant: float) -> None:
        """
            Adds a dihedral angle constraint to the UFF force field; if relative == True, then minDihedralDeg and maxDihedralDeg are intended as relative to the current dihedral angle.
        
        """
class MMFFMolProperties(Boost.Python.instance):
    """
    MMFF molecular properties
    """
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def GetMMFFAngleBendParams(self, mol: typing.Any, idx1: int, idx2: int, idx3: int) -> typing.Any:
        """
            Retrieves MMFF angle bend parameters for atoms with indexes idx1, idx2, idx3 as a (angleType, ka, theta0) tuple, or None if no parameters could be found
        
        """
    def GetMMFFAngleTerm(self) -> bool:
        """
            Returns whether the angle term is included in the MMFF equation.
        
        """
    def GetMMFFAtomType(self, idx: int) -> int:
        """
            Retrieves MMFF atom type for atom with index idx
        
        """
    def GetMMFFBondStretchParams(self, mol: typing.Any, idx1: int, idx2: int) -> typing.Any:
        """
            Retrieves MMFF bond stretch parameters for atoms with indexes idx1, idx2 as a (bondType, kb, r0) tuple, or None if no parameters could be found
        
        """
    def GetMMFFBondTerm(self) -> bool:
        """
            Returns whether the bond term is included in the MMFF equation.
        
        """
    def GetMMFFDielectricConstant(self) -> float:
        """
            Returns the currently configured MMFF dielectric constant.
        
        """
    def GetMMFFDielectricModel(self) -> int:
        """
            Returns the currently configured MMFF dielectric model (1: constant; 2: distance-dependent).
        
        """
    def GetMMFFEleTerm(self) -> bool:
        """
            Returns whether the electrostatic term is included in the MMFF equation.
        
        """
    def GetMMFFFormalCharge(self, idx: int) -> float:
        """
            Retrieves MMFF formal charge for atom with index idx
        
        """
    def GetMMFFOopBendParams(self, mol: typing.Any, idx1: int, idx2: int, idx3: int, idx4: int) -> typing.Any:
        """
            Retrieves MMFF out-of-plane bending force constant for atoms with indexes idx1, idx2, idx3, idx4 as a koop float value
        
        """
    def GetMMFFOopTerm(self) -> bool:
        """
            Returns whether the out-of-plane bend term is included in the MMFF equation.
        
        """
    def GetMMFFPartialCharge(self, idx: int) -> float:
        """
            Retrieves MMFF partial charge for atom with index idx
        
        """
    def GetMMFFStretchBendParams(self, mol: typing.Any, idx1: int, idx2: int, idx3: int) -> typing.Any:
        """
            Retrieves MMFF stretch-bend parameters for atoms with indexes idx1, idx2, idx3 as a (stretchBendType, kbaIJK, kbaKJI) tuple, or None if no parameters could be found
        
        """
    def GetMMFFStretchBendTerm(self) -> bool:
        """
            Returns whether the stretch-bend term is included in the MMFF equation.
        
        """
    def GetMMFFTorsionParams(self, mol: typing.Any, idx1: int, idx2: int, idx3: int, idx4: int) -> typing.Any:
        """
            Retrieves MMFF torsion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a (torsionType, V1, V2, V3) tuple, or None if no parameters could be found
        
        """
    def GetMMFFTorsionTerm(self) -> bool:
        """
            Returns whether the torsional term is included in the MMFF equation.
        
        """
    def GetMMFFVariant(self) -> str:
        """
            Returns the currently configured MMFF variant ("MMFF94" or "MMFF94s").
        
        """
    def GetMMFFVdWParams(self, idx1: int, idx2: int) -> typing.Any:
        """
            Retrieves MMFF van der Waals parameters for atoms with indexes idx1, idx2 as a (R_ij_starUnscaled, epsilonUnscaled, R_ij_star, epsilon) tuple, or None if no parameters could be found
        
        """
    def GetMMFFVdWTerm(self) -> bool:
        """
            Returns whether the Van der Waals term is included in the MMFF equation.
        
        """
    def SetMMFFAngleTerm(self, state: bool = True) -> None:
        """
            Sets the angle term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFBondTerm(self, state: bool = True) -> None:
        """
            Sets the bond term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFDielectricConstant(self, dielConst: float = 1.0) -> None:
        """
            Sets the DielConst MMFF property (defaults to 1.0)
        
        """
    def SetMMFFDielectricModel(self, dielModel: int = 1) -> None:
        """
            Sets the DielModel MMFF property (1: constant; 2: distance-dependent; defaults to constant)
        
        """
    def SetMMFFEleTerm(self, state: bool = True) -> None:
        """
            Sets the electrostatic term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFOopTerm(self, state: bool = True) -> None:
        """
            Sets the out-of-plane bend term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFStretchBendTerm(self, state: bool = True) -> None:
        """
            Sets the stretch-bend term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFTorsionTerm(self, state: bool = True) -> None:
        """
            Sets the torsional term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFVariant(self, mmffVariant: str = 'MMFF94') -> None:
        """
            Sets the MMFF variant to be used ("MMFF94" or "MMFF94s"; defaults to "MMFF94")
        
        """
    def SetMMFFVdWTerm(self, state: bool = True) -> None:
        """
            Sets the Van der Waals term to be included in the MMFF equation (defaults to True)
        
        """
    def SetMMFFVerbosity(self, verbosity: int = 0) -> None:
        """
            Sets the MMFF verbosity (0: none; 1: low; 2: high; defaults to 0)
        
        """
