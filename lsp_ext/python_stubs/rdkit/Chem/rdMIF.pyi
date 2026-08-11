"""
Module containing functions for calculating molecular interaction fields (MIFs)
  NOTE: This functionality is experimental and the API and/or results may change in future releases.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['CalculateDescriptors', 'ConstructGrid', 'Coulomb', 'CoulombDielectric', 'HBond', 'Hydrophilic', 'MMFFVdWaals', 'ReadFromCubeFile', 'UFFVdWaals', 'WriteToCubeFile']
class Coulomb(Boost.Python.instance):
    """
    Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
            vacuum (no dielectric).
    
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
                    vacuum (no dielectric).
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    - electrostatic potential in [kJ mol^-1]
            
        
        """
    @typing.overload
    def __init__(self, mol: Mol, confId: int = -1, probeCharge: float = 1.0, absVal: bool = False, chargeKey: str = '_GasteigerCharge', softcoreParam: float = 0.0, cutoff: float = 1.0) -> None:
        ...
    @typing.overload
    def __init__(self, charges: typing.Any, positions: typing.Any, probeCharge: float = 1.0, absVal: bool = False, softcoreParam: float = 0.0, cutoff: float = 1.0) -> typing.Any:
        ...
class CoulombDielectric(Boost.Python.instance):
    """
    Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
            by taking a distance-dependent dielectric into account.
            Same energy term as used in GRID MIFs.
            References:
            - J. Med. Chem. 1985, 28, 849.
            - J. Comp. Chem. 1983, 4, 187.
    
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
                    by taking a distance-dependent dielectric into account.
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    - electrostatic potential in [kJ mol^-1]
            
        
        """
    @typing.overload
    def __init__(self, mol: Mol, confId: int = -1, probeCharge: float = 1.0, absVal: bool = False, chargeKey: str = '_GasteigerCharge', softcoreParam: float = 0.0, cutoff: float = 1.0, epsilon: float = 80.0, xi: float = 4.0) -> None:
        ...
    @typing.overload
    def __init__(self, charges: typing.Any, positions: typing.Any, probeCharge: float = 1.0, absVal: bool = False, softcoreParam: float = 0.0, cutoff: float = 1.0, epsilon: float = 80.0, xi: float = 4.0) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, pklString: str) -> None:
        ...
class HBond(Boost.Python.instance):
    """
    Class for calculation of hydrogen bonding energy between a probe and a molecule.
    
            Similar to GRID hydrogen bonding descriptors.
            References:
            - J.Med.Chem. 1989, 32, 1083.
            - J.Med.Chem. 1993, 36, 140.
            - J.Med.Chem. 1993, 36, 148.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the hydrogen bonding energy between probe and molecule in
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    hydrogen bonding energy in [kJ mol^-1]
            
        
        """
    def __init__(self, mol: Mol, confId: int = -1, probeAtomType: str = 'OH', fixed: bool = True, cutoff: float = 1.0) -> None:
        ...
class Hydrophilic(Boost.Python.instance):
    """
    Class for calculation of a hydrophilic potential of a molecule at a point.
    
            The interaction energy of hydrogen and oxygen of water is calculated at each point as a 
            hydrogen bond interaction (either OH or O probe). The favored interaction is returned.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the hydrophilic field energy at a point.
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    hydrophilic field energy in [kJ mol^-1]
            
        
        """
    def __init__(self, mol: Mol, confId: int = -1, fixed: bool = True, cutoff: float = 1.0) -> None:
        ...
class MMFFVdWaals(Boost.Python.instance):
    """
    Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the MMFF forcefield.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the van der Waals interaction between molecule and a probe at a gridpoint.
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    - van der Waals potential in [kJ mol^-1]
            
        
        """
    def __init__(self, mol: Mol, confId: int = -1, probeAtomType: int = 6, scaling: bool = False, cutoff: float = 1.0) -> None:
        ...
class UFFVdWaals(Boost.Python.instance):
    """
    Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the UFF forcefield.
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
            Calculates the van der Waals interaction between molecule and a probe at a gridpoint.
            
                    ARGUMENTS:
                    - x, y, z:   coordinates of probe position for energy calculation
                    - threshold: maximal distance until which interactions are calculated
                    RETURNS:
                    - van der Waals potential in [kJ mol^-1]
            
        
        """
    def __init__(self, mol: Mol, confId: int = -1, probeAtomType: str = 'O_3', cutoff: float = 1.0) -> None:
        ...
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: Coulomb, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: CoulombDielectric, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: MMFFVdWaals, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: UFFVdWaals, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: HBond, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
@typing.overload
def CalculateDescriptors(grid: UniformRealValueGrid3D, descriptor: Hydrophilic, threshold: float = -1.0) -> None:
    """
        Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D which get the MIF values
                - descriptor:  Descriptor class which is used to calculate values
        
    
    """
def ConstructGrid(mol: Mol, confId: int = -1, margin: float = 5.0, spacing: float = 0.5) -> UniformRealValueGrid3D:
    """
        Constructs a UniformRealValueGrid3D (3D grid with real values at gridpoints) fitting to a molecule.
        
                ARGUMENTS:
                - mol:     molecule of interest
                - confId:  the ID of the conformer to be used (defaults to -1)
                - margin:  minimum distance of molecule to surface of grid [A] (defaults to 5.0 A)
                - spacing: grid spacing [A] (defaults to 0.5 A)
        
    
    """
def ReadFromCubeFile(filename: str) -> tuple:
    """
        Reads Grid from a file in Gaussian CUBE format.
        
                ARGUMENTS:
                - filename:  filename of file to be read
                RETURNS:
                a tuple where the first element is the grid and
                the second element is the molecule object associated to the grid
                (only atoms and coordinates, no bonds;
                None if no molecule was associated to the grid)
        
    
    """
def WriteToCubeFile(grid: UniformRealValueGrid3D, filename: str, mol: Mol = None, confId: int = -1) -> None:
    """
        Writes Grid to a file in Gaussian CUBE format.
        
                ARGUMENTS:
                - grid:      UniformRealValueGrid3D to be stored
                - filename:  filename of file to be written
                - mol:       associated molecule (defaults to None)
                - confId:    the ID of the conformer to be used (defaults to -1)
        
    
    """
