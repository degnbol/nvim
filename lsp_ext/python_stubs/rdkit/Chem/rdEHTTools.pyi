"""
Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next.
"""
from __future__ import annotations
import typing
__all__: list[str] = ['EHTResults', 'RunMol']
class EHTResults(Boost.Python.instance):
    """
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
    def GetAtomicCharges(self) -> typing.Any:
        """
            returns the calculated atomic charges
        
        """
    def GetHamiltonian(self) -> typing.Any:
        """
            returns the symmetric Hamiltonian matrix
        
        """
    def GetOrbitalEnergies(self) -> typing.Any:
        """
            returns the energies of the molecular orbitals as a vector
        
        """
    def GetOverlapMatrix(self) -> typing.Any:
        """
            returns the symmetric overlap matrix
        
        """
    def GetReducedChargeMatrix(self) -> typing.Any:
        """
            returns the reduced charge matrix
        
        """
    def GetReducedOverlapPopulationMatrix(self) -> typing.Any:
        """
            returns the reduced overlap population matrix
        
        """
    @property
    def fermiEnergy(*args, **kwargs):
        ...
    @property
    def numElectrons(*args, **kwargs):
        ...
    @property
    def numOrbitals(*args, **kwargs):
        ...
    @property
    def totalEnergy(*args, **kwargs):
        ...
def RunMol(mol: Mol, confId: int = -1, keepOverlapAndHamiltonianMatrices: bool = False) -> tuple:
    """
        Runs an extended Hueckel calculation for a molecule.
        The molecule should have at least one conformation
        
        ARGUMENTS:
           - mol: molecule to use
           - confId: (optional) conformation to use
           - keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap 
             and hamiltonian matrices in the EHTResults object
        
        RETURNS: a 2-tuple:
           - a boolean indicating whether or not the calculation succeeded
           - an EHTResults object with the results
        
    
    """
