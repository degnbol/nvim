# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
Module containing functions to compute atomic coordinates in 3D using distance geometry
"""
from __future__ import annotations
import typing
__all__: list[str] = ['BAD_DOUBLE_BOND_STEREO', 'CHECK_CHIRAL_CENTERS', 'CHECK_CHIRAL_CENTERS2', 'CHECK_TETRAHEDRAL_CENTERS', 'CLASH', 'DG', 'ETDG', 'ETDGv2', 'ETKDG', 'ETKDGv2', 'ETKDGv3', 'ETK_MINIMIZATION', 'EXCEEDED_TIMEOUT', 'EmbedFailureCauses', 'EmbedMolecule', 'EmbedMultipleConfs', 'EmbedParameters', 'EmbedParametersToJSON', 'FINAL_CENTER_IN_VOLUME', 'FINAL_CHIRAL_BOUNDS', 'FIRST_MINIMIZATION', 'GetExperimentalTorsions', 'GetMoleculeBoundsMatrix', 'INITIAL_COORDS', 'KDG', 'KTERM_VIOLATION', 'LINEAR_DOUBLE_BOND', 'MINIMIZATION', 'MINIMIZE_FOURTH_DIMENSION', 'OrderedEmbedFailureCauses', 'srETKDGv3']
class EmbedFailureCauses(Boost.Python.enum):
    BAD_DOUBLE_BOND_STEREO: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
    CHECK_CHIRAL_CENTERS: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
    CHECK_CHIRAL_CENTERS2: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS2
    CHECK_TETRAHEDRAL_CENTERS: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
    CLASH: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CLASH
    ETK_MINIMIZATION: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
    EXCEEDED_TIMEOUT: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.EXCEEDED_TIMEOUT
    FINAL_CENTER_IN_VOLUME: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
    FINAL_CHIRAL_BOUNDS: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
    FIRST_MINIMIZATION: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
    INITIAL_COORDS: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
    KTERM_VIOLATION: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.KTERM_VIOLATION
    LINEAR_DOUBLE_BOND: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
    MINIMIZATION: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZATION
    MINIMIZE_FOURTH_DIMENSION: typing.ClassVar[EmbedFailureCauses]  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'INITIAL_COORDS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS, 'FIRST_MINIMIZATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION, 'CHECK_TETRAHEDRAL_CENTERS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS, 'CHECK_CHIRAL_CENTERS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS, 'MINIMIZE_FOURTH_DIMENSION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION, 'ETK_MINIMIZATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION, 'FINAL_CHIRAL_BOUNDS': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS, 'FINAL_CENTER_IN_VOLUME': rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME, 'LINEAR_DOUBLE_BOND': rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND, 'BAD_DOUBLE_BOND_STEREO': rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO, 'CHECK_CHIRAL_CENTERS2': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS2, 'EXCEEDED_TIMEOUT': rdkit.Chem.rdDistGeom.EmbedFailureCauses.EXCEEDED_TIMEOUT, 'MINIMIZATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZATION, 'KTERM_VIOLATION': rdkit.Chem.rdDistGeom.EmbedFailureCauses.KTERM_VIOLATION, 'CLASH': rdkit.Chem.rdDistGeom.EmbedFailureCauses.CLASH}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS, 1: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION, 2: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS, 3: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS, 4: rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION, 5: rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION, 6: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS, 7: rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME, 8: rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND, 9: rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO, 10: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS2, 11: rdkit.Chem.rdDistGeom.EmbedFailureCauses.EXCEEDED_TIMEOUT, 12: rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZATION, 13: rdkit.Chem.rdDistGeom.EmbedFailureCauses.KTERM_VIOLATION, 14: rdkit.Chem.rdDistGeom.EmbedFailureCauses.CLASH}
class EmbedParameters(Boost.Python.instance):
    """
    Parameters controlling embedding
    """
    __instance_size__: typing.ClassVar[int] = 232
    @staticmethod
    def SetCoordMap(arg1: EmbedParameters, self: dict) -> None:
        """
            sets the coordmap to be used
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def GetFailureCounts(self) -> tuple:
        """
            returns the counts of each failure type
        
        """
    def SetBoundsMat(self, boundsMatArg: typing.Any) -> None:
        """
            set the distance-bounds matrix to be used (no triangle smoothing will be done on this) from a Numpy array
        
        """
    def SetCPCI(self, CPCIdict: dict) -> None:
        """
            set the customised pairwise Columb-like interaction to atom pairs.used during structural minimisation stage
        
        """
    def __init__(self) -> None:
        ...
    @property
    def ETversion(self) -> int:
        """version of the experimental torsion-angle preferences (default: 2)"""
    @ETversion.setter
    def ETversion(self, value: int) -> None: ...
    @property
    def basinThresh(self) -> float:
        """set the basin threshold for the DGeom force field. (default: 5.0)"""
    @basinThresh.setter
    def basinThresh(self, value: float) -> None: ...
    @property
    def boundsMatForceScaling(self) -> float:
        """scale the weights of the atom pair distance restraints relative to the other types of restraints (default: 1.0)"""
    @boundsMatForceScaling.setter
    def boundsMatForceScaling(self, value: float) -> None: ...
    @property
    def boxSizeMult(self) -> float:
        """determines the size of the box used for random coordinates (default: 2.0)"""
    @boxSizeMult.setter
    def boxSizeMult(self, value: float) -> None: ...
    @property
    def clearConfs(self) -> bool:
        """clear all existing conformations on the molecule (default: True)"""
    @clearConfs.setter
    def clearConfs(self, value: bool) -> None: ...
    @property
    def embedFragmentsSeparately(self) -> bool:
        """split the molecule into fragments and embed them separately (default: True)"""
    @embedFragmentsSeparately.setter
    def embedFragmentsSeparately(self, value: bool) -> None: ...
    @property
    def enableSequentialRandomSeeds(self) -> bool:
        """handle random number seeds so that conformer generation can be restarted (default: False)"""
    @enableSequentialRandomSeeds.setter
    def enableSequentialRandomSeeds(self, value: bool) -> None: ...
    @property
    def enforceChirality(self) -> bool:
        """enforce correct chirilaty if chiral centers are present (default: True)"""
    @enforceChirality.setter
    def enforceChirality(self, value: bool) -> None: ...
    @property
    def forceTransAmides(self) -> bool:
        """constrain amide bonds to be trans (default: True)"""
    @forceTransAmides.setter
    def forceTransAmides(self, value: bool) -> None: ...
    @property
    def ignoreSmoothingFailures(self) -> bool:
        """try and embed the molecule if if triangle smoothing of the bounds matrix fails (default: False)"""
    @ignoreSmoothingFailures.setter
    def ignoreSmoothingFailures(self, value: bool) -> None: ...
    @property
    def maxIterations(self) -> int:
        """maximum number of embedding attempts to use for a single conformation (default: 0)"""
    @maxIterations.setter
    def maxIterations(self, value: int) -> None: ...
    @property
    def numThreads(self) -> int:
        """number of threads to use when embedding multiple conformations (default: 1)"""
    @numThreads.setter
    def numThreads(self, value: int) -> None: ...
    @property
    def numZeroFail(self) -> int:
        """fail embedding if we have at least this many zero eigenvalues (default: 1)"""
    @numZeroFail.setter
    def numZeroFail(self, value: int) -> None: ...
    @property
    def onlyHeavyAtomsForRMS(self) -> bool:
        """Only consider heavy atoms when doing RMS filtering (default: True)"""
    @onlyHeavyAtomsForRMS.setter
    def onlyHeavyAtomsForRMS(self, value: bool) -> None: ...
    @property
    def optimizerForceTol(self) -> float:
        """the tolerance to be used during the distance-geometry force field minimization (default: 0.001)"""
    @optimizerForceTol.setter
    def optimizerForceTol(self, value: float) -> None: ...
    @property
    def pruneRmsThresh(self) -> float:
        """used to filter multiple conformations: keep only conformations that are at least this far apart from each other (default: -1.0)"""
    @pruneRmsThresh.setter
    def pruneRmsThresh(self, value: float) -> None: ...
    @property
    def randNegEig(self) -> bool:
        """if the embedding yields a negative eigenvalue, pick coordinates that correspond to this component at random (default: True)"""
    @randNegEig.setter
    def randNegEig(self, value: bool) -> None: ...
    @property
    def randomSeed(self) -> int:
        """seed for the random number generator (default: -1)"""
    @randomSeed.setter
    def randomSeed(self, value: int) -> None: ...
    @property
    def symmetrizeConjugatedTerminalGroupsForPruning(self) -> bool:
        """symmetrize terminal conjugated groups for RMSD pruning (default: True)"""
    @symmetrizeConjugatedTerminalGroupsForPruning.setter
    def symmetrizeConjugatedTerminalGroupsForPruning(self, value: bool) -> None: ...
    @property
    def timeout(self) -> int:
        """maximum time in seconds to generate a conformer for a single molecule fragment. If set to 0, no timeout is set (default: 0)"""
    @timeout.setter
    def timeout(self, value: int) -> None: ...
    @property
    def trackFailures(self) -> bool:
        """keep track of which checks during the embedding process fail (default: False)"""
    @trackFailures.setter
    def trackFailures(self, value: bool) -> None: ...
    @property
    def useBasicKnowledge(self) -> bool:
        """impose basic-knowledge constraints such as flat rings (default: False)"""
    @useBasicKnowledge.setter
    def useBasicKnowledge(self, value: bool) -> None: ...
    @property
    def useExpTorsionAnglePrefs(self) -> bool:
        """impose experimental torsion angle preferences (default: False)"""
    @useExpTorsionAnglePrefs.setter
    def useExpTorsionAnglePrefs(self, value: bool) -> None: ...
    @property
    def useLegacyImplementation(self) -> bool:
        """Whether to use the combined minimization approach (default: True)"""
    @useLegacyImplementation.setter
    def useLegacyImplementation(self, value: bool) -> None: ...
    @property
    def useMacrocycle14config(self) -> bool:
        """use the 1-4 distance bounds from ETKDGv3 (default: False)"""
    @useMacrocycle14config.setter
    def useMacrocycle14config(self, value: bool) -> None: ...
    @property
    def useMacrocycleTorsions(self) -> bool:
        """impose macrocycle torsion angle preferences (default: False)"""
    @useMacrocycleTorsions.setter
    def useMacrocycleTorsions(self, value: bool) -> None: ...
    @property
    def useRandomCoords(self) -> bool:
        """start the embedding from random coordinates instead of using eigenvalues of the distance matrix (default: False)"""
    @useRandomCoords.setter
    def useRandomCoords(self, value: bool) -> None: ...
    @property
    def useSmallRingTorsions(self) -> bool:
        """impose small ring torsion angle preferences (default: False)"""
    @useSmallRingTorsions.setter
    def useSmallRingTorsions(self, value: bool) -> None: ...
    @property
    def useSymmetryForPruning(self) -> bool:
        """use molecule symmetry when doing the RMSD pruning. Note that this option automatically also sets onlyHeavyAtomsForRMS to true. (default: True)"""
    @useSymmetryForPruning.setter
    def useSymmetryForPruning(self, value: bool) -> None: ...
    @property
    def verbose(self) -> bool:
        """be verbose about configuration (default: False)"""
    @verbose.setter
    def verbose(self, value: bool) -> None: ...
def DG() -> EmbedParameters:
    """
        Returns an EmbedParameters object for plain distance geometry.
    
    """
def ETDG() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETDG method.
    
    """
def ETDGv2() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETDG method - version 2.
    
    """
def ETKDG() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETKDG method - version 1.
    
    """
def ETKDGv2() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETKDG method - version 2.
    
    """
def ETKDGv3() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETKDG method - version 3 (macrocycles).
    
    """
@typing.overload
def EmbedMolecule(mol: Mol, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, ETversion: int = 2, useMacrocycle14config: bool = True) -> int:
    """
        Use distance geometry to obtain initial 
         coordinates for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - maxAttempts : maximum number of embedding attempts to use for a single conformation 
            - randomSeed : provide a seed for the random number generator 
                           so that the same coordinates can be obtained 
                           for a molecule on multiple runs. If -1, the 
                           RNG will not be seeded. 
            - clearConfs : clear all existing conformations on the molecule
            - useRandomCoords : Start the embedding from random coordinates instead of
                                using eigenvalues of the distance matrix.
            - boxSizeMult :  Determines the size of the box that is used for
                             random coordinates. If this is a positive number, the 
                             side length will equal the largest element of the distance
                             matrix times boxSizeMult. If this is a negative number,
                             the side length will equal -boxSizeMult (i.e. independent
                             of the elements of the distance matrix).
            - randNegEig : If the embedding yields a negative eigenvalue, 
                           pick coordinates that correspond 
                           to this component at random 
            - numZeroFail : fail embedding if we have at least this many zero eigenvalues 
            - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                         require some atoms to have fixed coordinates in the resulting 
                         conformation.
            - forceTol : tolerance to be used during the force-field minimization with 
                         the distance geometry force field.
            - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                         of the bounds matrix fails.
            - enforceChirality : enforce the correct chirality if chiral centers are present.
            - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
            - useBasicKnowledge : impose basic knowledge such as flat rings
            - printExpTorsionAngles : print the output from the experimental torsion angles
            - useMacrocycleTorsions : use additional torsion profiles for macrocycles
            - ETversion : version of the standard torsion definitions to use. NOTE for both
                          ETKDGv2 and ETKDGv3 this should be 2 since ETKDGv3 uses the ETKDGv2
                          definitions for standard torsions
            - useMacrocycle14config : use the 1-4 distance bounds from ETKDGv3
        
         RETURNS:
        
            ID of the new conformation added to the molecule or -1 if the embedding fails.
        
        
    
    """
@typing.overload
def EmbedMolecule(mol: Mol, params: EmbedParameters) -> int:
    """
        Use distance geometry to obtain intial 
         coordinates for a molecule
        
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - params : an EmbedParameters object 
        
         RETURNS:
        
            ID of the new conformation added to the molecule or -1 if the embedding fails. 
        
        
    
    """
@typing.overload
def EmbedMultipleConfs(mol: Mol, numConfs: int = 10, maxAttempts: int = 0, randomSeed: int = -1, clearConfs: bool = True, useRandomCoords: bool = False, boxSizeMult: float = 2.0, randNegEig: bool = True, numZeroFail: int = 1, pruneRmsThresh: float = -1.0, coordMap: dict = {}, forceTol: float = 0.001, ignoreSmoothingFailures: bool = False, enforceChirality: bool = True, numThreads: int = 1, useExpTorsionAnglePrefs: bool = True, useBasicKnowledge: bool = True, printExpTorsionAngles: bool = False, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, ETversion: int = 2, useMacrocycle14config: bool = True) -> typing.Sequence[int]:
    """
        Use distance geometry to obtain multiple sets of 
         coordinates for a molecule
         
         ARGUMENTS:
        
          - mol : the molecule of interest
          - numConfs : the number of conformers to generate 
          - maxAttempts : maximum number of embedding attempts to use for a single conformation 
          - randomSeed : provide a seed for the random number generator 
                         so that the same coordinates can be obtained 
                         for a molecule on multiple runs. If -1, the 
                         RNG will not be seeded. 
          - clearConfs : clear all existing conformations on the molecule
          - useRandomCoords : Start the embedding from random coordinates instead of
                              using eigenvalues of the distance matrix.
          - boxSizeMult    Determines the size of the box that is used for
                           random coordinates. If this is a positive number, the 
                           side length will equal the largest element of the distance
                           matrix times boxSizeMult. If this is a negative number,
                           the side length will equal -boxSizeMult (i.e. independent
                           of the elements of the distance matrix).
          - randNegEig : If the embedding yields a negative eigenvalue, 
                         pick coordinates that correspond 
                         to this component at random 
          - numZeroFail : fail embedding if we have at least this many zero eigenvalues 
          - pruneRmsThresh : Retain only the conformations out of 'numConfs' 
                            after embedding that are at least 
                            this far apart from each other. 
                            RMSD is computed on the heavy atoms. 
                            Pruning is greedy; i.e. the first embedded conformation
                            is retained and from then on only those that are at
                            least pruneRmsThresh away from all retained conformations
                            are kept. The pruning is done after embedding and 
                            bounds violation minimization. No pruning by default.
          - coordMap : a dictionary mapping atom IDs->coordinates. Use this to 
                       require some atoms to have fixed coordinates in the resulting 
                       conformation.
          - forceTol : tolerance to be used during the force-field minimization with 
                       the distance geometry force field.
          - ignoreSmoothingFailures : try to embed the molecule even if triangle smoothing
                       of the bounds matrix fails.
          - enforceChirality : enforce the correct chirality if chiral centers are present.
          - numThreads : number of threads to use while embedding. This only has an effect if the RDKit
                       was built with multi-thread support.
                      If set to zero, the max supported by the system will be used.
          - useExpTorsionAnglePrefs : impose experimental torsion angle preferences
          - useBasicKnowledge : impose basic knowledge such as flat rings
          - printExpTorsionAngles : print the output from the experimental torsion angles
         RETURNS:
        
            Iterator which yields new conformation IDs 
        
        
    
    """
@typing.overload
def EmbedMultipleConfs(mol: Mol, numConfs: int, params: EmbedParameters) -> typing.Sequence[int]:
    """
        Use distance geometry to obtain multiple sets of 
         coordinates for a molecule
         
         ARGUMENTS:
        
          - mol : the molecule of interest
          - numConfs : the number of conformers to generate 
          - params : an EmbedParameters object 
         RETURNS:
        
            Iterator which yields new conformation IDs 
        
        
    
    """
def EmbedParametersToJSON(embedParameters: EmbedParameters) -> str:
    """
        Returns json string containing embedParameters attributes
          
          ARGUMENTS:
        
            - embedParameters : the Params object you want serialized
          RETURNS:
        
            The Params object as json string
          
        
    
    """
@typing.overload
def GetExperimentalTorsions(mol: Mol, useExpTorsionAnglePrefs: bool = True, useSmallRingTorsions: bool = False, useMacrocycleTorsions: bool = True, useBasicKnowledge: bool = True, ETversion: int = 2, printExpTorsionAngles: bool = False) -> tuple:
    """
        returns information about the bonds corresponding to experimental torsions
    
    """
@typing.overload
def GetExperimentalTorsions(mol: Mol, embedParams: EmbedParameters) -> tuple:
    """
        returns information about the bonds corresponding to experimental torsions
    
    """
@typing.overload
def GetMoleculeBoundsMatrix(mol: Mol, set15bounds: bool = True, scaleVDW: bool = False, doTriangleSmoothing: bool = True, useMacrocycle14config: bool = False, forceTransAmides: bool = True, set14bounds: bool = True, set13bounds: bool = True) -> typing.Any:
    """
        Returns the distance bounds matrix for a molecule
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - set15bounds : set bounds for 1-5 atom distances based on 
                            topology (otherwise stop at 1-4s)
            - scaleVDW : ignored 
            - doTriangleSmoothing : run triangle smoothing on the bounds 
                         matrix before returning it 
            - useMacrocycle14config : use 1-4 distance bound heuristics for macrocycles
            - forceTransAmides : constrain amide bonds to be trans
            - set14bounds : set bounds for 1-4 atom distances based on 
                            topology
            - set13bounds : set bounds for 1-3 atom distances based on 
                            topology
         RETURNS:
        
            the bounds matrix as a Numeric array with lower bounds in 
            the lower triangle and upper bounds in the upper triangle
        
        
    
    """
@typing.overload
def GetMoleculeBoundsMatrix(mol: Mol, embedParams: EmbedParameters, doTriangleSmoothing: bool = True, scaleVDW: bool = False, set15bounds: bool = True, set14bounds: bool = True, set13bounds: bool = True) -> typing.Any:
    """
        Returns the distance bounds matrix for a molecule
         
         ARGUMENTS:
        
            - mol : the molecule of interest
            - set15bounds : set bounds for 1-5 atom distances based on 
                            topology (otherwise stop at 1-4s)
            - scaleVDW : ignored 
            - doTriangleSmoothing : run triangle smoothing on the bounds 
                         matrix before returning it 
            - useMacrocycle14config : use 1-4 distance bound heuristics for macrocycles
            - forceTransAmides : constrain amide bonds to be trans
            - set14bounds : set bounds for 1-4 atom distances based on 
                            topology
            - set13bounds : set bounds for 1-3 atom distances based on 
                            topology
         RETURNS:
        
            the bounds matrix as a Numeric array with lower bounds in 
            the lower triangle and upper bounds in the upper triangle
        
        
    
    """
def KDG() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the KDG method.
    
    """
def OrderedEmbedFailureCauses(legacyImplementation: bool = True) -> list:
    """
    """
def srETKDGv3() -> EmbedParameters:
    """
        Returns an EmbedParameters object for the ETKDG method - version 3 (small rings).
    
    """
BAD_DOUBLE_BOND_STEREO: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.BAD_DOUBLE_BOND_STEREO
CHECK_CHIRAL_CENTERS: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS
CHECK_CHIRAL_CENTERS2: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_CHIRAL_CENTERS2
CHECK_TETRAHEDRAL_CENTERS: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CHECK_TETRAHEDRAL_CENTERS
CLASH: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.CLASH
ETK_MINIMIZATION: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.ETK_MINIMIZATION
EXCEEDED_TIMEOUT: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.EXCEEDED_TIMEOUT
FINAL_CENTER_IN_VOLUME: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CENTER_IN_VOLUME
FINAL_CHIRAL_BOUNDS: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FINAL_CHIRAL_BOUNDS
FIRST_MINIMIZATION: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.FIRST_MINIMIZATION
INITIAL_COORDS: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.INITIAL_COORDS
KTERM_VIOLATION: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.KTERM_VIOLATION
LINEAR_DOUBLE_BOND: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.LINEAR_DOUBLE_BOND
MINIMIZATION: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZATION
MINIMIZE_FOURTH_DIMENSION: EmbedFailureCauses  # value = rdkit.Chem.rdDistGeom.EmbedFailureCauses.MINIMIZE_FOURTH_DIMENSION
