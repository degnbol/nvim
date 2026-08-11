# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing a C++ implementation of the FMCS algorithm
"""
from __future__ import annotations
import typing
__all__: list[str] = ['AtomCompare', 'BondCompare', 'FindMCS', 'MCSAcceptance', 'MCSAtomCompare', 'MCSAtomCompareParameters', 'MCSBondCompare', 'MCSBondCompareParameters', 'MCSFinalMatchCheck', 'MCSParameters', 'MCSProgress', 'MCSProgressData', 'MCSResult', 'RingCompare']
class AtomCompare(Boost.Python.enum):
    CompareAny: typing.ClassVar[AtomCompare]  # value = rdkit.Chem.rdFMCS.AtomCompare.CompareAny
    CompareAnyHeavyAtom: typing.ClassVar[AtomCompare]  # value = rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom
    CompareElements: typing.ClassVar[AtomCompare]  # value = rdkit.Chem.rdFMCS.AtomCompare.CompareElements
    CompareIsotopes: typing.ClassVar[AtomCompare]  # value = rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'CompareAny': rdkit.Chem.rdFMCS.AtomCompare.CompareAny, 'CompareElements': rdkit.Chem.rdFMCS.AtomCompare.CompareElements, 'CompareIsotopes': rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes, 'CompareAnyHeavyAtom': rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdFMCS.AtomCompare.CompareAny, 1: rdkit.Chem.rdFMCS.AtomCompare.CompareElements, 2: rdkit.Chem.rdFMCS.AtomCompare.CompareIsotopes, 3: rdkit.Chem.rdFMCS.AtomCompare.CompareAnyHeavyAtom}
class BondCompare(Boost.Python.enum):
    CompareAny: typing.ClassVar[BondCompare]  # value = rdkit.Chem.rdFMCS.BondCompare.CompareAny
    CompareOrder: typing.ClassVar[BondCompare]  # value = rdkit.Chem.rdFMCS.BondCompare.CompareOrder
    CompareOrderExact: typing.ClassVar[BondCompare]  # value = rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'CompareAny': rdkit.Chem.rdFMCS.BondCompare.CompareAny, 'CompareOrder': rdkit.Chem.rdFMCS.BondCompare.CompareOrder, 'CompareOrderExact': rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdFMCS.BondCompare.CompareAny, 1: rdkit.Chem.rdFMCS.BondCompare.CompareOrder, 2: rdkit.Chem.rdFMCS.BondCompare.CompareOrderExact}
class MCSAcceptance(Boost.Python.instance):
    """
    Base class. Subclass and override MCSAcceptance.__call__() to define a custom boolean callback function. Returning True will cause the MCS candidate to be accepted, False to be rejected
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self) -> bool:
        """
            override to implement a custom MCS acceptance callback
        
        """
    def __init__(self) -> None:
        ...
class MCSAtomCompare(Boost.Python.instance):
    """
    Base class. Subclass and override MCSAtomCompare.__call__() to define custom atom compare functions, then set MCSParameters.AtomTyper to an instance of the subclass
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def CheckAtomCharge(self, parameters: MCSAtomCompareParameters, mol1: Mol, atom1: int, mol2: Mol, atom2: int) -> bool:
        """
            Return True if both atoms have the same formal charge
        
        """
    def CheckAtomChirality(self, parameters: MCSAtomCompareParameters, mol1: Mol, atom1: int, mol2: Mol, atom2: int) -> bool:
        """
            Return True if both atoms have, or have not, a chiral tag
        
        """
    def CheckAtomRingMatch(self, parameters: MCSAtomCompareParameters, mol1: Mol, atom1: int, mol2: Mol, atom2: int) -> bool:
        """
            Return True if both atoms are, or are not, in a ring
        
        """
    def __call__(self, parameters: MCSAtomCompareParameters, mol1: Mol, atom1: int, mol2: Mol, atom2: int) -> bool:
        """
            override to implement custom atom comparison
        
        """
    def __init__(self) -> None:
        ...
class MCSAtomCompareParameters(Boost.Python.instance):
    """
    Parameters controlling how atom-atom matching is done
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def CompleteRingsOnly(self) -> bool:
        """results cannot include lone ring atoms (default: False)"""
    @CompleteRingsOnly.setter
    def CompleteRingsOnly(self, value: bool) -> None: ...
    @property
    def MatchChiralTag(self) -> bool:
        """include atom chirality in the match (default: False)"""
    @MatchChiralTag.setter
    def MatchChiralTag(self, value: bool) -> None: ...
    @property
    def MatchFormalCharge(self) -> bool:
        """include formal charge in the match (default: False)"""
    @MatchFormalCharge.setter
    def MatchFormalCharge(self, value: bool) -> None: ...
    @property
    def MatchIsotope(self) -> bool:
        """use isotope atom queries in MCSResults (default: False)"""
    @MatchIsotope.setter
    def MatchIsotope(self, value: bool) -> None: ...
    @property
    def MatchValences(self) -> bool:
        """include atom valences in the match (default: False)"""
    @MatchValences.setter
    def MatchValences(self, value: bool) -> None: ...
    @property
    def MaxDistance(self) -> float:
        """Require atoms to be within this many angstroms in 3D (default: -1.0)"""
    @MaxDistance.setter
    def MaxDistance(self, value: float) -> None: ...
    @property
    def RingMatchesRingOnly(self) -> bool:
        """ring atoms are only allowed to match other ring atoms (default: False)"""
    @RingMatchesRingOnly.setter
    def RingMatchesRingOnly(self, value: bool) -> None: ...
class MCSBondCompare(Boost.Python.instance):
    """
    Base class. Subclass and override MCSBondCompare.__call__() to define custom bond compare functions, then set MCSParameters.BondTyper to an instance of the subclass
    """
    __instance_size__: typing.ClassVar[int] = 64
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def CheckBondRingMatch(self, parameters: MCSBondCompareParameters, mol1: Mol, bond1: int, mol2: Mol, bond2: int) -> bool:
        """
            Return True if both bonds are, or are not, part of a ring
        
        """
    def CheckBondStereo(self, parameters: MCSBondCompareParameters, mol1: Mol, bond1: int, mol2: Mol, bond2: int) -> bool:
        """
            Return True if both bonds have, or have not, a stereo descriptor
        
        """
    def __call__(self, parameters: MCSBondCompareParameters, mol1: Mol, bond1: int, mol2: Mol, bond2: int) -> bool:
        """
            override to implement custom bond comparison
        
        """
    def __init__(self) -> None:
        ...
class MCSBondCompareParameters(Boost.Python.instance):
    """
    Parameters controlling how bond-bond matching is done
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def CompleteRingsOnly(self) -> bool:
        """results cannot include partial rings (default: False)"""
    @CompleteRingsOnly.setter
    def CompleteRingsOnly(self, value: bool) -> None: ...
    @property
    def MatchFusedRings(self) -> bool:
        """enforce check on ring fusion, i.e. alpha-methylnaphthalene won't match beta-methylnaphtalene, but decalin will match cyclodecane unless MatchFusedRingsStrict is True (default: False)"""
    @MatchFusedRings.setter
    def MatchFusedRings(self, value: bool) -> None: ...
    @property
    def MatchFusedRingsStrict(self) -> bool:
        """only enforced if MatchFusedRings is True; the ring fusion must be the same in both query and target, i.e. decalin won't match cyclodecane (default: False)"""
    @MatchFusedRingsStrict.setter
    def MatchFusedRingsStrict(self, value: bool) -> None: ...
    @property
    def MatchStereo(self) -> bool:
        """include bond stereo in the comparison (default: False)"""
    @MatchStereo.setter
    def MatchStereo(self, value: bool) -> None: ...
    @property
    def RingMatchesRingOnly(self) -> bool:
        """ring bonds are only allowed to match other ring bonds (default: False)"""
    @RingMatchesRingOnly.setter
    def RingMatchesRingOnly(self, value: bool) -> None: ...
class MCSFinalMatchCheck(Boost.Python.instance):
    """
    Base class. Subclass and override MCSFinalMatchCheck.__call__() to define a custom boolean callback function. Returning True will cause the growing seed to be accepted, False to be rejected
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self) -> bool:
        """
            override to implement a custom seed final match checker callback
        
        """
    def __init__(self) -> None:
        ...
class MCSParameters(Boost.Python.instance):
    """
    Parameters controlling how the MCS is constructed
    """
    __instance_size__: typing.ClassVar[int] = 168
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __setattr__(arg1: typing.Any, arg2: str, arg3: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
    @property
    def AtomCompareParameters(*args, **kwargs):
        """
        parameters for comparing atoms
        """
    @AtomCompareParameters.setter
    def AtomCompareParameters(*args, **kwargs):
        ...
    @property
    def AtomTyper(*args, **kwargs):
        """
        atom typer to be used. Must be one of the members of the rdFMCS.AtomCompare class or an instance of a user-defined subclass of rdFMCS.MCSAtomCompare
        """
    @AtomTyper.setter
    def AtomTyper(*args, **kwargs):
        ...
    @property
    def BondCompareParameters(*args, **kwargs):
        """
        parameters for comparing bonds
        """
    @BondCompareParameters.setter
    def BondCompareParameters(*args, **kwargs):
        ...
    @property
    def BondTyper(*args, **kwargs):
        """
        bond typer to be used. Must be one of the members of the rdFMCS.BondCompare class or an instance of a user-defined subclass of rdFMCS.MCSBondCompare
        """
    @BondTyper.setter
    def BondTyper(*args, **kwargs):
        ...
    @property
    def FinalMatchChecker(*args, **kwargs):
        """
        seed final match checker callback class. Must be a user-defined subclass of rdFMCS.MCSFinalMatchCheck
        """
    @FinalMatchChecker.setter
    def FinalMatchChecker(*args, **kwargs):
        ...
    @property
    def InitialSeed(self) -> str:
        """SMILES string to be used as the seed of the MCS (default: '')"""
    @InitialSeed.setter
    def InitialSeed(self, value: str) -> None: ...
    @property
    def MaximizeBonds(self) -> bool:
        """toggles maximizing the number of bonds (instead of the number of atoms) (default: True)"""
    @MaximizeBonds.setter
    def MaximizeBonds(self, value: bool) -> None: ...
    @property
    def ProgressCallback(*args, **kwargs):
        """
        progress callback class. Must be a user-defined subclass of rdFMCS.Progress
        """
    @ProgressCallback.setter
    def ProgressCallback(*args, **kwargs):
        ...
    @property
    def ShouldAcceptMCS(*args, **kwargs):
        """
        MCS acceptance callback class. Must be a user-defined subclass of rdFMCS.MCSAcceptance
        """
    @ShouldAcceptMCS.setter
    def ShouldAcceptMCS(*args, **kwargs):
        ...
    @property
    def StoreAll(self) -> bool:
        """toggles storage of degenerate MCSs (default: False)"""
    @StoreAll.setter
    def StoreAll(self, value: bool) -> None: ...
    @property
    def Threshold(self) -> float:
        """fraction of the dataset that must contain the MCS (default: 1.0)"""
    @Threshold.setter
    def Threshold(self, value: float) -> None: ...
    @property
    def Timeout(self) -> int:
        """timeout (in seconds) for the calculation (default: 0)"""
    @Timeout.setter
    def Timeout(self, value: int) -> None: ...
    @property
    def Verbose(self) -> bool:
        """toggles verbose mode (default: False)"""
    @Verbose.setter
    def Verbose(self, value: bool) -> None: ...
class MCSProgress(Boost.Python.instance):
    """
    Base class. Subclass and override MCSProgress.__call__() to define a custom callback function
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __call__(self, stat: typing.Any, parameters: typing.Any) -> bool:
        """
            override to implement a custom progress callback
        
        """
    def __init__(self) -> None:
        ...
class MCSProgressData(Boost.Python.instance):
    """
    Information about the MCS progress
    """
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def numAtoms(self) -> int:
        """number of atoms in MCS (default: 0)"""
    @property
    def numBonds(self) -> int:
        """number of bonds in MCS (default: 0)"""
    @property
    def seedProcessed(self) -> int:
        """number of processed seeds (default: 0)"""
class MCSResult(Boost.Python.instance):
    """
    used to return MCS results
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
    @property
    def canceled(*args, **kwargs):
        """
        if True, the MCS calculation did not finish
        """
    @property
    def degenerateSmartsQueryMolDict(*args, **kwargs):
        """
        Dictionary collecting all degenerate (SMARTS, queryMol) pairs (empty if MCSParameters.StoreAll is False)
        """
    @property
    def numAtoms(*args, **kwargs):
        """
        number of atoms in MCS
        """
    @property
    def numBonds(*args, **kwargs):
        """
        number of bonds in MCS
        """
    @property
    def queryMol(*args, **kwargs):
        """
        query molecule for the MCS
        """
    @property
    def smartsString(*args, **kwargs):
        """
        SMARTS string for the MCS
        """
class RingCompare(Boost.Python.enum):
    IgnoreRingFusion: typing.ClassVar[RingCompare]  # value = rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion
    PermissiveRingFusion: typing.ClassVar[RingCompare]  # value = rdkit.Chem.rdFMCS.RingCompare.PermissiveRingFusion
    StrictRingFusion: typing.ClassVar[RingCompare]  # value = rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'IgnoreRingFusion': rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion, 'PermissiveRingFusion': rdkit.Chem.rdFMCS.RingCompare.PermissiveRingFusion, 'StrictRingFusion': rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.rdFMCS.RingCompare.IgnoreRingFusion, 1: rdkit.Chem.rdFMCS.RingCompare.PermissiveRingFusion, 2: rdkit.Chem.rdFMCS.RingCompare.StrictRingFusion}
@typing.overload
def FindMCS(mols: typing.Any, maximizeBonds: bool = True, threshold: float = 1.0, timeout: int = 3600, verbose: bool = False, matchValences: bool = False, ringMatchesRingOnly: bool = False, completeRingsOnly: bool = False, matchChiralTag: bool = False, atomCompare: AtomCompare = ..., bondCompare: BondCompare = ..., ringCompare: RingCompare = ..., seedSmarts: str = '') -> MCSResult:
    """
        Find the MCS for a set of molecules
    
    """
@typing.overload
def FindMCS(mols: typing.Any, parameters: MCSParameters) -> MCSResult:
    """
        Find the MCS for a set of molecules
    
    """
