"""
Module containing tools for normalizing molecules defined by SMARTS patterns
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['AllowedAtomsValidation', 'CHARGE_CORRECTIONS', 'CanonicalTautomer', 'ChargeCorrection', 'ChargeParent', 'ChargeParentInPlace', 'Cleanup', 'CleanupInPlace', 'CleanupParameters', 'DisallowedAtomsValidation', 'DisallowedRadicalValidation', 'DisconnectOrganometallics', 'DisconnectOrganometallicsInPlace', 'FeaturesValidation', 'FragmentParent', 'FragmentParentInPlace', 'FragmentRemover', 'FragmentRemoverFromData', 'FragmentValidation', 'GetDefaultTautomerScoreSubstructs', 'GetV1TautomerEnumerator', 'Is2DValidation', 'IsotopeParent', 'IsotopeParentInPlace', 'IsotopeValidation', 'LargestFragmentChooser', 'Layout2DValidation', 'MOL_SPTR_VECT', 'MetalDisconnector', 'MetalDisconnectorOptions', 'MolVSValidation', 'NeutralValidation', 'NoAtomValidation', 'Normalize', 'NormalizeInPlace', 'Normalizer', 'NormalizerFromData', 'NormalizerFromParams', 'Pipeline', 'PipelineLog', 'PipelineLogEntry', 'PipelineOptions', 'PipelineResult', 'PipelineStage', 'PipelineStatus', 'RDKitValidation', 'Reionize', 'ReionizeInPlace', 'Reionizer', 'ReionizerFromData', 'RemoveFragments', 'RemoveFragmentsInPlace', 'ScoreHeteroHs', 'ScoreRings', 'ScoreSubstructs', 'SmilesTautomerMap', 'StandardizeSmiles', 'StereoParent', 'StereoParentInPlace', 'StereoValidation', 'SubstructTerm', 'SubstructTermVector', 'SuperParent', 'SuperParentInPlace', 'Tautomer', 'TautomerEnumerator', 'TautomerEnumeratorCallback', 'TautomerEnumeratorResult', 'TautomerEnumeratorStatus', 'TautomerParent', 'TautomerParentInPlace', 'Uncharger', 'UpdateParamsFromJSON', 'ValidateSmiles', 'ValidationMethod', 'map_indexing_suite_SmilesTautomerMap_entry']
class AllowedAtomsValidation(ValidationMethod):
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, arg1: typing.Any) -> typing.Any:
        ...
class ChargeCorrection(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 80
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, name: str, smarts: str, charge: int) -> None:
        ...
    @property
    def Charge(*args, **kwargs):
        ...
    @Charge.setter
    def Charge(*args, **kwargs):
        ...
    @property
    def Name(*args, **kwargs):
        ...
    @Name.setter
    def Name(*args, **kwargs):
        ...
    @property
    def Smarts(*args, **kwargs):
        ...
    @Smarts.setter
    def Smarts(*args, **kwargs):
        ...
class CleanupParameters(Boost.Python.instance):
    """
    Parameters controlling molecular standardization
    """
    __instance_size__: typing.ClassVar[int] = 272
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
    def acidbaseFile(self) -> str:
        """file containing the acid and base definitions (default: '')"""
    @acidbaseFile.setter
    def acidbaseFile(self, value: str) -> None: ...
    @property
    def doCanonical(self) -> bool:
        """apply atom-order dependent normalizations (like uncharging) in a canonical order (default: True)"""
    @doCanonical.setter
    def doCanonical(self, value: bool) -> None: ...
    @property
    def fragmentFile(self) -> str:
        """file containing the acid and base definitions (default: '')"""
    @fragmentFile.setter
    def fragmentFile(self, value: str) -> None: ...
    @property
    def largestFragmentChooserCountHeavyAtomsOnly(self) -> bool:
        """whether LargestFragmentChooser should only count heavy atoms (defaults to False) (default: False)"""
    @largestFragmentChooserCountHeavyAtomsOnly.setter
    def largestFragmentChooserCountHeavyAtomsOnly(self, value: bool) -> None: ...
    @property
    def largestFragmentChooserUseAtomCount(self) -> bool:
        """Whether LargestFragmentChooser should use atom count as main criterion before MW (defaults to True) (default: True)"""
    @largestFragmentChooserUseAtomCount.setter
    def largestFragmentChooserUseAtomCount(self, value: bool) -> None: ...
    @property
    def maxRestarts(self) -> int:
        """maximum number of restarts (default: 200)"""
    @maxRestarts.setter
    def maxRestarts(self, value: int) -> None: ...
    @property
    def maxTautomers(self) -> int:
        """maximum number of tautomers to generate (defaults to 1000) (default: 1000)"""
    @maxTautomers.setter
    def maxTautomers(self, value: int) -> None: ...
    @property
    def maxTransforms(self) -> int:
        """maximum number of transforms to apply during tautomer enumeration (defaults to 1000) (default: 1000)"""
    @maxTransforms.setter
    def maxTransforms(self, value: int) -> None: ...
    @property
    def normalizationsFile(self) -> str:
        """file containing the normalization transformations (default: '')"""
    @normalizationsFile.setter
    def normalizationsFile(self, value: str) -> None: ...
    @property
    def preferOrganic(self) -> bool:
        """prefer organic fragments to inorganic ones when deciding what to keep (default: False)"""
    @preferOrganic.setter
    def preferOrganic(self, value: bool) -> None: ...
    @property
    def tautomerReassignStereo(self) -> bool:
        """call AssignStereochemistry on all generated tautomers (defaults to True) (default: True)"""
    @tautomerReassignStereo.setter
    def tautomerReassignStereo(self, value: bool) -> None: ...
    @property
    def tautomerRemoveBondStereo(self) -> bool:
        """remove stereochemistry from double bonds involved in tautomerism (defaults to True) (default: True)"""
    @tautomerRemoveBondStereo.setter
    def tautomerRemoveBondStereo(self, value: bool) -> None: ...
    @property
    def tautomerRemoveIsotopicHs(self) -> bool:
        """remove isotopic Hs from centers involved in tautomerism (defaults to True) (default: True)"""
    @tautomerRemoveIsotopicHs.setter
    def tautomerRemoveIsotopicHs(self, value: bool) -> None: ...
    @property
    def tautomerRemoveSp3Stereo(self) -> bool:
        """remove stereochemistry from sp3 centers involved in tautomerism (defaults to True) (default: True)"""
    @tautomerRemoveSp3Stereo.setter
    def tautomerRemoveSp3Stereo(self, value: bool) -> None: ...
    @property
    def tautomerTransformsFile(self) -> str:
        """file containing the tautomer transformations (default: '')"""
    @tautomerTransformsFile.setter
    def tautomerTransformsFile(self, value: str) -> None: ...
class DisallowedAtomsValidation(ValidationMethod):
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, arg1: typing.Any) -> typing.Any:
        ...
class DisallowedRadicalValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class FeaturesValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, allowEnhancedStereo: bool = False, allowAromaticBondType: bool = False, allowDativeBondType: bool = False, allowQueries: bool = False, allowDummmies: bool = False, allowAtomAliases: bool = False) -> None:
        ...
    @property
    def allowAromaticBondType(self) -> bool:
        """default: False"""
    @allowAromaticBondType.setter
    def allowAromaticBondType(self, value: bool) -> None: ...
    @property
    def allowAtomAliases(self) -> bool:
        """default: False"""
    @allowAtomAliases.setter
    def allowAtomAliases(self, value: bool) -> None: ...
    @property
    def allowDativeBondType(self) -> bool:
        """default: False"""
    @allowDativeBondType.setter
    def allowDativeBondType(self, value: bool) -> None: ...
    @property
    def allowDummies(self) -> bool:
        """default: False"""
    @allowDummies.setter
    def allowDummies(self, value: bool) -> None: ...
    @property
    def allowEnhancedStereo(self) -> bool:
        """default: False"""
    @allowEnhancedStereo.setter
    def allowEnhancedStereo(self, value: bool) -> None: ...
    @property
    def allowQueries(self) -> bool:
        """default: False"""
    @allowQueries.setter
    def allowQueries(self, value: bool) -> None: ...
class FragmentRemover(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, fragmentFilename: str = '', leave_last: bool = True, skip_if_all_match: bool = False) -> None:
        ...
    def remove(self, mol: Mol) -> rdkit.Chem.Mol:
        """
        """
    def removeInPlace(self, mol: Mol) -> None:
        """
            modifies the molecule in place
        
        """
class FragmentValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class Is2DValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, threshold: float = 0.001) -> None:
        ...
    @property
    def threshold(self) -> float:
        """default: 0.001"""
    @threshold.setter
    def threshold(self, value: float) -> None: ...
class IsotopeValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, strict: bool = False) -> None:
        ...
    @property
    def strict(self) -> bool:
        """default: False"""
    @strict.setter
    def strict(self, value: bool) -> None: ...
class LargestFragmentChooser(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self, preferOrganic: bool = False) -> None:
        ...
    @typing.overload
    def __init__(self, params: CleanupParameters) -> None:
        ...
    def choose(self, mol: Mol) -> rdkit.Chem.Mol:
        """
        """
    def chooseInPlace(self, mol: Mol) -> None:
        """
        """
class Layout2DValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 64
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, clashLimit: float = 0.15, bondLengthLimit: float = 25.0, allowLongBondsInRings: bool = True, allowAtomBondClashExemption: bool = True, minMedianBondLength: float = False) -> None:
        ...
    @property
    def allowAtomBondClashExemption(self) -> bool:
        """default: True"""
    @allowAtomBondClashExemption.setter
    def allowAtomBondClashExemption(self, value: bool) -> None: ...
    @property
    def allowLongBondsInRings(self) -> bool:
        """default: True"""
    @allowLongBondsInRings.setter
    def allowLongBondsInRings(self, value: bool) -> None: ...
    @property
    def bondLengthLimit(self) -> float:
        """default: 25.0"""
    @bondLengthLimit.setter
    def bondLengthLimit(self, value: float) -> None: ...
    @property
    def clashLimit(self) -> float:
        """default: 0.15"""
    @clashLimit.setter
    def clashLimit(self, value: float) -> None: ...
    @property
    def minMedianBondLength(self) -> float:
        """default: 0.0"""
    @minMedianBondLength.setter
    def minMedianBondLength(self, value: float) -> None: ...
class MOL_SPTR_VECT(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class MetalDisconnector(Boost.Python.instance):
    """
    a class to disconnect metals that are defined as covalently bonded to non-metals
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def Disconnect(self, mol: Mol) -> rdkit.Chem.Mol:
        """
            performs the disconnection
        
        """
    def DisconnectInPlace(self, mol: Mol) -> None:
        """
            performs the disconnection, modifies the input molecule
        
        """
    def SetMetalNof(self, mol: Mol) -> None:
        """
            Set the query molecule defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine.
        
        """
    def SetMetalNon(self, mol: Mol) -> None:
        """
            Set the query molecule defining the metals to disconnect from other inorganic elements.
        
        """
    def __init__(self, options: typing.Any = None) -> None:
        ...
    @property
    def MetalNof(self) -> str:
        """SMARTS defining the metals to disconnect if attached to Nitrogen, Oxygen or Fluorine (default: '[#3,#11,#19,#37,#55,#87,#4,#12,#20,#38,#56,#88,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#13,#31,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#49,#50,#72,#73,#74,#75,#76,#77,#78,#79,#80,#81,#82,#83]~[#7,#8,F]')"""
    @property
    def MetalNon(self) -> str:
        """SMARTS defining the metals to disconnect other inorganic elements (default: '[#13,#21,#22,#23,#24,#25,#26,#27,#28,#29,#30,#39,#40,#41,#42,#43,#44,#45,#46,#47,#48,#72,#73,#74,#75,#76,#77,#78,#79]~[B,C,#14,P,#33,#51,S,#34,#52,Cl,Br,I,#85]')"""
class MetalDisconnectorOptions(Boost.Python.instance):
    """
    Metal Disconnector Options
    """
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def adjustCharges(self) -> bool:
        """Whether to adjust charges on ligand atoms.  Default true. (default: True)"""
    @adjustCharges.setter
    def adjustCharges(self, value: bool) -> None: ...
    @property
    def removeHapticDummies(self) -> bool:
        """Whether to remove the dummy atoms representing haptic bonds.  Such dummies are bonded to the metal with a bond that has the MolFileBondEndPts prop set.  Default false. (default: False)"""
    @removeHapticDummies.setter
    def removeHapticDummies(self, value: bool) -> None: ...
    @property
    def splitAromaticC(self) -> bool:
        """Whether to split metal-aromatic C bonds.  Default false. (default: False)"""
    @splitAromaticC.setter
    def splitAromaticC(self, value: bool) -> None: ...
    @property
    def splitGrignards(self) -> bool:
        """Whether to split Grignard-type complexes. Default false. (default: False)"""
    @splitGrignards.setter
    def splitGrignards(self, value: bool) -> None: ...
class MolVSValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg1: typing.Any) -> typing.Any:
        ...
class NeutralValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class NoAtomValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class Normalizer(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, normalizeFilename: str, maxRestarts: int) -> None:
        ...
    def normalize(self, mol: Mol) -> rdkit.Chem.Mol:
        """
        """
    def normalizeInPlace(self, mol: Mol) -> None:
        """
            modifies the input molecule
        
        """
class Pipeline(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 240
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def run(arg1: Pipeline, arg2: str) -> PipelineResult:
        """
        """
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg1: PipelineOptions) -> None:
        ...
class PipelineLog(Boost.Python.instance):
    @staticmethod
    def __init__(*args, **kwargs):
        """
        Raises an exception
        This class cannot be instantiated from Python
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class PipelineLogEntry(Boost.Python.instance):
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
    def detail(*args, **kwargs):
        ...
    @property
    def status(*args, **kwargs):
        ...
class PipelineOptions(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 168
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    @property
    def allowAromaticBondType(self) -> bool:
        """default: False"""
    @allowAromaticBondType.setter
    def allowAromaticBondType(self, value: bool) -> None: ...
    @property
    def allowAtomBondClashExemption(self) -> bool:
        """default: True"""
    @allowAtomBondClashExemption.setter
    def allowAtomBondClashExemption(self, value: bool) -> None: ...
    @property
    def allowDativeBondType(self) -> bool:
        """default: False"""
    @allowDativeBondType.setter
    def allowDativeBondType(self, value: bool) -> None: ...
    @property
    def allowEmptyMolecules(self) -> bool:
        """default: False"""
    @allowEmptyMolecules.setter
    def allowEmptyMolecules(self, value: bool) -> None: ...
    @property
    def allowEnhancedStereo(self) -> bool:
        """default: False"""
    @allowEnhancedStereo.setter
    def allowEnhancedStereo(self, value: bool) -> None: ...
    @property
    def allowLongBondsInRings(self) -> bool:
        """default: True"""
    @allowLongBondsInRings.setter
    def allowLongBondsInRings(self, value: bool) -> None: ...
    @property
    def atomClashLimit(self) -> float:
        """default: 0.03"""
    @atomClashLimit.setter
    def atomClashLimit(self, value: float) -> None: ...
    @property
    def bondLengthLimit(self) -> float:
        """default: 100.0"""
    @bondLengthLimit.setter
    def bondLengthLimit(self, value: float) -> None: ...
    @property
    def is2DZeroThreshold(self) -> float:
        """default: 0.001"""
    @is2DZeroThreshold.setter
    def is2DZeroThreshold(self, value: float) -> None: ...
    @property
    def metalNof(self) -> str:
        """default: '[Li,Na,K,Rb,Cs,Fr]~[#7,#8,F]'"""
    @metalNof.setter
    def metalNof(self, value: str) -> None: ...
    @property
    def metalNon(self) -> str:
        """default: ''"""
    @metalNon.setter
    def metalNon(self, value: str) -> None: ...
    @property
    def minMedianBondLength(self) -> float:
        """default: 0.001"""
    @minMedianBondLength.setter
    def minMedianBondLength(self, value: float) -> None: ...
    @property
    def normalizerData(self) -> str:
        """default: '// Name\tSMIRKS\nNitro to N+(O-)=O\t[N,P,As,Sb;X3:1](=[O,S,Se,Te:2])=[O,S,Se,Te:3]>>[*+1:1]([*-1:2])=[*:3]\nSulfone to S(=O)(=O)\t[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])\nPyridine oxide to n+O-\t[nH0+0:1]=[OH0+0:2]>>[n+:1][O-:2]\nAzide to N=N+=N-\t[*:1][N:2]=[N:3]#[N:4]>>[*:1][N:2]=[N+:3]=[N-:4]\nDiazo/azo to =N+=N-\t[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]\n[SH](=O)(=O) to S(=O)O\t[c,C,N,O,F,Cl,Br,I:1][SH+0:2](=[O:3])=[O:4]>>[*:1][*:2]([*:3])=[*:4]\nPhosphate to P(O-)=O\t[O-:1][P+;D4:2][O,S,Se,Te;-1:3]>>[O+0:1]=[P+0;D5:2][*-1:3]\nGeneralized phosphate to P(X-)=Y\t[S,Se,Te;-1:1][P+;D4:2][S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]\nC/S+N to C/S=N+\t[C,S&!$([S+]-[O-]);X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\nP+N to P=N+\t[P;X4+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]\nRecombine 1,3-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]\nRecombine 1,3-separated charges\t[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]\nRecombine 1,3-separated charges\t[N,O,P,S;-1:1]-[a+0:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]\nRecombine 1,5-separated charges\t[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\nRecombine 1,5-separated charges\t[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]\nRecombine 1,5-separated charges\t[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]\nNormalize 1,3 conjugated cation\t[N;+0!H0:1]@-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\nNormalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\nNormalize 1,3 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]\nNormalize 1,3 conjugated cation\t[n;+0!H0:1]:[c:2]=[N!$(*~[N,O,P,S;-1]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]\nNormalize 1,5 conjugated cation\t[N;+0!H0:1]@-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\nNormalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]@-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\nNormalize 1,5 conjugated cation\t[N,O!$(*N);+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]\nNormalize 1,5 conjugated cation\t[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*~[N,O,P,S;-1]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]\nCharge normalization\t[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]\nCharge recombination\t[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]\n'"""
    @normalizerData.setter
    def normalizerData(self, value: str) -> None: ...
    @property
    def normalizerMaxRestarts(self) -> int:
        """default: 200"""
    @normalizerMaxRestarts.setter
    def normalizerMaxRestarts(self, value: int) -> None: ...
    @property
    def outputV2000(self) -> bool:
        """default: False"""
    @outputV2000.setter
    def outputV2000(self, value: bool) -> None: ...
    @property
    def reportAllFailures(self) -> bool:
        """default: True"""
    @reportAllFailures.setter
    def reportAllFailures(self, value: bool) -> None: ...
    @property
    def scaledMedianBondLength(self) -> float:
        """default: 1.0"""
    @scaledMedianBondLength.setter
    def scaledMedianBondLength(self, value: float) -> None: ...
    @property
    def strictParsing(self) -> bool:
        """default: False"""
    @strictParsing.setter
    def strictParsing(self, value: bool) -> None: ...
class PipelineResult(Boost.Python.instance):
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
    def inputMolData(*args, **kwargs):
        ...
    @property
    def log(*args, **kwargs):
        ...
    @property
    def outputMolData(*args, **kwargs):
        ...
    @property
    def parentMolData(*args, **kwargs):
        ...
    @property
    def stage(*args, **kwargs):
        ...
    @property
    def status(*args, **kwargs):
        ...
class PipelineStage(Boost.Python.enum):
    COMPLETED: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.COMPLETED
    PARSING_INPUT: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PARSING_INPUT
    PREPARE_FOR_STANDARDIZATION: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_STANDARDIZATION
    PREPARE_FOR_VALIDATION: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_VALIDATION
    SERIALIZING_OUTPUT: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.SERIALIZING_OUTPUT
    STANDARDIZATION: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.STANDARDIZATION
    VALIDATION: typing.ClassVar[PipelineStage]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.VALIDATION
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'PARSING_INPUT': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PARSING_INPUT, 'PREPARE_FOR_VALIDATION': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_VALIDATION, 'VALIDATION': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.VALIDATION, 'PREPARE_FOR_STANDARDIZATION': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_STANDARDIZATION, 'STANDARDIZATION': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.STANDARDIZATION, 'SERIALIZING_OUTPUT': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.SERIALIZING_OUTPUT, 'COMPLETED': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.COMPLETED}
    values: typing.ClassVar[dict]  # value = {1: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PARSING_INPUT, 2: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_VALIDATION, 3: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.VALIDATION, 4: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.PREPARE_FOR_STANDARDIZATION, 5: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.STANDARDIZATION, 9: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.SERIALIZING_OUTPUT, 10: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStage.COMPLETED}
class PipelineStatus(Boost.Python.enum):
    BASIC_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.BASIC_VALIDATION_ERROR
    CHARGE_STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.CHARGE_STANDARDIZATION_ERROR
    FEATURES_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FEATURES_VALIDATION_ERROR
    FRAGMENTS_REMOVED: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENTS_REMOVED
    FRAGMENT_STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENT_STANDARDIZATION_ERROR
    INPUT_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.INPUT_ERROR
    IS2D_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.IS2D_VALIDATION_ERROR
    LAYOUT2D_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.LAYOUT2D_VALIDATION_ERROR
    METALS_DISCONNECTED: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METALS_DISCONNECTED
    METAL_STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METAL_STANDARDIZATION_ERROR
    NORMALIZATION_APPLIED: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZATION_APPLIED
    NORMALIZER_STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZER_STANDARDIZATION_ERROR
    NO_EVENT: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NO_EVENT
    OUTPUT_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.OUTPUT_ERROR
    PIPELINE_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PIPELINE_ERROR
    PREPARE_FOR_STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_STANDARDIZATION_ERROR
    PREPARE_FOR_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_VALIDATION_ERROR
    PROTONATION_CHANGED: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PROTONATION_CHANGED
    STANDARDIZATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STANDARDIZATION_ERROR
    STEREO_VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STEREO_VALIDATION_ERROR
    STRUCTURE_MODIFICATION: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STRUCTURE_MODIFICATION
    VALIDATION_ERROR: typing.ClassVar[PipelineStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.VALIDATION_ERROR
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'NO_EVENT': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NO_EVENT, 'INPUT_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.INPUT_ERROR, 'PREPARE_FOR_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_VALIDATION_ERROR, 'FEATURES_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FEATURES_VALIDATION_ERROR, 'BASIC_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.BASIC_VALIDATION_ERROR, 'IS2D_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.IS2D_VALIDATION_ERROR, 'LAYOUT2D_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.LAYOUT2D_VALIDATION_ERROR, 'STEREO_VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STEREO_VALIDATION_ERROR, 'VALIDATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.VALIDATION_ERROR, 'PREPARE_FOR_STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_STANDARDIZATION_ERROR, 'METAL_STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METAL_STANDARDIZATION_ERROR, 'NORMALIZER_STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZER_STANDARDIZATION_ERROR, 'FRAGMENT_STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENT_STANDARDIZATION_ERROR, 'CHARGE_STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.CHARGE_STANDARDIZATION_ERROR, 'STANDARDIZATION_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STANDARDIZATION_ERROR, 'OUTPUT_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.OUTPUT_ERROR, 'PIPELINE_ERROR': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PIPELINE_ERROR, 'METALS_DISCONNECTED': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METALS_DISCONNECTED, 'NORMALIZATION_APPLIED': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZATION_APPLIED, 'FRAGMENTS_REMOVED': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENTS_REMOVED, 'PROTONATION_CHANGED': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PROTONATION_CHANGED, 'STRUCTURE_MODIFICATION': rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STRUCTURE_MODIFICATION}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NO_EVENT, 1: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.INPUT_ERROR, 2: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_VALIDATION_ERROR, 4: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FEATURES_VALIDATION_ERROR, 8: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.BASIC_VALIDATION_ERROR, 16: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.IS2D_VALIDATION_ERROR, 32: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.LAYOUT2D_VALIDATION_ERROR, 64: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STEREO_VALIDATION_ERROR, 124: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.VALIDATION_ERROR, 128: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PREPARE_FOR_STANDARDIZATION_ERROR, 256: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METAL_STANDARDIZATION_ERROR, 512: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZER_STANDARDIZATION_ERROR, 1024: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENT_STANDARDIZATION_ERROR, 2048: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.CHARGE_STANDARDIZATION_ERROR, 3840: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STANDARDIZATION_ERROR, 4096: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.OUTPUT_ERROR, 8191: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PIPELINE_ERROR, 8388608: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.METALS_DISCONNECTED, 16777216: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.NORMALIZATION_APPLIED, 33554432: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.FRAGMENTS_REMOVED, 67108864: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.PROTONATION_CHANGED, 125829120: rdkit.Chem.MolStandardize.rdMolStandardize.PipelineStatus.STRUCTURE_MODIFICATION}
class RDKitValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, allowEmptyMolecules: bool = False) -> None:
        ...
    @property
    def allowEmptyMolecules(self) -> bool:
        """default: False"""
    @allowEmptyMolecules.setter
    def allowEmptyMolecules(self, value: bool) -> None: ...
class Reionizer(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, acidbaseFile: str) -> None:
        ...
    @typing.overload
    def __init__(self, acidbaseFile: str, ccs: typing.Any) -> None:
        ...
    def reionize(self, mol: Mol) -> rdkit.Chem.Mol:
        """
        """
    def reionizeInPlace(self, mol: Mol) -> None:
        """
            modifies the input molecule
        
        """
class SmilesTautomerMap(Boost.Python.instance):
    """
    maps SMILES strings to the respective Tautomer objects
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
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def items(self) -> tuple:
        """
        """
    def keys(self) -> tuple:
        """
        """
    def values(self) -> tuple:
        """
        """
class StereoValidation(ValidationMethod):
    __instance_size__: typing.ClassVar[int] = 32
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
class SubstructTerm(Boost.Python.instance):
    """
    Sets the score of this particular tautomer substructure, higher scores are more preferable
    Aromatic rings score 100, all carbon aromatic rings score 250
    """
    __instance_size__: typing.ClassVar[int] = 624
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, name: str, smarts: str, score: int) -> None:
        ...
    @property
    def name(*args, **kwargs):
        ...
    @property
    def score(*args, **kwargs):
        ...
    @property
    def smarts(*args, **kwargs):
        ...
class SubstructTermVector(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 48
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __contains__(self, item: typing.Any) -> bool:
        """
        """
    def __delitem__(self, item: typing.Any) -> None:
        """
        """
    def __getitem__(self, item: typing.Any) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    def __setitem__(self, item: typing.Any, value: typing.Any) -> None:
        """
        """
    def append(self, item: typing.Any) -> None:
        """
        """
    def extend(self, other: typing.Any) -> None:
        """
        """
class Tautomer(Boost.Python.instance):
    """
    used to hold the aromatic and kekulized versions of each tautomer
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
    def kekulized(*args, **kwargs):
        """
        kekulized version of the tautomer
        """
    @property
    def tautomer(*args, **kwargs):
        """
        aromatic version of the tautomer
        """
class TautomerEnumerator(Boost.Python.instance):
    tautomerScoreVersion: typing.ClassVar[str] = '1.0.0'
    @staticmethod
    def ScoreTautomer(mol: Mol) -> int:
        """
            returns the score for a tautomer using the default scoring scheme.
        
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def Canonicalize(self, mol: Mol) -> rdkit.Chem.Mol:
        """
            Returns the canonical tautomer for a molecule.
            
              The default scoring scheme is inspired by the publication:
              M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
              https://doi.org/10.1007/s10822-010-9346-4
            
              Note that the canonical tautomer is very likely not the most stable tautomer
              for any given conditions. The default scoring rules are designed to produce
              "reasonable" tautomers, but the primary concern is that the results are
              canonical: you always get the same canonical tautomer for a molecule
              regardless of what the input tautomer or atom ordering were.
        
        """
    @typing.overload
    def Canonicalize(self, mol: Mol, scoreFunc: typing.Any) -> rdkit.Chem.Mol:
        """
            picks the canonical tautomer from an iterable of molecules using a custom scoring function
        
        """
    def Enumerate(self, mol: Mol) -> TautomerEnumeratorResult:
        """
            Generates the tautomers for a molecule.
                         
              The enumeration rules are inspired by the publication:
              M. Sitzmann et al., “Tautomerism in Large Databases.”, JCAMD 24:521 (2010)
              https://doi.org/10.1007/s10822-010-9346-4
              
              Note: the definitions used here are that the atoms modified during
              tautomerization are the atoms at the beginning and end of each tautomer
              transform (the H "donor" and H "acceptor" in the transform) and the bonds
              modified during transformation are any bonds whose order is changed during
              the tautomer transform (these are the bonds between the "donor" and the
              "acceptor").
        
        """
    def GetCallback(self) -> typing.Any:
        """
            Get the TautomerEnumeratorCallback subclass instance,
            or None if none was set.
        
        """
    def GetMaxTautomers(self) -> int:
        """
            returns the maximum number of tautomers to be generated.
        
        """
    def GetMaxTransforms(self) -> int:
        """
            returns the maximum number of transformations to be applied.
        
        """
    def GetReassignStereo(self) -> bool:
        """
            returns whether AssignStereochemistry will be called on each tautomer generated by the Enumerate() method.
        
        """
    def GetRemoveBondStereo(self) -> bool:
        """
            returns whether stereochemistry information will be removed from double bonds involved in tautomerism.
        
        """
    def GetRemoveSp3Stereo(self) -> bool:
        """
            returns whether stereochemistry information will be removed from sp3 atoms involved in tautomerism.
        
        """
    @typing.overload
    def PickCanonical(self, iterable: typing.Any) -> rdkit.Chem.Mol:
        """
            picks the canonical tautomer from an iterable of molecules
        
        """
    @typing.overload
    def PickCanonical(self, iterable: typing.Any, scoreFunc: typing.Any) -> rdkit.Chem.Mol:
        """
            returns the canonical tautomer for a molecule using a custom scoring function
        
        """
    def SetCallback(self, callback: typing.Any) -> None:
        """
            Pass an instance of a class derived from
            TautomerEnumeratorCallback, which must implement the
            __call__() method.
        
        """
    def SetMaxTautomers(self, maxTautomers: int) -> None:
        """
            set the maximum number of tautomers to be generated.
        
        """
    def SetMaxTransforms(self, maxTransforms: int) -> None:
        """
            set the maximum number of transformations to be applied. This limit is usually hit earlier than the maxTautomers limit and leads to a more linear scaling of CPU time with increasing number of tautomeric centers (see Sitzmann et al.).
        
        """
    def SetReassignStereo(self, reassignStereo: bool) -> None:
        """
            set to True if you wish AssignStereochemistry to be called on each tautomer generated by the Enumerate() method. This defaults to True.
        
        """
    def SetRemoveBondStereo(self, removeBondStereo: bool) -> None:
        """
            set to True if you wish stereochemistry information to be removed from double bonds involved in tautomerism. This means that enols will lose their E/Z stereochemistry after going through tautomer enumeration because of the keto-enolic tautomerism. This defaults to True in the RDKit and also in the workflow described by Sitzmann et al.
        
        """
    def SetRemoveSp3Stereo(self, removeSp3Stereo: bool) -> None:
        """
            set to True if you wish stereochemistry information to be removed from sp3 atoms involved in tautomerism. This means that S-aminoacids will lose their stereochemistry after going through tautomer enumeration because of the amido-imidol tautomerism. This defaults to True in RDKit, and to False in the workflow described by Sitzmann et al.
        
        """
    @typing.overload
    def __init__(self) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, arg1: CleanupParameters) -> typing.Any:
        ...
    @typing.overload
    def __init__(self, arg1: TautomerEnumerator) -> typing.Any:
        ...
class TautomerEnumeratorCallback(Boost.Python.instance):
    """
    Create a derived class from this abstract base class and
        implement the __call__() method.
        The __call__() method is called in the innermost loop of the
        algorithm, and provides a mechanism to monitor or stop
        its progress.
    
        To have your callback called, pass an instance of your
        derived class to TautomerEnumerator.SetCallback()
    """
    __instance_size__: typing.ClassVar[int] = 56
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @typing.overload
    def __call__(self, mol: Mol, res: typing.Any) -> bool:
        """
            This must be implemented in the derived class. Return True if the tautomer enumeration should continue; False if the tautomer enumeration should stop.
            
        
        """
    @typing.overload
    def __call__(self, arg1: Mol, arg2: typing.Any) -> None:
        """
        """
    def __init__(self) -> None:
        ...
class TautomerEnumeratorResult(Boost.Python.instance):
    """
    used to return tautomer enumeration results
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
    def __call__(self) -> MOL_SPTR_VECT:
        """
            tautomers generated by the enumerator
        
        """
    def __getitem__(self, pos: int) -> rdkit.Chem.Mol:
        """
        """
    def __iter__(self) -> typing.Any:
        """
        """
    def __len__(self) -> int:
        """
        """
    @property
    def modifiedAtoms(*args, **kwargs):
        """
        tuple of atom indices modified by the transforms
        """
    @property
    def modifiedBonds(*args, **kwargs):
        """
        tuple of bond indices modified by the transforms
        """
    @property
    def smiles(*args, **kwargs):
        """
        SMILES of tautomers generated by the enumerator
        """
    @property
    def smilesTautomerMap(*args, **kwargs):
        """
        dictionary mapping SMILES strings to the respective Tautomer objects
        """
    @property
    def status(*args, **kwargs):
        """
        whether the enumeration completed or not; see TautomerEnumeratorStatus for possible values
        """
    @property
    def tautomers(*args, **kwargs):
        """
        tautomers generated by the enumerator
        """
class TautomerEnumeratorStatus(Boost.Python.enum):
    Canceled: typing.ClassVar[TautomerEnumeratorStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled
    Completed: typing.ClassVar[TautomerEnumeratorStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed
    MaxTautomersReached: typing.ClassVar[TautomerEnumeratorStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached
    MaxTransformsReached: typing.ClassVar[TautomerEnumeratorStatus]  # value = rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached
    __slots__: typing.ClassVar[tuple] = tuple()
    names: typing.ClassVar[dict]  # value = {'Completed': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed, 'MaxTautomersReached': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached, 'MaxTransformsReached': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached, 'Canceled': rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled}
    values: typing.ClassVar[dict]  # value = {0: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Completed, 1: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTautomersReached, 2: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.MaxTransformsReached, 3: rdkit.Chem.MolStandardize.rdMolStandardize.TautomerEnumeratorStatus.Canceled}
class Uncharger(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self, canonicalOrder: bool = True, force: bool = False, protonationOnly: bool = False) -> None:
        ...
    def uncharge(self, mol: Mol) -> rdkit.Chem.Mol:
        """
        """
    def unchargeInPlace(self, mol: Mol) -> None:
        """
            modifies the input molecule
        
        """
class ValidationMethod(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 40
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
    def validate(self, mol: Mol, reportAllFailures: bool = False) -> list:
        """
        """
class map_indexing_suite_SmilesTautomerMap_entry(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 104
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    @staticmethod
    def __repr__(arg1: map_indexing_suite_SmilesTautomerMap_entry) -> typing.Any:
        """
        """
    def __init__(self) -> None:
        ...
    def data(self) -> Tautomer:
        """
        """
    def key(self) -> str:
        """
        """
def CHARGE_CORRECTIONS() -> typing.Any:
    """
    """
def CanonicalTautomer(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Returns the canonical tautomer for the molecule
    
    """
def ChargeParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        Returns the uncharged version of the largest fragment
    
    """
@typing.overload
def ChargeParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the charge parent in place
    
    """
@typing.overload
def ChargeParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the chargeparent in place for multiple molecules
    
    """
def Cleanup(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Standardizes a molecule
    
    """
@typing.overload
def CleanupInPlace(mol: Mol, params: typing.Any = None) -> None:
    """
        Standardizes a molecule in place
    
    """
@typing.overload
def CleanupInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None) -> None:
    """
        Standardizes multiple molecules in place
    
    """
def DisconnectOrganometallics(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Returns the molecule disconnected using the organometallics rules.
    
    """
def DisconnectOrganometallicsInPlace(mol: Mol, params: typing.Any = None) -> None:
    """
        Disconnects the molecule using the organometallics rules, modifies the input molecule
    
    """
def FragmentParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        Returns the largest fragment after doing a cleanup
    
    """
@typing.overload
def FragmentParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the largest fragment in place
    
    """
@typing.overload
def FragmentParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the largest fragment in place for multiple molecules
    
    """
def FragmentRemoverFromData(fragmentData: str, leave_last: bool = True, skip_if_all_match: bool = False) -> FragmentRemover:
    """
        creates a FragmentRemover from a string containing parameter data
    
    """
def GetDefaultTautomerScoreSubstructs() -> SubstructTermVector:
    """
        Return the default tautomer substructure scoring terms
    
    """
def GetV1TautomerEnumerator() -> TautomerEnumerator:
    """
        return a TautomerEnumerator using v1 of the enumeration rules
    
    """
def IsotopeParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        removes all isotopes specifications from the given molecule
    
    """
@typing.overload
def IsotopeParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the isotope parent in place
    
    """
@typing.overload
def IsotopeParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the isotope parent in place for multiple molecules
    
    """
def Normalize(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Applies a series of standard transformations to correct functional groups and recombine charges
    
    """
@typing.overload
def NormalizeInPlace(mol: Mol, params: typing.Any = None) -> None:
    """
        Applies a series of standard transformations to correct functional groups and recombine charges, modifies the input molecule
    
    """
@typing.overload
def NormalizeInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None) -> None:
    """
        Normalizes multiple molecules in place
    
    """
def NormalizerFromData(paramData: str, params: CleanupParameters) -> Normalizer:
    """
        creates a Normalizer from a string containing normalization SMARTS
    
    """
def NormalizerFromParams(params: CleanupParameters) -> Normalizer:
    """
        creates a Normalizer from CleanupParameters
    
    """
def Reionize(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Ensures the strongest acid groups are charged first
    
    """
@typing.overload
def ReionizeInPlace(mol: Mol, params: typing.Any = None) -> None:
    """
        Ensures the strongest acid groups are charged first, modifies the input molecule
    
    """
@typing.overload
def ReionizeInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None) -> None:
    """
        Reionizes multiple molecules in place
    
    """
def ReionizerFromData(paramData: str, chargeCorrections: typing.Any = []) -> Reionizer:
    """
        creates a reionizer from a string containing parameter data and a list of charge corrections
    
    """
def RemoveFragments(mol: Mol, params: typing.Any = None) -> rdkit.Chem.Mol:
    """
        Removes fragments from the molecule
    
    """
@typing.overload
def RemoveFragmentsInPlace(mol: Mol, params: typing.Any = None) -> None:
    """
        Removes fragments from the molecule, modifies the input molecule
    
    """
@typing.overload
def RemoveFragmentsInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None) -> None:
    """
        Removes fragments from multiple molecules in place
    
    """
def ScoreHeteroHs(mol: Mol) -> int:
    """
        scores the number of heteroHs of the tautomer for canonicalization
        This gives a negative penalty to hydrogens attached to S,P, Se and Te
    
    """
def ScoreRings(mol: Mol) -> int:
    """
        scores the ring system of the tautomer for canonicalization
        Aromatic rings score 100, all carbon aromatic rings score 250
    
    """
def ScoreSubstructs(mol: Mol, terms: SubstructTermVector) -> int:
    """
        scores the tautomer substructures
    
    """
def StandardizeSmiles(smiles: str) -> str:
    """
        Convenience function for standardizing a SMILES
    
    """
def StereoParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        Generates the largest fragment in place for multiple molecules
    
    """
@typing.overload
def StereoParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the stereo parent in place
    
    """
@typing.overload
def StereoParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the stereo parent in place for multiple molecules
    
    """
def SuperParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        Returns the super parent. The super parent is the fragment, charge, isotope, stereo, and tautomer parent of the molecule.
    
    """
@typing.overload
def SuperParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the super parent in place
    
    """
@typing.overload
def SuperParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the super parent in place for multiple molecules
    
    """
def TautomerParent(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> rdkit.Chem.Mol:
    """
        Returns the tautomer parent of a given molecule. The fragment parent is the standardized canonical tautomer of the molecule
    
    """
@typing.overload
def TautomerParentInPlace(mol: Mol, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the tautomer parent in place
    
    """
@typing.overload
def TautomerParentInPlace(mols: typing.Any, numThreads: int, params: typing.Any = None, skipStandardize: bool = False) -> None:
    """
        Generates the tautomer parent in place for multiple molecules
    
    """
def UpdateParamsFromJSON(params: CleanupParameters, json: str) -> None:
    """
        updates the cleanup parameters from the provided JSON string
    
    """
def ValidateSmiles(mol: str) -> list:
    """
    """
