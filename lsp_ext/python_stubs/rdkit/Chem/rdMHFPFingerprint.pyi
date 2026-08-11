# fix_pybind_stubs: rdkit 2026.3.5
from __future__ import annotations
import typing
__all__: list[str] = ['MHFPEncoder']
class MHFPEncoder(Boost.Python.instance):
    __instance_size__: typing.ClassVar[int] = 96
    @staticmethod
    def Distance(a: _vectj, b: _vectj) -> float:
        """
        """
    @staticmethod
    def __reduce__(*args, **kwargs):
        ...
    def CreateShinglingFromMol(self, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.
        
        """
    def CreateShinglingFromSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> _vectNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE:
        """
            Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.
        
        """
    def EncodeMol(self, mol: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> typing.Sequence[int]:
        """
            Creates a MHFP vector from an RDKit Mol instance.
        
        """
    def EncodeMolsBulk(self, mols: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> typing.Any:
        """
            Creates a MHFP vector from a list of RDKit Mol instances.
        
        """
    def EncodeSECFPMol(self, smiles: Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect:
        """
            Creates a SECFP binary vector from an RDKit Mol instance.
        
        """
    def EncodeSECFPMolsBulk(self, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> typing.Any:
        """
            Creates a SECFP binary vector from a list of RDKit Mol instances.
        
        """
    def EncodeSECFPSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> ExplicitBitVect:
        """
            Creates a SECFP binary vector from a SMILES string.
        
        """
    def EncodeSECFPSmilesBulk(self, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1, length: int = 2048) -> typing.Any:
        """
            Creates a SECFP binary vector from a list of SMILES strings.
        
        """
    def EncodeSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> typing.Sequence[int]:
        """
            Creates a MHFP vector from a SMILES string.
        
        """
    def EncodeSmilesBulk(self, smiles: list, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = False, min_radius: int = 1) -> typing.Any:
        """
            Creates a MHFP vector from a list of SMILES strings.
        
        """
    def FromArray(self, vec: list) -> typing.Sequence[int]:
        """
            Creates a MHFP vector from a list of unsigned integers.
        
        """
    def FromStringArray(self, vec: list) -> typing.Sequence[int]:
        """
            Creates a MHFP vector from a list of arbitrary strings.
        
        """
    def __init__(self, n_permutations: int, seed: int) -> None:
        ...
