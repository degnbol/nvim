# fix_pybind_stubs: rdkit 2026.3.5
"""
Module containing RDKit functionality for querying molecules.
"""
from __future__ import annotations
import rdkit.Chem
import typing
__all__: list[str] = ['AAtomQueryAtom', 'AHAtomQueryAtom', 'AtomNumEqualsQueryAtom', 'AtomNumGreaterQueryAtom', 'AtomNumLessQueryAtom', 'ExplicitDegreeEqualsQueryAtom', 'ExplicitDegreeGreaterQueryAtom', 'ExplicitDegreeLessQueryAtom', 'ExplicitValenceEqualsQueryAtom', 'ExplicitValenceGreaterQueryAtom', 'ExplicitValenceLessQueryAtom', 'FormalChargeEqualsQueryAtom', 'FormalChargeGreaterQueryAtom', 'FormalChargeLessQueryAtom', 'HCountEqualsQueryAtom', 'HCountGreaterQueryAtom', 'HCountLessQueryAtom', 'HasBitVectPropWithValueQueryAtom', 'HasBoolPropWithValueQueryAtom', 'HasBoolPropWithValueQueryBond', 'HasChiralTagQueryAtom', 'HasDoublePropWithValueQueryAtom', 'HasDoublePropWithValueQueryBond', 'HasIntPropWithValueQueryAtom', 'HasIntPropWithValueQueryBond', 'HasPropQueryAtom', 'HasPropQueryBond', 'HasStringPropWithValueQueryAtom', 'HasStringPropWithValueQueryBond', 'HybridizationEqualsQueryAtom', 'HybridizationGreaterQueryAtom', 'HybridizationLessQueryAtom', 'InNRingsEqualsQueryAtom', 'InNRingsGreaterQueryAtom', 'InNRingsLessQueryAtom', 'IsAliphaticQueryAtom', 'IsAromaticQueryAtom', 'IsBridgeheadQueryAtom', 'IsInRingQueryAtom', 'IsUnsaturatedQueryAtom', 'IsotopeEqualsQueryAtom', 'IsotopeGreaterQueryAtom', 'IsotopeLessQueryAtom', 'MAtomQueryAtom', 'MHAtomQueryAtom', 'MassEqualsQueryAtom', 'MassGreaterQueryAtom', 'MassLessQueryAtom', 'MinRingSizeEqualsQueryAtom', 'MinRingSizeGreaterQueryAtom', 'MinRingSizeLessQueryAtom', 'MissingChiralTagQueryAtom', 'NonHydrogenDegreeEqualsQueryAtom', 'NonHydrogenDegreeGreaterQueryAtom', 'NonHydrogenDegreeLessQueryAtom', 'NumAliphaticHeteroatomNeighborsEqualsQueryAtom', 'NumAliphaticHeteroatomNeighborsGreaterQueryAtom', 'NumAliphaticHeteroatomNeighborsLessQueryAtom', 'NumHeteroatomNeighborsEqualsQueryAtom', 'NumHeteroatomNeighborsGreaterQueryAtom', 'NumHeteroatomNeighborsLessQueryAtom', 'NumRadicalElectronsEqualsQueryAtom', 'NumRadicalElectronsGreaterQueryAtom', 'NumRadicalElectronsLessQueryAtom', 'QAtomQueryAtom', 'QHAtomQueryAtom', 'ReplaceAtomWithQueryAtom', 'RingBondCountEqualsQueryAtom', 'RingBondCountGreaterQueryAtom', 'RingBondCountLessQueryAtom', 'TotalDegreeEqualsQueryAtom', 'TotalDegreeGreaterQueryAtom', 'TotalDegreeLessQueryAtom', 'TotalValenceEqualsQueryAtom', 'TotalValenceGreaterQueryAtom', 'TotalValenceLessQueryAtom', 'XAtomQueryAtom', 'XHAtomQueryAtom']
def AAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when AAtom is True.
    
    """
def AHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when AHAtom is True.
    
    """
def AtomNumEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
    
    """
def AtomNumGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where AtomNum is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def AtomNumLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where AtomNum is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def ExplicitDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
    
    """
def ExplicitDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def ExplicitDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def ExplicitValenceEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
    
    """
def ExplicitValenceGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def ExplicitValenceLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where ExplicitValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def FormalChargeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
    
    """
def FormalChargeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where FormalCharge is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def FormalChargeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where FormalCharge is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def HCountEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.
    
    """
def HCountGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where HCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def HCountLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where HCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def HasBitVectPropWithValueQueryAtom(propname: str, val: ExplicitBitVect, negate: bool = False, tolerance: float = 0) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' has the specified explicit bit vector value.  The Tolerance is the allowed Tanimoto difference
    
    """
def HasBoolPropWithValueQueryAtom(propname: str, val: bool, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' has the specified boolean value.
    
    """
def HasBoolPropWithValueQueryBond(propname: str, val: bool, negate: bool = False) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' has the specified boolean value.
    
    """
def HasChiralTagQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when HasChiralTag is True.
    
    """
def HasDoublePropWithValueQueryAtom(propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' has the specified value +- tolerance
    
    """
def HasDoublePropWithValueQueryBond(propname: str, val: float, negate: bool = False, tolerance: float = 0.0) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' has the specified value +- tolerance
    
    """
def HasIntPropWithValueQueryAtom(propname: str, val: int, negate: bool = False, tolerance: int = 0) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' has the specified int value.
    
    """
def HasIntPropWithValueQueryBond(propname: str, val: int, negate: bool = False, tolerance: int = 0) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' has the specified int value.
    
    """
def HasPropQueryAtom(propname: str, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' exists in the atom.
    
    """
@typing.overload
def HasPropQueryBond(propname: str, negate: bool = False) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' exists in the bond.
    
    """
@typing.overload
def HasPropQueryBond(propname: str, negate: bool = False) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' exists in the bond.
    
    """
@typing.overload
def HasPropQueryBond(propname: str, negate: bool = False) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' exists in the bond.
    
    """
def HasStringPropWithValueQueryAtom(propname: str, val: str, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches when the property 'propname' has the specified string value.
    
    """
def HasStringPropWithValueQueryBond(propname: str, val: str, negate: bool = False) -> rdkit.Chem.QueryBond:
    """
        Returns a QueryBond that matches when the property 'propname' has the specified string value.
    
    """
def HybridizationEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
    
    """
def HybridizationGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Hybridization is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def HybridizationLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Hybridization is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def InNRingsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
    
    """
def InNRingsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where InNRings is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def InNRingsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where InNRings is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def IsAliphaticQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when IsAliphatic is True.
    
    """
def IsAromaticQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when IsAromatic is True.
    
    """
def IsBridgeheadQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when IsBridgehead is True.
    
    """
def IsInRingQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when IsInRing is True.
    
    """
def IsUnsaturatedQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when IsUnsaturated is True.
    
    """
def IsotopeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
    
    """
def IsotopeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Isotope is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def IsotopeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Isotope is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def MAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when MAtom is True.
    
    """
def MHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when MHAtom is True.
    
    """
def MassEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.
    
    """
def MassGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Mass is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def MassLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where Mass is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def MinRingSizeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
    
    """
def MinRingSizeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where MinRingSize is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def MinRingSizeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where MinRingSize is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def MissingChiralTagQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when MissingChiralTag is True.
    
    """
def NonHydrogenDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
    
    """
def NonHydrogenDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NonHydrogenDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NonHydrogenDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumAliphaticHeteroatomNeighborsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
    
    """
def NumAliphaticHeteroatomNeighborsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumAliphaticHeteroatomNeighborsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumAliphaticHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumHeteroatomNeighborsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
    
    """
def NumHeteroatomNeighborsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumHeteroatomNeighborsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumHeteroatomNeighbors is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumRadicalElectronsEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
    
    """
def NumRadicalElectronsGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def NumRadicalElectronsLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where NumRadicalElectrons is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def QAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when QAtom is True.
    
    """
def QHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when QHAtom is True.
    
    """
def ReplaceAtomWithQueryAtom(mol: Mol, atom: Atom) -> rdkit.Chem.Atom:
    """
        Changes the given atom in the molecule to
        a query atom and returns the atom which can then be modified, for example
        with additional query constraints added.  The new atom is otherwise a copy
        of the old.
        If the atom already has a query, nothing will be changed.
    
    """
def RingBondCountEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
    
    """
def RingBondCountGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where RingBondCount is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def RingBondCountLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where RingBondCount is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def TotalDegreeEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
    
    """
def TotalDegreeGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalDegree is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def TotalDegreeLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalDegree is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def TotalValenceEqualsQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
    
    """
def TotalValenceGreaterQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalValence is equal to the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def TotalValenceLessQueryAtom(val: int, negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms where TotalValence is less than the target value.
        NOTE: the direction of comparison is reversed relative to the C++ API
    
    """
def XAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when XAtom is True.
    
    """
def XHAtomQueryAtom(negate: bool = False) -> rdkit.Chem.QueryAtom:
    """
        Returns a QueryAtom that matches atoms when XHAtom is True.
    
    """
