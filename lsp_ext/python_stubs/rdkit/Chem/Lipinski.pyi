# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Calculation of Lipinski parameters for molecules

"""
from __future__ import annotations
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import rdkit.Chem.rdchem
__all__: list[str] = ['Chem', 'HAcceptorSmarts', 'HDonorSmarts', 'HeavyAtomCount', 'HeteroatomSmarts', 'NHOHSmarts', 'NOCountSmarts', 'RotatableBondSmarts', 'nm', 'rdMolDescriptors', 'txt']
def HeavyAtomCount(mol):
    """
     Number of heavy atoms a molecule.
    """
HAcceptorSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
HDonorSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
HeteroatomSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
NHOHSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
NOCountSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
RotatableBondSmarts: rdkit.Chem.rdchem.Mol  # value = <rdkit.Chem.rdchem.Mol object>
_bulkConvert: tuple = ('CalcFractionCSP3', 'CalcNumAromaticRings', 'CalcNumSaturatedRings', 'CalcNumAromaticHeterocycles', 'CalcNumAromaticCarbocycles', 'CalcNumSaturatedHeterocycles', 'CalcNumSaturatedCarbocycles', 'CalcNumAliphaticRings', 'CalcNumAliphaticHeterocycles', 'CalcNumAliphaticCarbocycles', 'CalcNumHeterocycles', 'CalcNumBridgeheadAtoms', 'CalcNumAmideBonds', 'CalcNumAtomStereoCenters', 'CalcNumHeterocycles', 'CalcNumUnspecifiedAtomStereoCenters', 'CalcNumSpiroAtoms', 'CalcPhi')
nm: str = 'Phi'
txt: str = 'CalcPhi'

# present at runtime, absent from the generated stub:
def FractionCSP3(x, y=...):
    r"""
    CalcFractionCSP3( (Mol)mol) -> float :
        returns the fraction of C atoms that are SP3 hybridized

        C++ signature :
            double CalcFractionCSP3(RDKit::ROMol)
    """
def NHOHCount(x):
    r"""
    Number of NHs or OHs
    """
def NOCount(x):
    r"""
    Number of Nitrogens and Oxygens
    """
def NumAliphaticCarbocycles(x, y=...):
    r"""
    CalcNumAliphaticCarbocycles( (Mol)mol) -> int :
        returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticCarbocycles(RDKit::ROMol)
    """
def NumAliphaticHeterocycles(x, y=...):
    r"""
    CalcNumAliphaticHeterocycles( (Mol)mol) -> int :
        returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticHeterocycles(RDKit::ROMol)
    """
def NumAliphaticRings(x, y=...):
    r"""
    CalcNumAliphaticRings( (Mol)mol) -> int :
        returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule

        C++ signature :
            unsigned int CalcNumAliphaticRings(RDKit::ROMol)
    """
def NumAmideBonds(x, y=...):
    r"""
    CalcNumAmideBonds( (Mol)mol) -> int :
        returns the number of amide bonds in a molecule

        C++ signature :
            unsigned int CalcNumAmideBonds(RDKit::ROMol)
    """
def NumAromaticCarbocycles(x, y=...):
    r"""
    CalcNumAromaticCarbocycles( (Mol)mol) -> int :
        returns the number of aromatic carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumAromaticCarbocycles(RDKit::ROMol)
    """
def NumAromaticHeterocycles(x, y=...):
    r"""
    CalcNumAromaticHeterocycles( (Mol)mol) -> int :
        returns the number of aromatic heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumAromaticHeterocycles(RDKit::ROMol)
    """
def NumAromaticRings(x, y=...):
    r"""
    CalcNumAromaticRings( (Mol)mol) -> int :
        returns the number of aromatic rings for a molecule

        C++ signature :
            unsigned int CalcNumAromaticRings(RDKit::ROMol)
    """
def NumAtomStereoCenters(x, y=...):
    r"""
    CalcNumAtomStereoCenters( (Mol)mol) -> int :
        Returns the total number of atomic stereocenters (specified and unspecified)

        C++ signature :
            unsigned int CalcNumAtomStereoCenters(RDKit::ROMol)
    """
def NumBridgeheadAtoms(x, y=...):
    r"""
    CalcNumBridgeheadAtoms( (Mol)mol [, (AtomPairsParameters)atoms=None]) -> int :
        Returns the number of bridgehead atoms (atoms shared between rings that share at least two bonds)

        C++ signature :
            unsigned int CalcNumBridgeheadAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
def NumHAcceptors(x):
    r"""
    Number of Hydrogen Bond Acceptors
    """
def NumHDonors(x):
    r"""
    Number of Hydrogen Bond Donors
    """
def NumHeteroatoms(x):
    r"""
    Number of Heteroatoms
    """
def NumHeterocycles(x, y=...):
    r"""
    CalcNumHeterocycles( (Mol)mol) -> int :
        returns the number of heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumHeterocycles(RDKit::ROMol)
    """
def NumRotatableBonds(x):
    r"""
    Number of Rotatable Bonds
    """
def NumSaturatedCarbocycles(x, y=...):
    r"""
    CalcNumSaturatedCarbocycles( (Mol)mol) -> int :
        returns the number of saturated carbocycles for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedCarbocycles(RDKit::ROMol)
    """
def NumSaturatedHeterocycles(x, y=...):
    r"""
    CalcNumSaturatedHeterocycles( (Mol)mol) -> int :
        returns the number of saturated heterocycles for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedHeterocycles(RDKit::ROMol)
    """
def NumSaturatedRings(x, y=...):
    r"""
    CalcNumSaturatedRings( (Mol)mol) -> int :
        returns the number of saturated rings for a molecule

        C++ signature :
            unsigned int CalcNumSaturatedRings(RDKit::ROMol)
    """
def NumSpiroAtoms(x, y=...):
    r"""
    CalcNumSpiroAtoms( (Mol)mol [, (AtomPairsParameters)atoms=None]) -> int :
        Returns the number of spiro atoms (atoms shared between rings that share exactly one atom)

        C++ signature :
            unsigned int CalcNumSpiroAtoms(RDKit::ROMol [,boost::python::api::object=None])
    """
def NumUnspecifiedAtomStereoCenters(x, y=...):
    r"""
    CalcNumUnspecifiedAtomStereoCenters( (Mol)mol) -> int :
        Returns the number of unspecified atomic stereocenters

        C++ signature :
            unsigned int CalcNumUnspecifiedAtomStereoCenters(RDKit::ROMol)
    """
def Phi(x, y=...):
    r"""
    CalcPhi( (Mol)mol) -> float :

        C++ signature :
            double CalcPhi(RDKit::ROMol)
    """
def RingCount(x):
    ...
