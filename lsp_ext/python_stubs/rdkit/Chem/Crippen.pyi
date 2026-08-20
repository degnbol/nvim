# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Atom-based calculation of LogP and MR using Crippen's approach


    Reference:
      S. A. Wildman and G. M. Crippen *JCICS* _39_ 868-873 (1999)


"""
from __future__ import annotations
import numpy as numpy
import os as os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDConfig
__all__: list[str] = ['Chem', 'RDConfig', 'defaultPatternFileName', 'numpy', 'os', 'rdMolDescriptors']
def _Init():
    ...
def _ReadPatts(fileName):
    """
     *Internal Use Only*
    
        parses the pattern list from the data file
    
      
    """
def _pyGetAtomContribs(mol, patts = None, order = None, verbose = 0, force = 0):
    """
     *Internal Use Only*
    
        calculates atomic contributions to the LogP and MR values
    
        if the argument *force* is not set, we'll use the molecules stored
        _crippenContribs value when possible instead of re-calculating.
    
      **Note:** Changes here affect the version numbers of MolLogP and MolMR
        as well as the VSA descriptors in Chem.MolSurf
    
      
    """
def _pyMolLogP(inMol, patts = None, order = None, verbose = 0, addHs = 1):
    """
     DEPRECATED
      
    """
def _pyMolMR(inMol, patts = None, order = None, verbose = 0, addHs = 1):
    """
     DEPRECATED
      
    """
_patternOrder: list = list()
_smartsPatterns: dict = {}
defaultPatternFileName: str = '/Users/runner/work/rdkit-pypi/rdkit-pypi/build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/share/RDKit/Data/Crippen.txt'

# present at runtime, absent from the generated stub:
def MolLogP(*x, **y):
    r"""
    Wildman-Crippen LogP value

    Uses an atom-based scheme based on the values in the paper:
       S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

    **Arguments**

      - inMol: a molecule

      - addHs: (optional) toggles adding of Hs to the molecule for the calculation.
        If true, hydrogens will be added to the molecule and used in the calculation.
    """
def MolMR(*x, **y):
    r"""
    Wildman-Crippen MR value

    Uses an atom-based scheme based on the values in the paper:
       S. A. Wildman and G. M. Crippen JCICS 39 868-873 (1999)

    **Arguments**

      - inMol: a molecule

      - addHs: (optional) toggles adding of Hs to the molecule for the calculation.
        If true, hydrogens will be added to the molecule and used in the calculation.
    """
