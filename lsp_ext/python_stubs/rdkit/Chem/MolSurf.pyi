# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Exposes functionality for MOE-like approximate molecular surface area
descriptors.

  The MOE-like VSA descriptors are also calculated here

"""
from __future__ import annotations
import bisect as bisect
import numpy as numpy
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdPartialCharges
import rdkit.Chem.rdchem
__all__: list[str] = ['Chem', 'Crippen', 'bisect', 'bondScaleFacts', 'chgBins', 'logpBins', 'mrBins', 'numpy', 'ptable', 'pyLabuteASA', 'pyPEOE_VSA_', 'pySMR_VSA_', 'pySlogP_VSA_', 'rdMolDescriptors', 'rdPartialCharges']
def _InstallDescriptors():
    ...
def _LabuteHelper(mol, includeHs = 1, force = 0):
    """
     *Internal Use Only*
        helper function for LabuteASA calculation
        returns an array of atomic contributions to the ASA
    
      **Note:** Changes here affect the version numbers of all ASA descriptors
    
      
    """
def _pyLabuteHelper(mol, includeHs = 1, force = 0):
    """
     *Internal Use Only*
        helper function for LabuteASA calculation
        returns an array of atomic contributions to the ASA
    
      **Note:** Changes here affect the version numbers of all ASA descriptors
    
      
    """
def _pyTPSA(mol, verbose = False):
    """
     DEPRECATED: this has been reimplmented in C++
       calculates the polar surface area of a molecule based upon fragments
    
       Algorithm in:
        P. Ertl, B. Rohde, P. Selzer
         Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
         Contributions and Its Application to the Prediction of Drug Transport
         Properties, J.Med.Chem. 43, 3714-3717, 2000
    
       Implementation based on the Daylight contrib program tpsa.c
      
    """
def _pyTPSAContribs(mol, verbose = False):
    """
     DEPRECATED: this has been reimplmented in C++
      calculates atomic contributions to a molecules TPSA
    
       Algorithm described in:
        P. Ertl, B. Rohde, P. Selzer
         Fast Calculation of Molecular Polar Surface Area as a Sum of Fragment-based
         Contributions and Its Application to the Prediction of Drug Transport
         Properties, J.Med.Chem. 43, 3714-3717, 2000
    
       Implementation based on the Daylight contrib program tpsa.c
    
       NOTE: The JMC paper describing the TPSA algorithm includes
       contributions from sulfur and phosphorus, however according to
       Peter Ertl (personal communication, 2010) the correlation of TPSA
       with various ADME properties is better if only contributions from
       oxygen and nitrogen are used. This matches the daylight contrib
       implementation.
    
      
    """
def pyLabuteASA(mol, includeHs = 1):
    """
     calculates Labute's Approximate Surface Area (ASA from MOE)
    
        Definition from P. Labute's article in the Journal of the Chemical Computing Group
        and J. Mol. Graph. Mod.  _18_ 464-477 (2000)
    
      
    """
def pyPEOE_VSA_(mol, bins = None, force = 1):
    """
     *Internal Use Only*
      
    """
def pySMR_VSA_(mol, bins = None, force = 1):
    """
     *Internal Use Only*
      
    """
def pySlogP_VSA_(mol, bins = None, force = 1):
    """
     *Internal Use Only*
      
    """
bondScaleFacts: list = [0.1, 0, 0.2, 0.3]
chgBins: list = [-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
logpBins: list = [-0.4, -0.2, 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
mrBins: list = [1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63, 3.8, 4.0]
ptable: rdkit.Chem.rdchem.PeriodicTable  # value = <rdkit.Chem.rdchem.PeriodicTable object>

# present at runtime, absent from the generated stub:
def LabuteASA(*x, **y):
    ...
def PEOE_VSA1(x, y=...):
    r"""
    MOE Charge VSA Descriptor 1 (-inf < x < -0.30)
    """
def PEOE_VSA10(x, y=...):
    r"""
    MOE Charge VSA Descriptor 10 ( 0.10 <= x <  0.15)
    """
def PEOE_VSA11(x, y=...):
    r"""
    MOE Charge VSA Descriptor 11 ( 0.15 <= x <  0.20)
    """
def PEOE_VSA12(x, y=...):
    r"""
    MOE Charge VSA Descriptor 12 ( 0.20 <= x <  0.25)
    """
def PEOE_VSA13(x, y=...):
    r"""
    MOE Charge VSA Descriptor 13 ( 0.25 <= x <  0.30)
    """
def PEOE_VSA14(x, y=...):
    r"""
    MOE Charge VSA Descriptor 14 ( 0.30 <= x < inf)
    """
def PEOE_VSA2(x, y=...):
    r"""
    MOE Charge VSA Descriptor 2 (-0.30 <= x < -0.25)
    """
def PEOE_VSA3(x, y=...):
    r"""
    MOE Charge VSA Descriptor 3 (-0.25 <= x < -0.20)
    """
def PEOE_VSA4(x, y=...):
    r"""
    MOE Charge VSA Descriptor 4 (-0.20 <= x < -0.15)
    """
def PEOE_VSA5(x, y=...):
    r"""
    MOE Charge VSA Descriptor 5 (-0.15 <= x < -0.10)
    """
def PEOE_VSA6(x, y=...):
    r"""
    MOE Charge VSA Descriptor 6 (-0.10 <= x < -0.05)
    """
def PEOE_VSA7(x, y=...):
    r"""
    MOE Charge VSA Descriptor 7 (-0.05 <= x <  0.00)
    """
def PEOE_VSA8(x, y=...):
    r"""
    MOE Charge VSA Descriptor 8 ( 0.00 <= x <  0.05)
    """
def PEOE_VSA9(x, y=...):
    r"""
    MOE Charge VSA Descriptor 9 ( 0.05 <= x <  0.10)
    """
from rdkit.Chem.rdMolDescriptors import PEOE_VSA_ as PEOE_VSA_
def SMR_VSA1(x, y=...):
    r"""
    MOE MR VSA Descriptor 1 (-inf < x <  1.29)
    """
def SMR_VSA10(x, y=...):
    r"""
    MOE MR VSA Descriptor 10 ( 4.00 <= x < inf)
    """
def SMR_VSA2(x, y=...):
    r"""
    MOE MR VSA Descriptor 2 ( 1.29 <= x <  1.82)
    """
def SMR_VSA3(x, y=...):
    r"""
    MOE MR VSA Descriptor 3 ( 1.82 <= x <  2.24)
    """
def SMR_VSA4(x, y=...):
    r"""
    MOE MR VSA Descriptor 4 ( 2.24 <= x <  2.45)
    """
def SMR_VSA5(x, y=...):
    r"""
    MOE MR VSA Descriptor 5 ( 2.45 <= x <  2.75)
    """
def SMR_VSA6(x, y=...):
    r"""
    MOE MR VSA Descriptor 6 ( 2.75 <= x <  3.05)
    """
def SMR_VSA7(x, y=...):
    r"""
    MOE MR VSA Descriptor 7 ( 3.05 <= x <  3.63)
    """
def SMR_VSA8(x, y=...):
    r"""
    MOE MR VSA Descriptor 8 ( 3.63 <= x <  3.80)
    """
def SMR_VSA9(x, y=...):
    r"""
    MOE MR VSA Descriptor 9 ( 3.80 <= x <  4.00)
    """
from rdkit.Chem.rdMolDescriptors import SMR_VSA_ as SMR_VSA_
def SlogP_VSA1(x, y=...):
    r"""
    MOE logP VSA Descriptor 1 (-inf < x < -0.40)
    """
def SlogP_VSA10(x, y=...):
    r"""
    MOE logP VSA Descriptor 10 ( 0.40 <= x <  0.50)
    """
def SlogP_VSA11(x, y=...):
    r"""
    MOE logP VSA Descriptor 11 ( 0.50 <= x <  0.60)
    """
def SlogP_VSA12(x, y=...):
    r"""
    MOE logP VSA Descriptor 12 ( 0.60 <= x < inf)
    """
def SlogP_VSA2(x, y=...):
    r"""
    MOE logP VSA Descriptor 2 (-0.40 <= x < -0.20)
    """
def SlogP_VSA3(x, y=...):
    r"""
    MOE logP VSA Descriptor 3 (-0.20 <= x <  0.00)
    """
def SlogP_VSA4(x, y=...):
    r"""
    MOE logP VSA Descriptor 4 ( 0.00 <= x <  0.10)
    """
def SlogP_VSA5(x, y=...):
    r"""
    MOE logP VSA Descriptor 5 ( 0.10 <= x <  0.15)
    """
def SlogP_VSA6(x, y=...):
    r"""
    MOE logP VSA Descriptor 6 ( 0.15 <= x <  0.20)
    """
def SlogP_VSA7(x, y=...):
    r"""
    MOE logP VSA Descriptor 7 ( 0.20 <= x <  0.25)
    """
def SlogP_VSA8(x, y=...):
    r"""
    MOE logP VSA Descriptor 8 ( 0.25 <= x <  0.30)
    """
def SlogP_VSA9(x, y=...):
    r"""
    MOE logP VSA Descriptor 9 ( 0.30 <= x <  0.40)
    """
from rdkit.Chem.rdMolDescriptors import SlogP_VSA_ as SlogP_VSA_
def TPSA(*x, **y):
    ...
