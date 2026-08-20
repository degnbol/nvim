# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Descriptors derived from a molecule's 3D structure

"""
from __future__ import annotations
from rdkit.Chem.Descriptors import _isCallable
from rdkit.Chem import rdMolDescriptors
__all__: list[str] = ['CalcMolDescriptors3D', 'descList', 'rdMolDescriptors']
def CalcMolDescriptors3D(mol, confId = -1):
    """
    
        Compute all 3D descriptors of a molecule
        
        Arguments:
        - mol: the molecule to work with
        - confId: conformer ID to work with. If not specified the default (-1) is used
        
        Return:
        
        dict
            A dictionary with decriptor names as keys and the descriptor values as values
    
        raises a ValueError 
            If the molecule does not have conformers
        
    """
def _setupDescriptors(namespace):
    ...
descList: list  # value = [('PMI1', <function <lambda> at 0x104ae7d80>), ('PMI2', <function <lambda> at 0x104ae76a0>), ('PMI3', <function <lambda> at 0x104ae7560>), ('NPR1', <function <lambda> at 0x104ae74c0>), ('NPR2', <function <lambda> at 0x104ae7420>), ('RadiusOfGyration', <function <lambda> at 0x104ae7380>), ('InertialShapeFactor', <function <lambda> at 0x104ae72e0>), ('Eccentricity', <function <lambda> at 0x104ae7240>), ('Asphericity', <function <lambda> at 0x104ae71a0>), ('SpherocityIndex', <function <lambda> at 0x104ae7100>), ('PBF', <function <lambda> at 0x104ae7060>)]

# present at runtime, absent from the generated stub:
def Asphericity(*x, **y):
    r"""
    molecular asphericity

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       https://doi.org/10.1002/9783527618279.ch37

       Definition:
         0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def Eccentricity(*x, **y):
    r"""
    molecular eccentricity

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       https://doi.org/10.1002/9783527618279.ch37

       Definition:
         sqrt(pm3**2 -pm1**2) / pm3

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def InertialShapeFactor(*x, **y):
    r"""
    Inertial shape factor

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       https://doi.org/10.1002/9783527618279.ch37

       Definition:
         pm2 / (pm1*pm3)

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def NPR1(*x, **y):
    r"""
    Normalized principal moments ratio 1 (=I1/I3)

        from Sauer and Schwarz JCIM 43:987-1003 (2003)
        https://doi.org/10.1021/ci025599w


    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def NPR2(*x, **y):
    r"""
    Normalized principal moments ratio 2 (=I2/I3)

        from Sauer and Schwarz JCIM 43:987-1003 (2003)
        https://doi.org/10.1021/ci025599w


    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def PBF(*x, **y):
    r"""
    Plane of best fit
      from: https://doi.org/10.1021/ci300293f

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use
    """
def PMI1(*x, **y):
    r"""
    First (smallest) principal moment of inertia


    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def PMI2(*x, **y):
    r"""
    Second principal moment of inertia

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def PMI3(*x, **y):
    r"""
    Third (largest) principal moment of inertia

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def RadiusOfGyration(*x, **y):
    r"""
    Radius of gyration

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       https://doi.org/10.1002/9783527618279.ch37

       Definition:
         for planar molecules: sqrt( sqrt(pm3*pm2)/MW )
         for nonplanar molecules: sqrt( 2*pi*pow(pm3*pm2*pm1,1/3)/MW )

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use

      - useAtomicMasses: (optional) toggles use of atomic masses in the
        calculation. Defaults to True
    """
def SpherocityIndex(*x, **y):
    r"""
    Molecular spherocityIndex

       from Todeschini and Consoni "Descriptors from Molecular Geometry"
       Handbook of Chemoinformatics
       https://doi.org/10.1002/9783527618279.ch37

       Definition:
         3 * pm1 / (pm1+pm2+pm3) where the moments are calculated without weights

    **Arguments**

      - inMol: a molecule

      - confId: (optional) the conformation ID to use
    """
