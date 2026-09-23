"""
 functions to match a bunch of fragment descriptors from a file

No user-servicable parts inside.  ;-)

"""
from __future__ import annotations
import os as os
from rdkit import Chem
from rdkit import RDConfig
__all__: list[str] = ['Chem', 'RDConfig', 'defaultPatternFileName', 'fn', 'fns', 'name', 'os']
def _CountMatches(mol, patt, unique = True):
    ...
def _LoadPatterns(fileName = None):
    ...
defaultPatternFileName: str = '/Users/runner/work/rdkit-pypi/rdkit-pypi/build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/share/RDKit/Data/FragmentDescriptors.csv'
fn = None
fns: list  # value = [('fr_C_O', <function _LoadPatterns.<locals>.<lambda> at 0x105677d80>), ('fr_C_O_noCOO', <function _LoadPatterns.<locals>.<lambda> at 0x105677c40>), ('fr_Al_OH', <function _LoadPatterns.<locals>.<lambda> at 0x105677ce0>), ('fr_Ar_OH', <function _LoadPatterns.<locals>.<lambda> at 0x105677ba0>), ('fr_methoxy', <function _LoadPatterns.<locals>.<lambda> at 0x105677b00>), ('fr_oxime', <function _LoadPatterns.<locals>.<lambda> at 0x105677a60>), ('fr_ester', <function _LoadPatterns.<locals>.<lambda> at 0x1056779c0>), ('fr_Al_COO', <function _LoadPatterns.<locals>.<lambda> at 0x105675440>), ('fr_Ar_COO', <function _LoadPatterns.<locals>.<lambda> at 0x105677920>), ('fr_COO', <function _LoadPatterns.<locals>.<lambda> at 0x105677880>), ('fr_COO2', <function _LoadPatterns.<locals>.<lambda> at 0x1056777e0>), ('fr_ketone', <function _LoadPatterns.<locals>.<lambda> at 0x105677740>), ('fr_ether', <function _LoadPatterns.<locals>.<lambda> at 0x1056776a0>), ('fr_phenol', <function _LoadPatterns.<locals>.<lambda> at 0x105677600>), ('fr_aldehyde', <function _LoadPatterns.<locals>.<lambda> at 0x105677560>), ('fr_quatN', <function _LoadPatterns.<locals>.<lambda> at 0x1056774c0>), ('fr_NH2', <function _LoadPatterns.<locals>.<lambda> at 0x105677420>), ('fr_NH1', <function _LoadPatterns.<locals>.<lambda> at 0x105677380>), ('fr_NH0', <function _LoadPatterns.<locals>.<lambda> at 0x1056772e0>), ('fr_Ar_N', <function _LoadPatterns.<locals>.<lambda> at 0x105677240>), ('fr_Ar_NH', <function _LoadPatterns.<locals>.<lambda> at 0x1056771a0>), ('fr_aniline', <function _LoadPatterns.<locals>.<lambda> at 0x105677100>), ('fr_Imine', <function _LoadPatterns.<locals>.<lambda> at 0x105677060>), ('fr_nitrile', <function _LoadPatterns.<locals>.<lambda> at 0x105676fc0>), ('fr_hdrzine', <function _LoadPatterns.<locals>.<lambda> at 0x105676f20>), ('fr_hdrzone', <function _LoadPatterns.<locals>.<lambda> at 0x105676e80>), ('fr_nitroso', <function _LoadPatterns.<locals>.<lambda> at 0x105676de0>), ('fr_N_O', <function _LoadPatterns.<locals>.<lambda> at 0x105676d40>), ('fr_nitro', <function _LoadPatterns.<locals>.<lambda> at 0x105676ca0>), ('fr_azo', <function _LoadPatterns.<locals>.<lambda> at 0x105676c00>), ('fr_diazo', <function _LoadPatterns.<locals>.<lambda> at 0x105676b60>), ('fr_azide', <function _LoadPatterns.<locals>.<lambda> at 0x105676ac0>), ('fr_amide', <function _LoadPatterns.<locals>.<lambda> at 0x105676a20>), ('fr_priamide', <function _LoadPatterns.<locals>.<lambda> at 0x105676980>), ('fr_amidine', <function _LoadPatterns.<locals>.<lambda> at 0x1056768e0>), ('fr_guanido', <function _LoadPatterns.<locals>.<lambda> at 0x105676840>), ('fr_Nhpyrrole', <function _LoadPatterns.<locals>.<lambda> at 0x1056767a0>), ('fr_imide', <function _LoadPatterns.<locals>.<lambda> at 0x10563ae80>), ('fr_isocyan', <function _LoadPatterns.<locals>.<lambda> at 0x10563b2e0>), ('fr_isothiocyan', <function _LoadPatterns.<locals>.<lambda> at 0x10563af20>), ('fr_thiocyan', <function _LoadPatterns.<locals>.<lambda> at 0x10563b060>), ('fr_halogen', <function _LoadPatterns.<locals>.<lambda> at 0x105638cc0>), ('fr_alkyl_halide', <function _LoadPatterns.<locals>.<lambda> at 0x1056396c0>), ('fr_sulfide', <function _LoadPatterns.<locals>.<lambda> at 0x105639800>), ('fr_SH', <function _LoadPatterns.<locals>.<lambda> at 0x105639760>), ('fr_C_S', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8b80>), ('fr_sulfone', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8c20>), ('fr_sulfonamd', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8cc0>), ('fr_prisulfonamd', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8d60>), ('fr_barbitur', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8e00>), ('fr_urea', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8ea0>), ('fr_term_acetylene', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8f40>), ('fr_imidazole', <function _LoadPatterns.<locals>.<lambda> at 0x1085c8fe0>), ('fr_furan', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9080>), ('fr_thiophene', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9120>), ('fr_thiazole', <function _LoadPatterns.<locals>.<lambda> at 0x1085c91c0>), ('fr_oxazole', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9260>), ('fr_pyridine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9300>), ('fr_piperdine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c93a0>), ('fr_piperzine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9440>), ('fr_morpholine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c94e0>), ('fr_lactam', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9580>), ('fr_lactone', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9620>), ('fr_tetrazole', <function _LoadPatterns.<locals>.<lambda> at 0x1085c96c0>), ('fr_epoxide', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9760>), ('fr_unbrch_alkane', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9800>), ('fr_bicyclic', <function _LoadPatterns.<locals>.<lambda> at 0x1085c98a0>), ('fr_benzene', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9940>), ('fr_phos_acid', <function _LoadPatterns.<locals>.<lambda> at 0x1085c99e0>), ('fr_phos_ester', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9a80>), ('fr_nitro_arom', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9b20>), ('fr_nitro_arom_nonortho', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9bc0>), ('fr_dihydropyridine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9c60>), ('fr_phenol_noOrthoHbond', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9d00>), ('fr_Al_OH_noTert', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9da0>), ('fr_benzodiazepine', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9e40>), ('fr_para_hydroxylation', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9ee0>), ('fr_allylic_oxid', <function _LoadPatterns.<locals>.<lambda> at 0x1085c9f80>), ('fr_aryl_methyl', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca020>), ('fr_Ndealkylation1', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca160>), ('fr_Ndealkylation2', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca0c0>), ('fr_alkyl_carbamate', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca200>), ('fr_ketone_Topliss', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca2a0>), ('fr_ArN', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca340>), ('fr_HOCCN', <function _LoadPatterns.<locals>.<lambda> at 0x1085ca3e0>)]
name: str = 'fr_HOCCN'

# present at runtime, absent from the generated stub:
def fr_Al_COO(mol, countUnique=..., pattern=...):
    r"""
    Number of aliphatic carboxylic acids
    """
def fr_Al_OH(mol, countUnique=..., pattern=...):
    r"""
    Number of aliphatic hydroxyl groups
    """
def fr_Al_OH_noTert(mol, countUnique=..., pattern=...):
    r"""
    Number of aliphatic hydroxyl groups excluding tert-OH
    """
def fr_ArN(mol, countUnique=..., pattern=...):
    r"""
    Number of N functional groups attached to aromatics
    """
def fr_Ar_COO(mol, countUnique=..., pattern=...):
    r"""
    Number of Aromatic carboxylic acide
    """
def fr_Ar_N(mol, countUnique=..., pattern=...):
    r"""
    Number of aromatic nitrogens
    """
def fr_Ar_NH(mol, countUnique=..., pattern=...):
    r"""
    Number of aromatic amines
    """
def fr_Ar_OH(mol, countUnique=..., pattern=...):
    r"""
    Number of aromatic hydroxyl groups
    """
def fr_COO(mol, countUnique=..., pattern=...):
    r"""
    Number of carboxylic acids
    """
def fr_COO2(mol, countUnique=..., pattern=...):
    r"""
    Number of carboxylic acids
    """
def fr_C_O(mol, countUnique=..., pattern=...):
    r"""
    Number of carbonyl O
    """
def fr_C_O_noCOO(mol, countUnique=..., pattern=...):
    r"""
    Number of carbonyl O, excluding COOH
    """
def fr_C_S(mol, countUnique=..., pattern=...):
    r"""
    Number of thiocarbonyl
    """
def fr_HOCCN(mol, countUnique=..., pattern=...):
    r"""
    Number of C(OH)CCN-Ctert-alkyl or  C(OH)CCNcyclic
    """
def fr_Imine(mol, countUnique=..., pattern=...):
    r"""
    Number of Imines
    """
def fr_NH0(mol, countUnique=..., pattern=...):
    r"""
    Number of Tertiary amines
    """
def fr_NH1(mol, countUnique=..., pattern=...):
    r"""
    Number of Secondary amines
    """
def fr_NH2(mol, countUnique=..., pattern=...):
    r"""
    Number of Primary amines
    """
def fr_N_O(mol, countUnique=..., pattern=...):
    r"""
    Number of hydroxylamine groups
    """
def fr_Ndealkylation1(mol, countUnique=..., pattern=...):
    r"""
    Number of XCCNR groups
    """
def fr_Ndealkylation2(mol, countUnique=..., pattern=...):
    r"""
    Number of tert-alicyclic amines (no heteroatoms, not quinine-like bridged N)
    """
def fr_Nhpyrrole(mol, countUnique=..., pattern=...):
    r"""
    Number of H-pyrrole nitrogens
    """
def fr_SH(mol, countUnique=..., pattern=...):
    r"""
    Number of thiol groups
    """
def fr_aldehyde(mol, countUnique=..., pattern=...):
    r"""
    Number of aldehydes
    """
def fr_alkyl_carbamate(mol, countUnique=..., pattern=...):
    r"""
    Number of alkyl carbamates (subject to hydrolysis)
    """
def fr_alkyl_halide(mol, countUnique=..., pattern=...):
    r"""
    Number of alkyl halides
    """
def fr_allylic_oxid(mol, countUnique=..., pattern=...):
    r"""
    Number of allylic oxidation sites excluding steroid dienone
    """
def fr_amide(mol, countUnique=..., pattern=...):
    r"""
    Number of amides
    """
def fr_amidine(mol, countUnique=..., pattern=...):
    r"""
    Number of amidine groups
    """
def fr_aniline(mol, countUnique=..., pattern=...):
    r"""
    Number of anilines
    """
def fr_aryl_methyl(mol, countUnique=..., pattern=...):
    r"""
    Number of aryl methyl sites for hydroxylation
    """
def fr_azide(mol, countUnique=..., pattern=...):
    r"""
    Number of azide groups
    """
def fr_azo(mol, countUnique=..., pattern=...):
    r"""
    Number of azo groups
    """
def fr_barbitur(mol, countUnique=..., pattern=...):
    r"""
    Number of barbiturate groups
    """
def fr_benzene(mol, countUnique=..., pattern=...):
    r"""
    Number of benzene rings
    """
def fr_benzodiazepine(mol, countUnique=..., pattern=...):
    r"""
    Number of benzodiazepines with no additional fused rings
    """
def fr_bicyclic(mol, countUnique=..., pattern=...):
    r"""
    Bicyclic
    """
def fr_diazo(mol, countUnique=..., pattern=...):
    r"""
    Number of diazo groups
    """
def fr_dihydropyridine(mol, countUnique=..., pattern=...):
    r"""
    Number of dihydropyridines
    """
def fr_epoxide(mol, countUnique=..., pattern=...):
    r"""
    Number of epoxide rings
    """
def fr_ester(mol, countUnique=..., pattern=...):
    r"""
    Number of esters
    """
def fr_ether(mol, countUnique=..., pattern=...):
    r"""
    Number of ether oxygens (including phenoxy)
    """
def fr_furan(mol, countUnique=..., pattern=...):
    r"""
    Number of furan rings
    """
def fr_guanido(mol, countUnique=..., pattern=...):
    r"""
    Number of guanidine groups
    """
def fr_halogen(mol, countUnique=..., pattern=...):
    r"""
    Number of halogens
    """
def fr_hdrzine(mol, countUnique=..., pattern=...):
    r"""
    Number of hydrazine groups
    """
def fr_hdrzone(mol, countUnique=..., pattern=...):
    r"""
    Number of hydrazone groups
    """
def fr_imidazole(mol, countUnique=..., pattern=...):
    r"""
    Number of imidazole rings
    """
def fr_imide(mol, countUnique=..., pattern=...):
    r"""
    Number of imide groups
    """
def fr_isocyan(mol, countUnique=..., pattern=...):
    r"""
    Number of isocyanates
    """
def fr_isothiocyan(mol, countUnique=..., pattern=...):
    r"""
    Number of isothiocyanates
    """
def fr_ketone(mol, countUnique=..., pattern=...):
    r"""
    Number of ketones
    """
def fr_ketone_Topliss(mol, countUnique=..., pattern=...):
    r"""
    Number of ketones excluding diaryl, a,b-unsat. dienones, heteroatom on Calpha
    """
def fr_lactam(mol, countUnique=..., pattern=...):
    r"""
    Number of beta lactams
    """
def fr_lactone(mol, countUnique=..., pattern=...):
    r"""
    Number of cyclic esters (lactones)
    """
def fr_methoxy(mol, countUnique=..., pattern=...):
    r"""
    Number of methoxy groups -OCH3
    """
def fr_morpholine(mol, countUnique=..., pattern=...):
    r"""
    Number of morpholine rings
    """
def fr_nitrile(mol, countUnique=..., pattern=...):
    r"""
    Number of nitriles
    """
def fr_nitro(mol, countUnique=..., pattern=...):
    r"""
    Number of nitro groups
    """
def fr_nitro_arom(mol, countUnique=..., pattern=...):
    r"""
    Number of nitro benzene ring substituents
    """
def fr_nitro_arom_nonortho(mol, countUnique=..., pattern=...):
    r"""
    Number of non-ortho nitro benzene ring substituents
    """
def fr_nitroso(mol, countUnique=..., pattern=...):
    r"""
    Number of nitroso groups, excluding NO2
    """
def fr_oxazole(mol, countUnique=..., pattern=...):
    r"""
    Number of oxazole rings
    """
def fr_oxime(mol, countUnique=..., pattern=...):
    r"""
    Number of oxime groups
    """
def fr_para_hydroxylation(mol, countUnique=..., pattern=...):
    r"""
    Number of para-hydroxylation sites
    """
def fr_phenol(mol, countUnique=..., pattern=...):
    r"""
    Number of phenols
    """
def fr_phenol_noOrthoHbond(mol, countUnique=..., pattern=...):
    r"""
    Number of phenolic OH excluding ortho intramolecular Hbond substituents
    """
def fr_phos_acid(mol, countUnique=..., pattern=...):
    r"""
    Number of phosphoric acid groups
    """
def fr_phos_ester(mol, countUnique=..., pattern=...):
    r"""
    Number of phosphoric ester groups
    """
def fr_piperdine(mol, countUnique=..., pattern=...):
    r"""
    Number of piperdine rings
    """
def fr_piperzine(mol, countUnique=..., pattern=...):
    r"""
    Number of piperzine rings
    """
def fr_priamide(mol, countUnique=..., pattern=...):
    r"""
    Number of primary amides
    """
def fr_prisulfonamd(mol, countUnique=..., pattern=...):
    r"""
    Number of primary sulfonamides
    """
def fr_pyridine(mol, countUnique=..., pattern=...):
    r"""
    Number of pyridine rings
    """
def fr_quatN(mol, countUnique=..., pattern=...):
    r"""
    Number of quaternary nitrogens
    """
def fr_sulfide(mol, countUnique=..., pattern=...):
    r"""
    Number of thioether
    """
def fr_sulfonamd(mol, countUnique=..., pattern=...):
    r"""
    Number of sulfonamides
    """
def fr_sulfone(mol, countUnique=..., pattern=...):
    r"""
    Number of sulfone groups
    """
def fr_term_acetylene(mol, countUnique=..., pattern=...):
    r"""
    Number of terminal acetylenes
    """
def fr_tetrazole(mol, countUnique=..., pattern=...):
    r"""
    Number of tetrazole rings
    """
def fr_thiazole(mol, countUnique=..., pattern=...):
    r"""
    Number of thiazole rings
    """
def fr_thiocyan(mol, countUnique=..., pattern=...):
    r"""
    Number of thiocyanates
    """
def fr_thiophene(mol, countUnique=..., pattern=...):
    r"""
    Number of thiophene rings
    """
def fr_unbrch_alkane(mol, countUnique=..., pattern=...):
    r"""
    Number of unbranched alkanes  of at least 4 members (excludes halogenated alkanes)
    """
def fr_urea(mol, countUnique=..., pattern=...):
    r"""
    Number of urea groups
    """
