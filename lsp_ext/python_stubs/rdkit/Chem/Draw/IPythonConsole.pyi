# fix_pybind_stubs: rdkit 2026.3.5 5beea910
import rdkit.Chem.Draw.rdMolDraw2D
from _typeshed import Incomplete
from rdkit.Chem.Draw.rdMolDraw2D import drawOptions as drawOptions

molSize: tuple
highlightSubstructs: bool
kekulizeStructures: bool
highlightByReactant: bool
ipython_useSVG: bool
ipython_showProperties: bool
ipython_maxProperties: int
ipython_3d: bool
molSize_3d: tuple
drawing_type_3d: str
bgcolor_3d: str
def addMolToView(mol, view, confId: int = ..., drawAs: Incomplete | None = ...): ...
def drawMol3D(m, view: Incomplete | None = ..., confId: int = ..., drawAs: Incomplete | None = ..., bgColor: Incomplete | None = ..., size: Incomplete | None = ...): ...
def drawMols3D(mols, view: Incomplete | None = ..., confIds: Incomplete | None = ..., drawAs: Incomplete | None = ..., bgColor: Incomplete | None = ..., size: Incomplete | None = ..., removeHs: bool = ..., colors: tuple = ...): ...
def listToLists(lst): ...
def display_pil_image(img): ...
def ShowMols(mols, maxMols: int = ..., **kwargs): ...
def DrawMorganBit(mol, bitId, bitInfo, drawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions = ..., **kwargs): ...
def DrawMorganBits(*args, drawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions = ..., **kwargs): ...
def DrawRDKitBit(mol, bitId, bitInfo, drawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions = ..., **kwargs): ...
def DrawRDKitBits(*args, drawOptions: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions = ..., **kwargs): ...
def EnableSubstructMatchRendering(): ...
def InstallIPythonRenderer(): ...
def DisableSubstructMatchRendering(): ...
def UninstallIPythonRenderer(): ...

# present at runtime, absent from the generated stub:
from _io import BytesIO as BytesIO
from rdkit import Chem as Chem
from rdkit.Chem import Draw as Draw
from IPython.core.display import HTML as HTML
import IPython as IPython
from PIL import Image as Image
from rdkit.Chem.Draw import InteractiveRenderer as InteractiveRenderer
from PIL.PngImagePlugin import PngInfo as PngInfo
from IPython.core.display import SVG as SVG
import base64 as base64
import copy as copy
from IPython import display as display
import html as html
from rdkit.Chem import rdChemReactions as rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D
from rdkit.Chem import rdchem as rdchem
import warnings as warnings
