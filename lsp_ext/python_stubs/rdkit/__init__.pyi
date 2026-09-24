# patch_stubs: rdkit 2026.3.6 58692ef9
from __future__ import annotations
import logging as logging
import sys as sys
from .rdBase import *
__all__: list[str] = ['VECT_WRAPS', 'VectIter', 'log_handler', 'logger', 'logging', 'name', 'object', 'rdBase', 'sys']
class VectIter:
    def __init__(self, vect):
        ...
    def __iter__(self):
        ...
    def __next__(self):
        ...
def __vect__iter__(vect):
    ...
VECT_WRAPS: set = {'VectorOfStringVectors', 'MatchTypeVect', 'VectSizeT', 'UnsignedLong_Vect'}
__version__: str = '2026.03.6'
log_handler: logging.StreamHandler  # value = <StreamHandler <stderr> (NOTSET)>
logger: logging.Logger  # value = <Logger rdkit (WARNING)>
name: str = '__file__'
object: str = '/Users/runner/work/rdkit-pypi/rdkit-pypi/build/temp.macosx-11.0-arm64-cpython-311/rdkit_install/lib/python3.11/site-packages/rdkit/rdBase.so'

# present at runtime, absent from the generated stub:
from . import Avalon as Avalon
from . import Chem as Chem
from . import DataManip as DataManip
from . import DataStructs as DataStructs
from . import Dbase as Dbase
from . import DistanceGeometry as DistanceGeometry
from . import ForceField as ForceField
from . import Geometry as Geometry
from . import ML as ML
from . import Numerics as Numerics
from . import RDConfig as RDConfig
from . import RDLogger as RDLogger
from . import RDPaths as RDPaths
from . import RDRandom as RDRandom
from . import SimDivFilters as SimDivFilters
from . import VLib as VLib
from . import rdBase as rdBase
from . import sping as sping
from . import utils as utils
