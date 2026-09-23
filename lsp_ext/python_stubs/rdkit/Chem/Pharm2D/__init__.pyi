"""
 module with functionality for 2D pharmacophores

"""
from __future__ import annotations
__all__: list[str] = ['DefaultSigFactory']
def DefaultSigFactory(fdefFile = None, minPointCount = 2, maxPointCount = 3, bins = [(2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 100)]):
    ...

# present at runtime, absent from the generated stub:
from . import Generate as Generate
from . import Gobbi_Pharm2D as Gobbi_Pharm2D
from . import Matcher as Matcher
from . import SigFactory as SigFactory
from . import Utils as Utils
