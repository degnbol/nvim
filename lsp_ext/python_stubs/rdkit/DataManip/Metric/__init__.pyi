# fix_pybind_stubs: rdkit 2026.3.5 5beea910
from __future__ import annotations
from rdkit import rdBase
from .rdMetricMatrixCalc import *
__all__: list[str] = ['rdBase', 'rdMetricMatrixCalc']

# present at runtime, absent from the generated stub:
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetEuclideanDistMat as GetEuclideanDistMat
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoDistMat as GetTanimotoDistMat
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat as GetTanimotoSimMat
from . import rdMetricMatrixCalc as rdMetricMatrixCalc
