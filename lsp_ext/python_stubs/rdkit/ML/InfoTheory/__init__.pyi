# fix_pybind_stubs: rdkit 2026.3.5 5beea910
"""
 Information Theory functionality

"""
from __future__ import annotations
from rdkit.ML.InfoTheory.rdInfoTheory import BitCorrMatGenerator
from rdkit.ML.InfoTheory.rdInfoTheory import InfoBitRanker
from rdkit.ML.InfoTheory.rdInfoTheory import InfoType
from .rdInfoTheory import *
__all__: list[str] = ['BIASCHISQUARE', 'BIASENTROPY', 'BitCorrMatGenerator', 'CHISQUARE', 'ENTROPY', 'InfoBitRanker', 'InfoType', 'rdInfoTheory']
BIASCHISQUARE: rdInfoTheory.InfoType  # value = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASCHISQUARE
BIASENTROPY: rdInfoTheory.InfoType  # value = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.BIASENTROPY
CHISQUARE: rdInfoTheory.InfoType  # value = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.CHISQUARE
ENTROPY: rdInfoTheory.InfoType  # value = rdkit.ML.InfoTheory.rdInfoTheory.InfoType.ENTROPY

# present at runtime, absent from the generated stub:
from rdkit.ML.InfoTheory.rdInfoTheory import ChiSquare as ChiSquare
from rdkit.ML.InfoTheory.rdInfoTheory import InfoEntropy as InfoEntropy
from rdkit.ML.InfoTheory.rdInfoTheory import InfoGain as InfoGain
from . import entropy as entropy
from . import rdInfoTheory as rdInfoTheory
