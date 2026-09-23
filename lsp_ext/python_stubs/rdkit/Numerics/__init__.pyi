"""
 A module for Numerics stuff	
 
"""
from __future__ import annotations
from .rdAlignment import *
__all__: list[str] = ['Alignment', 'rdAlignment']
Alignment = rdAlignment

# present at runtime, absent from the generated stub:
from . import rdAlignment as rdAlignment
