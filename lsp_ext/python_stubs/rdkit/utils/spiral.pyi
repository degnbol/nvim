from __future__ import annotations
import math as math
import numpy
from numpy import __array_namespace_info__
import numpy._pytesttester
from numpy import all
from numpy import allclose
from numpy import amax
from numpy import amin
from numpy import angle
from numpy import any
from numpy import append
from numpy import apply_along_axis
from numpy import apply_over_axes
from numpy import arange
from numpy import argmax
from numpy import argmin
from numpy import argpartition
from numpy import argsort
from numpy import argwhere
from numpy import around
from numpy import array
from numpy import array2string
from numpy import array_equal
from numpy import array_equiv
from numpy import array_repr
from numpy import array_split
from numpy import array_str
from numpy import asanyarray
from numpy import asarray
from numpy import asarray_chkfinite
from numpy import ascontiguousarray
from numpy import asfortranarray
from numpy import asmatrix
from numpy import astype
from numpy import atleast_1d
from numpy import atleast_2d
from numpy import atleast_3d
from numpy import average
from numpy import bartlett
from numpy import base_repr
from numpy import binary_repr
from numpy import bincount
from numpy import blackman
from numpy import block
from numpy import bmat
from numpy import bool as bool_
from numpy import bool
from numpy import broadcast
from numpy import broadcast_arrays
from numpy import broadcast_shapes
from numpy import broadcast_to
from numpy import busday_count
from numpy import busday_offset
from numpy import busdaycalendar
from numpy import bytes_
from numpy import can_cast
from numpy import char
from numpy import character
from numpy import choose
from numpy import clip
from numpy import clongdouble
from numpy import column_stack
from numpy import common_type
from numpy import complex128 as cdouble
from numpy import complex128
from numpy import complex64
from numpy import complex64 as csingle
from numpy import complexfloating
from numpy import compress
from numpy import concatenate as concat
from numpy import concatenate
from numpy import convolve
from numpy import copy
from numpy import copyto
from numpy import core
from numpy import corrcoef
from numpy import correlate
from numpy import count_nonzero
from numpy import cov
from numpy import cross
from numpy import ctypeslib
from numpy import cumprod
from numpy import cumsum
from numpy import cumulative_prod
from numpy import cumulative_sum
from numpy import datetime64
from numpy import datetime_as_string
from numpy import datetime_data
from numpy import delete
from numpy import diag
from numpy import diag_indices
from numpy import diag_indices_from
from numpy import diagflat
from numpy import diagonal
from numpy import diff
from numpy import digitize
from numpy import dot
from numpy import dsplit
from numpy import dstack
from numpy import dtype
from numpy import dtypes
from numpy import ediff1d
from numpy import einsum
from numpy import einsum_path
from numpy import empty
from numpy import empty_like
from numpy import errstate
from numpy import exceptions
from numpy import expand_dims
from numpy import extract
from numpy import eye
from numpy import f2py
from numpy import fft
from numpy import fill_diagonal
from numpy import finfo
from numpy import fix
from numpy import flatiter
from numpy import flatnonzero
from numpy import flexible
from numpy import flip
from numpy import fliplr
from numpy import flipud
from numpy import float16
from numpy import float16 as half
from numpy import float32
from numpy import float32 as single
from numpy import float64 as double
from numpy import float64
from numpy import floating
from numpy import format_float_positional
from numpy import format_float_scientific
from numpy import from_dlpack
from numpy import frombuffer
from numpy import fromfile
from numpy import fromfunction
from numpy import fromiter
from numpy import frompyfunc
from numpy import fromregex
from numpy import fromstring
from numpy import full
from numpy import full_like
from numpy import generic
from numpy import genfromtxt
from numpy import geomspace
from numpy import get_include
from numpy import get_printoptions
from numpy import getbufsize
from numpy import geterr
from numpy import geterrcall
from numpy import gradient
from numpy import hamming
from numpy import hanning
from numpy import histogram
from numpy import histogram2d
from numpy import histogram_bin_edges
from numpy import histogramdd
from numpy import hsplit
from numpy import hstack
from numpy import i0
from numpy import identity
from numpy import iinfo
from numpy import imag
from numpy import indices
from numpy import inexact
from numpy import info
from numpy import inner
from numpy import insert
from numpy import int16 as short
from numpy import int16
from numpy import int32 as intc
from numpy import int32
from numpy import int64 as long
from numpy import int64
from numpy import int64 as int_
from numpy import int64 as intp
from numpy import int8
from numpy import int8 as byte
from numpy import integer
from numpy import interp
from numpy import intersect1d
from numpy import is_busday
from numpy import isclose
from numpy import iscomplex
from numpy import iscomplexobj
from numpy import isdtype
from numpy import isfortran
from numpy import isin
from numpy import isneginf
from numpy import isposinf
from numpy import isreal
from numpy import isrealobj
from numpy import isscalar
from numpy import issubdtype
from numpy import iterable
from numpy import ix_
from numpy import kaiser
from numpy import kron
from numpy import lexsort
from numpy import lib
import numpy.lib._index_tricks_impl
from numpy.lib import scimath as emath
from numpy import linalg
from numpy import linspace
from numpy import load
from numpy import loadtxt
from numpy import logspace
from numpy import longdouble
from numpy import longlong
from numpy import ma
from numpy import mask_indices
from numpy import matrix
from numpy import matrix_transpose
from numpy import max
from numpy import may_share_memory
from numpy import mean
from numpy import median
from numpy import memmap
from numpy import meshgrid
from numpy import min
from numpy import min_scalar_type
from numpy import mintypecode
from numpy import moveaxis
from numpy import nan_to_num
from numpy import nanargmax
from numpy import nanargmin
from numpy import nancumprod
from numpy import nancumsum
from numpy import nanmax
from numpy import nanmean
from numpy import nanmedian
from numpy import nanmin
from numpy import nanpercentile
from numpy import nanprod
from numpy import nanquantile
from numpy import nanstd
from numpy import nansum
from numpy import nanvar
from numpy import ndarray
from numpy import ndenumerate
from numpy import ndim
from numpy import ndindex
from numpy import nditer
from numpy import nested_iters
from numpy import nonzero
from numpy import number
from numpy import object_
from numpy import ones
from numpy import ones_like
from numpy import outer
from numpy import packbits
from numpy import pad
from numpy import partition
from numpy import percentile
from numpy import piecewise
from numpy import place
from numpy import poly
from numpy import poly1d
from numpy import polyadd
from numpy import polyder
from numpy import polydiv
from numpy import polyfit
from numpy import polyint
from numpy import polymul
from numpy import polynomial
from numpy import polysub
from numpy import polyval
from numpy import printoptions
from numpy import prod
from numpy import promote_types
from numpy import ptp
from numpy import put
from numpy import put_along_axis
from numpy import putmask
from numpy import quantile
from numpy import random
from numpy import ravel
from numpy import ravel_multi_index
from numpy import real
from numpy import real_if_close
from numpy import rec
from numpy.rec import recarray
from numpy import record
from numpy import repeat
from numpy import require
from numpy import reshape
from numpy import resize
from numpy import result_type
from numpy import roll
from numpy import rollaxis
from numpy import roots
from numpy import rot90
from numpy import round
from numpy import row_stack
from numpy import save
from numpy import savetxt
from numpy import savez
from numpy import savez_compressed
from numpy import searchsorted
from numpy import select
from numpy import set_printoptions
from numpy import setbufsize
from numpy import setdiff1d
from numpy import seterr
from numpy import seterrcall
from numpy import setxor1d
from numpy import shape
from numpy import shares_memory
from numpy import show_config
from numpy import show_runtime
from numpy import signedinteger
from numpy import sinc
from numpy import size
from numpy import sort
from numpy import sort_complex
from numpy import split
from numpy import squeeze
from numpy import stack
from numpy import std
from numpy import str_
from numpy import strings
from numpy import sum
from numpy import swapaxes
from numpy import take
from numpy import take_along_axis
from numpy import tensordot
from numpy import testing
from numpy import tile
from numpy import timedelta64
from numpy import trace
from numpy import transpose
from numpy import transpose as permute_dims
from numpy import trapezoid
from numpy import tri
from numpy import tril
from numpy import tril_indices
from numpy import tril_indices_from
from numpy import trim_zeros
from numpy import triu
from numpy import triu_indices
from numpy import triu_indices_from
from numpy import typename
from numpy import typing
from numpy import ufunc
from numpy import uint16
from numpy import uint16 as ushort
from numpy import uint32
from numpy import uint32 as uintc
from numpy import uint64
from numpy import uint64 as ulong
from numpy import uint64 as uint
from numpy import uint64 as uintp
from numpy import uint8 as ubyte
from numpy import uint8
from numpy import ulonglong
from numpy import union1d
from numpy import unique
from numpy import unique_all
from numpy import unique_counts
from numpy import unique_inverse
from numpy import unique_values
from numpy import unpackbits
from numpy import unravel_index
from numpy import unsignedinteger
from numpy import unstack
from numpy import unwrap
from numpy import vander
from numpy import var
from numpy import vdot
from numpy import vectorize
from numpy import void
from numpy import vsplit
from numpy import vstack
from numpy import where
from numpy import zeros
from numpy import zeros_like
from rdkit.sping import pid
__all__: list[str] = ['DrawSpiral', 'False_', 'ScalarType', 'True_', 'abs', 'absolute', 'acos', 'acosh', 'add', 'all', 'allclose', 'amax', 'amin', 'angle', 'any', 'append', 'apply_along_axis', 'apply_over_axes', 'arange', 'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan', 'arctan2', 'arctanh', 'argmax', 'argmin', 'argpartition', 'argsort', 'argwhere', 'around', 'array', 'array2string', 'array_equal', 'array_equiv', 'array_repr', 'array_split', 'array_str', 'asanyarray', 'asarray', 'asarray_chkfinite', 'ascontiguousarray', 'asfortranarray', 'asin', 'asinh', 'asmatrix', 'astype', 'atan', 'atan2', 'atanh', 'atleast_1d', 'atleast_2d', 'atleast_3d', 'average', 'bartlett', 'base_repr', 'binary_repr', 'bincount', 'bitwise_and', 'bitwise_count', 'bitwise_invert', 'bitwise_left_shift', 'bitwise_not', 'bitwise_or', 'bitwise_right_shift', 'bitwise_xor', 'blackman', 'block', 'bmat', 'bool', 'bool_', 'broadcast', 'broadcast_arrays', 'broadcast_shapes', 'broadcast_to', 'busday_count', 'busday_offset', 'busdaycalendar', 'byte', 'bytes_', 'c_', 'can_cast', 'cbrt', 'cdouble', 'ceil', 'char', 'character', 'choose', 'clip', 'clongdouble', 'column_stack', 'common_type', 'complex128', 'complex64', 'complexfloating', 'compress', 'concat', 'concatenate', 'conj', 'conjugate', 'convolve', 'copy', 'copysign', 'copyto', 'core', 'corrcoef', 'correlate', 'cos', 'cosh', 'count_nonzero', 'cov', 'cross', 'csingle', 'ctypeslib', 'cumprod', 'cumsum', 'cumulative_prod', 'cumulative_sum', 'datetime64', 'datetime_as_string', 'datetime_data', 'deg2rad', 'degrees', 'delete', 'diag', 'diag_indices', 'diag_indices_from', 'diagflat', 'diagonal', 'diff', 'digitize', 'divide', 'divmod', 'dot', 'double', 'dsplit', 'dstack', 'dtype', 'dtypes', 'e', 'ediff1d', 'einsum', 'einsum_path', 'emath', 'empty', 'empty_like', 'equal', 'errstate', 'euler_gamma', 'exceptions', 'exp', 'exp2', 'expand_dims', 'expm1', 'extract', 'eye', 'f2py', 'fabs', 'fft', 'fill_diagonal', 'finfo', 'fix', 'flatiter', 'flatnonzero', 'flexible', 'flip', 'fliplr', 'flipud', 'float16', 'float32', 'float64', 'float_power', 'floating', 'floor', 'floor_divide', 'fmax', 'fmin', 'fmod', 'format_float_positional', 'format_float_scientific', 'frexp', 'from_dlpack', 'frombuffer', 'fromfile', 'fromfunction', 'fromiter', 'frompyfunc', 'fromregex', 'fromstring', 'full', 'full_like', 'gcd', 'generic', 'genfromtxt', 'geomspace', 'get_include', 'get_printoptions', 'getbufsize', 'geterr', 'geterrcall', 'gradient', 'greater', 'greater_equal', 'half', 'hamming', 'hanning', 'heaviside', 'histogram', 'histogram2d', 'histogram_bin_edges', 'histogramdd', 'hsplit', 'hstack', 'hypot', 'i0', 'identity', 'iinfo', 'imag', 'index_exp', 'indices', 'inexact', 'inf', 'info', 'inner', 'insert', 'int16', 'int32', 'int64', 'int8', 'int_', 'intc', 'integer', 'interp', 'intersect1d', 'intp', 'invert', 'is_busday', 'isclose', 'iscomplex', 'iscomplexobj', 'isdtype', 'isfinite', 'isfortran', 'isin', 'isinf', 'isnan', 'isnat', 'isneginf', 'isposinf', 'isreal', 'isrealobj', 'isscalar', 'issubdtype', 'iterable', 'ix_', 'kaiser', 'kron', 'lcm', 'ldexp', 'left_shift', 'less', 'less_equal', 'lexsort', 'lib', 'linalg', 'linspace', 'little_endian', 'load', 'loadtxt', 'log', 'log10', 'log1p', 'log2', 'logaddexp', 'logaddexp2', 'logical_and', 'logical_not', 'logical_or', 'logical_xor', 'logspace', 'long', 'longdouble', 'longlong', 'ma', 'mask_indices', 'math', 'matmul', 'matrix', 'matrix_transpose', 'matvec', 'max', 'maximum', 'may_share_memory', 'mean', 'median', 'memmap', 'meshgrid', 'mgrid', 'min', 'min_scalar_type', 'minimum', 'mintypecode', 'mod', 'modf', 'moveaxis', 'multiply', 'nan', 'nan_to_num', 'nanargmax', 'nanargmin', 'nancumprod', 'nancumsum', 'nanmax', 'nanmean', 'nanmedian', 'nanmin', 'nanpercentile', 'nanprod', 'nanquantile', 'nanstd', 'nansum', 'nanvar', 'ndarray', 'ndenumerate', 'ndim', 'ndindex', 'nditer', 'negative', 'nested_iters', 'newaxis', 'nextafter', 'nonzero', 'not_equal', 'number', 'object_', 'ogrid', 'ones', 'ones_like', 'outer', 'packbits', 'pad', 'partition', 'percentile', 'permute_dims', 'pi', 'pid', 'piecewise', 'place', 'poly', 'poly1d', 'polyadd', 'polyder', 'polydiv', 'polyfit', 'polyint', 'polymul', 'polynomial', 'polysub', 'polyval', 'positive', 'pow', 'power', 'printoptions', 'prod', 'promote_types', 'ptp', 'put', 'put_along_axis', 'putmask', 'quantile', 'r_', 'rad2deg', 'radians', 'random', 'ravel', 'ravel_multi_index', 'real', 'real_if_close', 'rec', 'recarray', 'reciprocal', 'record', 'remainder', 'repeat', 'require', 'reshape', 'resize', 'result_type', 'right_shift', 'rint', 'roll', 'rollaxis', 'roots', 'rot90', 'round', 'row_stack', 's_', 'save', 'savetxt', 'savez', 'savez_compressed', 'sctypeDict', 'searchsorted', 'select', 'set_printoptions', 'setbufsize', 'setdiff1d', 'seterr', 'seterrcall', 'setxor1d', 'shape', 'shares_memory', 'short', 'show_config', 'show_runtime', 'sign', 'signbit', 'signedinteger', 'sin', 'sinc', 'single', 'sinh', 'size', 'sort', 'sort_complex', 'spacing', 'split', 'sqrt', 'square', 'squeeze', 'stack', 'std', 'str_', 'strings', 'subtract', 'sum', 'swapaxes', 'take', 'take_along_axis', 'tan', 'tanh', 'tensordot', 'test', 'testing', 'tile', 'timedelta64', 'trace', 'transpose', 'trapezoid', 'tri', 'tril', 'tril_indices', 'tril_indices_from', 'trim_zeros', 'triu', 'triu_indices', 'triu_indices_from', 'true_divide', 'trunc', 'typecodes', 'typename', 'typing', 'ubyte', 'ufunc', 'uint', 'uint16', 'uint32', 'uint64', 'uint8', 'uintc', 'uintp', 'ulong', 'ulonglong', 'union1d', 'unique', 'unique_all', 'unique_counts', 'unique_inverse', 'unique_values', 'unpackbits', 'unravel_index', 'unsignedinteger', 'unstack', 'unwrap', 'ushort', 'vander', 'var', 'vdot', 'vecdot', 'vecmat', 'vectorize', 'void', 'vsplit', 'vstack', 'where', 'zeros', 'zeros_like']
def DrawSpiral(canvas, startColor, endColor, startRadius, endRadius, nLoops, degsPerSlice = 70, degsPerStep = 1, startAngle = 0, centerPos = None, dir = 1):
    ...
False_: numpy.bool  # value = np.False_
ScalarType: tuple = (int, float, complex, bool, bytes, str, memoryview, numpy.bool, numpy.complex64, numpy.complex128, numpy.clongdouble, numpy.float16, numpy.float32, numpy.float64, numpy.longdouble, numpy.int8, numpy.int16, numpy.int32, numpy.int64, numpy.longlong, numpy.datetime64, numpy.timedelta64, numpy.object_, numpy.bytes_, numpy.str_, numpy.uint8, numpy.uint16, numpy.uint32, numpy.uint64, numpy.ulonglong, numpy.void)
True_: numpy.bool  # value = np.True_
__version__: str = '2.4.6'
abs: numpy.ufunc  # value = <ufunc 'absolute'>
r"""
absolute(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Calculate the absolute value element-wise.

``np.abs`` is a shorthand for this function.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
absolute : ndarray
    An ndarray containing the absolute value of
    each element in `x`.  For complex input, ``a + ib``, the
    absolute value is :math:`\sqrt{ a^2 + b^2 }`.
    This is a scalar if `x` is a scalar.

Examples
--------
>>> import numpy as np
>>> x = np.array([-1.2, 1.2])
>>> np.absolute(x)
array([ 1.2,  1.2])
>>> np.absolute(1.2 + 1j)
1.5620499351813308

Plot the function over ``[-10, 10]``:

>>> import matplotlib.pyplot as plt

>>> x = np.linspace(start=-10, stop=10, num=101)
>>> plt.plot(x, np.absolute(x))
>>> plt.show()

Plot the function over the complex plane:

>>> xx = x + 1j * x[:, np.newaxis]
>>> plt.imshow(np.abs(xx), extent=[-10, 10, -10, 10], cmap='gray')
>>> plt.show()

The `abs` function can be used as a shorthand for ``np.absolute`` on
ndarrays.

>>> x = np.array([-1.2, 1.2])
>>> abs(x)
array([1.2, 1.2])
"""
absolute: numpy.ufunc  # value = <ufunc 'absolute'>
r"""
absolute(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Calculate the absolute value element-wise.

``np.abs`` is a shorthand for this function.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
absolute : ndarray
    An ndarray containing the absolute value of
    each element in `x`.  For complex input, ``a + ib``, the
    absolute value is :math:`\sqrt{ a^2 + b^2 }`.
    This is a scalar if `x` is a scalar.

Examples
--------
>>> import numpy as np
>>> x = np.array([-1.2, 1.2])
>>> np.absolute(x)
array([ 1.2,  1.2])
>>> np.absolute(1.2 + 1j)
1.5620499351813308

Plot the function over ``[-10, 10]``:

>>> import matplotlib.pyplot as plt

>>> x = np.linspace(start=-10, stop=10, num=101)
>>> plt.plot(x, np.absolute(x))
>>> plt.show()

Plot the function over the complex plane:

>>> xx = x + 1j * x[:, np.newaxis]
>>> plt.imshow(np.abs(xx), extent=[-10, 10, -10, 10], cmap='gray')
>>> plt.show()

The `abs` function can be used as a shorthand for ``np.absolute`` on
ndarrays.

>>> x = np.array([-1.2, 1.2])
>>> abs(x)
array([1.2, 1.2])
"""
acos: numpy.ufunc  # value = <ufunc 'arccos'>
"""
arccos(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Trigonometric inverse cosine, element-wise.

The inverse of `cos` so that, if ``y = cos(x)``, then ``x = arccos(y)``.

Parameters
----------
x : array_like
    `x`-coordinate on the unit circle.
    For real arguments, the domain is [-1, 1].
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    The angle of the ray intersecting the unit circle at the given
    `x`-coordinate in radians [0, pi].
    This is a scalar if `x` is a scalar.

See Also
--------
cos, arctan, arcsin, emath.arccos

Notes
-----
`arccos` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that ``cos(z) = x``. The convention is to return
the angle `z` whose real part lies in `[0, pi]`.

For real-valued input data types, `arccos` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arccos` is a complex analytic function that
has branch cuts ``[-inf, -1]`` and `[1, inf]` and is continuous from
above on the former and from below on the latter.

The inverse `cos` is also known as `acos` or cos^-1.

References
----------
M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
10th printing, 1964, pp. 79.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------
>>> import numpy as np

We expect the arccos of 1 to be 0, and of -1 to be pi:

>>> np.arccos([1, -1])
array([ 0.        ,  3.14159265])

Plot arccos:

>>> import matplotlib.pyplot as plt
>>> x = np.linspace(-1, 1, num=100)
>>> plt.plot(x, np.arccos(x))
>>> plt.axis('tight')
>>> plt.show()
"""
acosh: numpy.ufunc  # value = <ufunc 'arccosh'>
"""
arccosh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic cosine, element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
arccosh : ndarray
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------

cosh, arcsinh, sinh, arctanh, tanh

Notes
-----
`arccosh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that `cosh(z) = x`. The convention is to return the
`z` whose imaginary part lies in ``[-pi, pi]`` and the real part in
``[0, inf]``.

For real-valued input data types, `arccosh` always returns real output.
For each value that cannot be expressed as a real number or infinity, it
yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arccosh` is a complex analytical function that
has a branch cut `[-inf, 1]` and is continuous from above on it.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arccosh

Examples
--------
>>> import numpy as np
>>> np.arccosh([np.e, 10.0])
array([ 1.65745445,  2.99322285])
>>> np.arccosh(1)
0.0
"""
add: numpy.ufunc  # value = <ufunc 'add'>
"""
add(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Add arguments element-wise.

Parameters
----------
x1, x2 : array_like
    The arrays to be added.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
add : ndarray or scalar
    The sum of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

Notes
-----
Equivalent to `x1` + `x2` in terms of array broadcasting.

Examples
--------
>>> import numpy as np
>>> np.add(1.0, 4.0)
5.0
>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> np.add(x1, x2)
array([[  0.,   2.,   4.],
       [  3.,   5.,   7.],
       [  6.,   8.,  10.]])

The ``+`` operator can be used as a shorthand for ``np.add`` on ndarrays.

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> x1 + x2
array([[ 0.,  2.,  4.],
       [ 3.,  5.,  7.],
       [ 6.,  8., 10.]])
"""
arccos: numpy.ufunc  # value = <ufunc 'arccos'>
"""
arccos(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Trigonometric inverse cosine, element-wise.

The inverse of `cos` so that, if ``y = cos(x)``, then ``x = arccos(y)``.

Parameters
----------
x : array_like
    `x`-coordinate on the unit circle.
    For real arguments, the domain is [-1, 1].
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    The angle of the ray intersecting the unit circle at the given
    `x`-coordinate in radians [0, pi].
    This is a scalar if `x` is a scalar.

See Also
--------
cos, arctan, arcsin, emath.arccos

Notes
-----
`arccos` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that ``cos(z) = x``. The convention is to return
the angle `z` whose real part lies in `[0, pi]`.

For real-valued input data types, `arccos` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arccos` is a complex analytic function that
has branch cuts ``[-inf, -1]`` and `[1, inf]` and is continuous from
above on the former and from below on the latter.

The inverse `cos` is also known as `acos` or cos^-1.

References
----------
M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
10th printing, 1964, pp. 79.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------
>>> import numpy as np

We expect the arccos of 1 to be 0, and of -1 to be pi:

>>> np.arccos([1, -1])
array([ 0.        ,  3.14159265])

Plot arccos:

>>> import matplotlib.pyplot as plt
>>> x = np.linspace(-1, 1, num=100)
>>> plt.plot(x, np.arccos(x))
>>> plt.axis('tight')
>>> plt.show()
"""
arccosh: numpy.ufunc  # value = <ufunc 'arccosh'>
"""
arccosh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic cosine, element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
arccosh : ndarray
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------

cosh, arcsinh, sinh, arctanh, tanh

Notes
-----
`arccosh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that `cosh(z) = x`. The convention is to return the
`z` whose imaginary part lies in ``[-pi, pi]`` and the real part in
``[0, inf]``.

For real-valued input data types, `arccosh` always returns real output.
For each value that cannot be expressed as a real number or infinity, it
yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arccosh` is a complex analytical function that
has a branch cut `[-inf, 1]` and is continuous from above on it.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arccosh

Examples
--------
>>> import numpy as np
>>> np.arccosh([np.e, 10.0])
array([ 1.65745445,  2.99322285])
>>> np.arccosh(1)
0.0
"""
arcsin: numpy.ufunc  # value = <ufunc 'arcsin'>
"""
arcsin(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse sine, element-wise.

Parameters
----------
x : array_like
    `y`-coordinate on the unit circle.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    The inverse sine of each element in `x`, in radians and in the
    closed interval ``[-pi/2, pi/2]``.
    This is a scalar if `x` is a scalar.

See Also
--------
sin, cos, arccos, tan, arctan, arctan2, emath.arcsin

Notes
-----
`arcsin` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that :math:`sin(z) = x`.  The convention is to
return the angle `z` whose real part lies in [-pi/2, pi/2].

For real-valued input data types, *arcsin* always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arcsin` is a complex analytic function that
has, by convention, the branch cuts [-inf, -1] and [1, inf]  and is
continuous from above on the former and from below on the latter.

The inverse sine is also known as `asin` or sin^{-1}.

References
----------
Abramowitz, M. and Stegun, I. A., *Handbook of Mathematical Functions*,
10th printing, New York: Dover, 1964, pp. 79ff.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------
>>> import numpy as np
>>> np.arcsin(1)     # pi/2
1.5707963267948966
>>> np.arcsin(-1)    # -pi/2
-1.5707963267948966
>>> np.arcsin(0)
0.0
"""
arcsinh: numpy.ufunc  # value = <ufunc 'arcsinh'>
"""
arcsinh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic sine element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

Notes
-----
`arcsinh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that `sinh(z) = x`. The convention is to return the
`z` whose imaginary part lies in `[-pi/2, pi/2]`.

For real-valued input data types, `arcsinh` always returns real output.
For each value that cannot be expressed as a real number or infinity, it
returns ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arcsinh` is a complex analytical function that
has branch cuts `[1j, infj]` and `[-1j, -infj]` and is continuous from
the right on the former and from the left on the latter.

The inverse hyperbolic sine is also known as `asinh` or ``sinh^-1``.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arcsinh

Examples
--------
>>> import numpy as np
>>> np.arcsinh(np.array([np.e, 10.0]))
array([ 1.72538256,  2.99822295])
"""
arctan: numpy.ufunc  # value = <ufunc 'arctan'>
"""
arctan(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Trigonometric inverse tangent, element-wise.

The inverse of tan, so that if ``y = tan(x)`` then ``x = arctan(y)``.

Parameters
----------
x : array_like
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Out has the same shape as `x`.  Its real part is in
    ``[-pi/2, pi/2]`` (``arctan(+/-inf)`` returns ``+/-pi/2``).
    This is a scalar if `x` is a scalar.

See Also
--------
arctan2 : The "four quadrant" arctan of the angle formed by (`x`, `y`)
    and the positive `x`-axis.
angle : Argument of complex values.

Notes
-----
`arctan` is a multi-valued function: for each `x` there are infinitely
many numbers `z` such that tan(`z`) = `x`.  The convention is to return
the angle `z` whose real part lies in [-pi/2, pi/2].

For real-valued input data types, `arctan` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arctan` is a complex analytic function that
has [``1j, infj``] and [``-1j, -infj``] as branch cuts, and is continuous
from the left on the former and from the right on the latter.

The inverse tangent is also known as `atan` or tan^{-1}.

References
----------
Abramowitz, M. and Stegun, I. A., *Handbook of Mathematical Functions*,
10th printing, New York: Dover, 1964, pp. 79.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------

We expect the arctan of 0 to be 0, and of 1 to be pi/4:

>>> import numpy as np
>>> np.arctan([0, 1])
array([ 0.        ,  0.78539816])

>>> np.pi/4
0.78539816339744828

Plot arctan:

>>> import matplotlib.pyplot as plt
>>> x = np.linspace(-10, 10)
>>> plt.plot(x, np.arctan(x))
>>> plt.axis('tight')
>>> plt.show()
"""
arctan2: numpy.ufunc  # value = <ufunc 'arctan2'>
"""
arctan2(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise arc tangent of ``x1/x2`` choosing the quadrant correctly.

The quadrant (i.e., branch) is chosen so that ``arctan2(x1, x2)`` is
the signed angle in radians between the ray ending at the origin and
passing through the point (1,0), and the ray ending at the origin and
passing through the point (`x2`, `x1`).  (Note the role reversal: the
"`y`-coordinate" is the first function parameter, the "`x`-coordinate"
is the second.)  By IEEE convention, this function is defined for
`x2` = +/-0 and for either or both of `x1` and `x2` = +/-inf (see
Notes for specific values).

This function is not defined for complex-valued arguments; for the
so-called argument of complex values, use `angle`.

Parameters
----------
x1 : array_like, real-valued
    `y`-coordinates.
x2 : array_like, real-valued
    `x`-coordinates.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    Array of angles in radians, in the range ``[-pi, pi]``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
arctan, tan, angle

Notes
-----
*arctan2* is identical to the `atan2` function of the underlying
C library.  The following special values are defined in the C
standard: [1]_

====== ====== ================
`x1`   `x2`   `arctan2(x1,x2)`
====== ====== ================
+/- 0  +0     +/- 0
+/- 0  -0     +/- pi
 > 0   +/-inf +0 / +pi
 < 0   +/-inf -0 / -pi
+/-inf +inf   +/- (pi/4)
+/-inf -inf   +/- (3*pi/4)
====== ====== ================

Note that +0 and -0 are distinct floating point numbers, as are +inf
and -inf.

References
----------
.. [1] ISO/IEC standard 9899:1999, "Programming language C."

Examples
--------

Consider four points in different quadrants:

>>> import numpy as np
>>> x = np.array([-1, +1, +1, -1])
>>> y = np.array([-1, -1, +1, +1])
>>> np.arctan2(y, x) * 180 / np.pi
array([-135.,  -45.,   45.,  135.])

Note the order of the parameters. `arctan2` is defined also when `x2` = 0
and at several other special points, obtaining values in
the range ``[-pi, pi]``:

>>> np.arctan2([1., -1.], [0., 0.])
array([ 1.57079633, -1.57079633])
>>> np.arctan2([0., 0., np.inf], [+0., -0., np.inf])
array([0.        , 3.14159265, 0.78539816])
"""
arctanh: numpy.ufunc  # value = <ufunc 'arctanh'>
"""
arctanh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic tangent element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
emath.arctanh

Notes
-----
`arctanh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that ``tanh(z) = x``. The convention is to return
the `z` whose imaginary part lies in `[-pi/2, pi/2]`.

For real-valued input data types, `arctanh` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arctanh` is a complex analytical function
that has branch cuts `[-1, -inf]` and `[1, inf]` and is continuous from
above on the former and from below on the latter.

The inverse hyperbolic tangent is also known as `atanh` or ``tanh^-1``.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arctanh

Examples
--------
>>> import numpy as np
>>> np.arctanh([0, -0.5])
array([ 0.        , -0.54930614])
"""
asin: numpy.ufunc  # value = <ufunc 'arcsin'>
"""
arcsin(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse sine, element-wise.

Parameters
----------
x : array_like
    `y`-coordinate on the unit circle.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    The inverse sine of each element in `x`, in radians and in the
    closed interval ``[-pi/2, pi/2]``.
    This is a scalar if `x` is a scalar.

See Also
--------
sin, cos, arccos, tan, arctan, arctan2, emath.arcsin

Notes
-----
`arcsin` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that :math:`sin(z) = x`.  The convention is to
return the angle `z` whose real part lies in [-pi/2, pi/2].

For real-valued input data types, *arcsin* always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arcsin` is a complex analytic function that
has, by convention, the branch cuts [-inf, -1] and [1, inf]  and is
continuous from above on the former and from below on the latter.

The inverse sine is also known as `asin` or sin^{-1}.

References
----------
Abramowitz, M. and Stegun, I. A., *Handbook of Mathematical Functions*,
10th printing, New York: Dover, 1964, pp. 79ff.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------
>>> import numpy as np
>>> np.arcsin(1)     # pi/2
1.5707963267948966
>>> np.arcsin(-1)    # -pi/2
-1.5707963267948966
>>> np.arcsin(0)
0.0
"""
asinh: numpy.ufunc  # value = <ufunc 'arcsinh'>
"""
arcsinh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic sine element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

Notes
-----
`arcsinh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that `sinh(z) = x`. The convention is to return the
`z` whose imaginary part lies in `[-pi/2, pi/2]`.

For real-valued input data types, `arcsinh` always returns real output.
For each value that cannot be expressed as a real number or infinity, it
returns ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arcsinh` is a complex analytical function that
has branch cuts `[1j, infj]` and `[-1j, -infj]` and is continuous from
the right on the former and from the left on the latter.

The inverse hyperbolic sine is also known as `asinh` or ``sinh^-1``.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arcsinh

Examples
--------
>>> import numpy as np
>>> np.arcsinh(np.array([np.e, 10.0]))
array([ 1.72538256,  2.99822295])
"""
atan: numpy.ufunc  # value = <ufunc 'arctan'>
"""
arctan(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Trigonometric inverse tangent, element-wise.

The inverse of tan, so that if ``y = tan(x)`` then ``x = arctan(y)``.

Parameters
----------
x : array_like
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Out has the same shape as `x`.  Its real part is in
    ``[-pi/2, pi/2]`` (``arctan(+/-inf)`` returns ``+/-pi/2``).
    This is a scalar if `x` is a scalar.

See Also
--------
arctan2 : The "four quadrant" arctan of the angle formed by (`x`, `y`)
    and the positive `x`-axis.
angle : Argument of complex values.

Notes
-----
`arctan` is a multi-valued function: for each `x` there are infinitely
many numbers `z` such that tan(`z`) = `x`.  The convention is to return
the angle `z` whose real part lies in [-pi/2, pi/2].

For real-valued input data types, `arctan` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arctan` is a complex analytic function that
has [``1j, infj``] and [``-1j, -infj``] as branch cuts, and is continuous
from the left on the former and from the right on the latter.

The inverse tangent is also known as `atan` or tan^{-1}.

References
----------
Abramowitz, M. and Stegun, I. A., *Handbook of Mathematical Functions*,
10th printing, New York: Dover, 1964, pp. 79.
https://personal.math.ubc.ca/~cbm/aands/page_79.htm

Examples
--------

We expect the arctan of 0 to be 0, and of 1 to be pi/4:

>>> import numpy as np
>>> np.arctan([0, 1])
array([ 0.        ,  0.78539816])

>>> np.pi/4
0.78539816339744828

Plot arctan:

>>> import matplotlib.pyplot as plt
>>> x = np.linspace(-10, 10)
>>> plt.plot(x, np.arctan(x))
>>> plt.axis('tight')
>>> plt.show()
"""
atan2: numpy.ufunc  # value = <ufunc 'arctan2'>
"""
arctan2(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise arc tangent of ``x1/x2`` choosing the quadrant correctly.

The quadrant (i.e., branch) is chosen so that ``arctan2(x1, x2)`` is
the signed angle in radians between the ray ending at the origin and
passing through the point (1,0), and the ray ending at the origin and
passing through the point (`x2`, `x1`).  (Note the role reversal: the
"`y`-coordinate" is the first function parameter, the "`x`-coordinate"
is the second.)  By IEEE convention, this function is defined for
`x2` = +/-0 and for either or both of `x1` and `x2` = +/-inf (see
Notes for specific values).

This function is not defined for complex-valued arguments; for the
so-called argument of complex values, use `angle`.

Parameters
----------
x1 : array_like, real-valued
    `y`-coordinates.
x2 : array_like, real-valued
    `x`-coordinates.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
angle : ndarray
    Array of angles in radians, in the range ``[-pi, pi]``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
arctan, tan, angle

Notes
-----
*arctan2* is identical to the `atan2` function of the underlying
C library.  The following special values are defined in the C
standard: [1]_

====== ====== ================
`x1`   `x2`   `arctan2(x1,x2)`
====== ====== ================
+/- 0  +0     +/- 0
+/- 0  -0     +/- pi
 > 0   +/-inf +0 / +pi
 < 0   +/-inf -0 / -pi
+/-inf +inf   +/- (pi/4)
+/-inf -inf   +/- (3*pi/4)
====== ====== ================

Note that +0 and -0 are distinct floating point numbers, as are +inf
and -inf.

References
----------
.. [1] ISO/IEC standard 9899:1999, "Programming language C."

Examples
--------

Consider four points in different quadrants:

>>> import numpy as np
>>> x = np.array([-1, +1, +1, -1])
>>> y = np.array([-1, -1, +1, +1])
>>> np.arctan2(y, x) * 180 / np.pi
array([-135.,  -45.,   45.,  135.])

Note the order of the parameters. `arctan2` is defined also when `x2` = 0
and at several other special points, obtaining values in
the range ``[-pi, pi]``:

>>> np.arctan2([1., -1.], [0., 0.])
array([ 1.57079633, -1.57079633])
>>> np.arctan2([0., 0., np.inf], [+0., -0., np.inf])
array([0.        , 3.14159265, 0.78539816])
"""
atanh: numpy.ufunc  # value = <ufunc 'arctanh'>
"""
arctanh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Inverse hyperbolic tangent element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Array of the same shape as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
emath.arctanh

Notes
-----
`arctanh` is a multivalued function: for each `x` there are infinitely
many numbers `z` such that ``tanh(z) = x``. The convention is to return
the `z` whose imaginary part lies in `[-pi/2, pi/2]`.

For real-valued input data types, `arctanh` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `arctanh` is a complex analytical function
that has branch cuts `[-1, -inf]` and `[1, inf]` and is continuous from
above on the former and from below on the latter.

The inverse hyperbolic tangent is also known as `atanh` or ``tanh^-1``.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 86.
       https://personal.math.ubc.ca/~cbm/aands/page_86.htm
.. [2] Wikipedia, "Inverse hyperbolic function",
       https://en.wikipedia.org/wiki/Arctanh

Examples
--------
>>> import numpy as np
>>> np.arctanh([0, -0.5])
array([ 0.        , -0.54930614])
"""
bitwise_and: numpy.ufunc  # value = <ufunc 'bitwise_and'>
"""
bitwise_and(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the bit-wise AND of two arrays element-wise.

Computes the bit-wise AND of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``&``.

Parameters
----------
x1, x2 : array_like
    Only integer and boolean types are handled.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_and
bitwise_or
bitwise_xor
binary_repr :
    Return the binary representation of the input number as a string.

Examples
--------
>>> import numpy as np

The number 13 is represented by ``00001101``.  Likewise, 17 is
represented by ``00010001``.  The bit-wise AND of 13 and 17 is
therefore ``000000001``, or 1:

>>> np.bitwise_and(13, 17)
1

>>> np.bitwise_and(14, 13)
12
>>> np.binary_repr(12)
'1100'
>>> np.bitwise_and([14,3], 13)
array([12,  1])

>>> np.bitwise_and([11,7], [4,25])
array([0, 1])
>>> np.bitwise_and(np.array([2,5,255]), np.array([3,14,16]))
array([ 2,  4, 16])
>>> np.bitwise_and([True, True], [False, True])
array([False,  True])

The ``&`` operator can be used as a shorthand for ``np.bitwise_and`` on
ndarrays.

>>> x1 = np.array([2, 5, 255])
>>> x2 = np.array([3, 14, 16])
>>> x1 & x2
array([ 2,  4, 16])
"""
bitwise_count: numpy.ufunc  # value = <ufunc 'bitwise_count'>
"""
bitwise_count(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Computes the number of 1-bits in the absolute value of ``x``.
Analogous to the builtin `int.bit_count` or ``popcount`` in C++.

Parameters
----------
x : array_like, unsigned int
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding number of 1-bits in the input.
    Returns uint8 for all integer types
    This is a scalar if `x` is a scalar.

References
----------
.. [1] https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel

.. [2] Wikipedia, "Hamming weight",
       https://en.wikipedia.org/wiki/Hamming_weight

.. [3] http://aggregate.ee.engr.uky.edu/MAGIC/#Population%20Count%20(Ones%20Count)

Examples
--------
>>> import numpy as np
>>> np.bitwise_count(1023)
np.uint8(10)
>>> a = np.array([2**i - 1 for i in range(16)])
>>> np.bitwise_count(a)
array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15],
      dtype=uint8)
"""
bitwise_invert: numpy.ufunc  # value = <ufunc 'invert'>
"""
invert(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute bit-wise inversion, or bit-wise NOT, element-wise.

Computes the bit-wise NOT of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``~``.

For signed integer inputs, the bit-wise NOT of the absolute value is
returned. In a two's-complement system, this operation effectively flips
all the bits, resulting in a representation that corresponds to the
negative of the input plus one. This is the most common method of
representing signed integers on computers [1]_. A N-bit two's-complement
system can represent every integer in the range :math:`-2^{N-1}` to
:math:`+2^{N-1}-1`.

Parameters
----------
x : array_like
    Only integer and boolean types are handled.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if `x` is a scalar.

See Also
--------
bitwise_and, bitwise_or, bitwise_xor
logical_not
binary_repr :
    Return the binary representation of the input number as a string.

Notes
-----
``numpy.bitwise_not`` is an alias for `invert`:

>>> np.bitwise_not is np.invert
True

References
----------
.. [1] Wikipedia, "Two's complement",
    https://en.wikipedia.org/wiki/Two's_complement

Examples
--------
>>> import numpy as np

We've seen that 13 is represented by ``00001101``.
The invert or bit-wise NOT of 13 is then:

>>> x = np.invert(np.array(13, dtype=np.uint8))
>>> x
np.uint8(242)
>>> np.binary_repr(x, width=8)
'11110010'

The result depends on the bit-width:

>>> x = np.invert(np.array(13, dtype=np.uint16))
>>> x
np.uint16(65522)
>>> np.binary_repr(x, width=16)
'1111111111110010'

When using signed integer types, the result is the bit-wise NOT of
the unsigned type, interpreted as a signed integer:

>>> np.invert(np.array([13], dtype=np.int8))
array([-14], dtype=int8)
>>> np.binary_repr(-14, width=8)
'11110010'

Booleans are accepted as well:

>>> np.invert(np.array([True, False]))
array([False,  True])

The ``~`` operator can be used as a shorthand for ``np.invert`` on
ndarrays.

>>> x1 = np.array([True, False])
>>> ~x1
array([False,  True])
"""
bitwise_left_shift: numpy.ufunc  # value = <ufunc 'left_shift'>
"""
left_shift(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Shift the bits of an integer to the left.

Bits are shifted to the left by appending `x2` 0s at the right of `x1`.
Since the internal representation of numbers is in binary format, this
operation is equivalent to multiplying `x1` by ``2**x2``.

Parameters
----------
x1 : array_like of integer type
    Input values.
x2 : array_like of integer type
    Number of zeros to append to `x1`. Has to be non-negative.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : array of integer type
    Return `x1` with bits shifted `x2` times to the left.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
right_shift : Shift the bits of an integer to the right.
binary_repr : Return the binary representation of the input number
    as a string.

Examples
--------
>>> import numpy as np
>>> np.binary_repr(5)
'101'
>>> np.left_shift(5, 2)
20
>>> np.binary_repr(20)
'10100'

>>> np.left_shift(5, [1,2,3])
array([10, 20, 40])

Note that the dtype of the second argument may change the dtype of the
result and can lead to unexpected results in some cases (see
:ref:`Casting Rules <ufuncs.casting>`):

>>> a = np.left_shift(np.uint8(255), np.int64(1))  # Expect 254
>>> print(a, type(a)) # Unexpected result due to upcasting
510 <class 'numpy.int64'>
>>> b = np.left_shift(np.uint8(255), np.uint8(1))
>>> print(b, type(b))
254 <class 'numpy.uint8'>

The ``<<`` operator can be used as a shorthand for ``np.left_shift`` on
ndarrays.

>>> x1 = 5
>>> x2 = np.array([1, 2, 3])
>>> x1 << x2
array([10, 20, 40])
"""
bitwise_not: numpy.ufunc  # value = <ufunc 'invert'>
"""
invert(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute bit-wise inversion, or bit-wise NOT, element-wise.

Computes the bit-wise NOT of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``~``.

For signed integer inputs, the bit-wise NOT of the absolute value is
returned. In a two's-complement system, this operation effectively flips
all the bits, resulting in a representation that corresponds to the
negative of the input plus one. This is the most common method of
representing signed integers on computers [1]_. A N-bit two's-complement
system can represent every integer in the range :math:`-2^{N-1}` to
:math:`+2^{N-1}-1`.

Parameters
----------
x : array_like
    Only integer and boolean types are handled.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if `x` is a scalar.

See Also
--------
bitwise_and, bitwise_or, bitwise_xor
logical_not
binary_repr :
    Return the binary representation of the input number as a string.

Notes
-----
``numpy.bitwise_not`` is an alias for `invert`:

>>> np.bitwise_not is np.invert
True

References
----------
.. [1] Wikipedia, "Two's complement",
    https://en.wikipedia.org/wiki/Two's_complement

Examples
--------
>>> import numpy as np

We've seen that 13 is represented by ``00001101``.
The invert or bit-wise NOT of 13 is then:

>>> x = np.invert(np.array(13, dtype=np.uint8))
>>> x
np.uint8(242)
>>> np.binary_repr(x, width=8)
'11110010'

The result depends on the bit-width:

>>> x = np.invert(np.array(13, dtype=np.uint16))
>>> x
np.uint16(65522)
>>> np.binary_repr(x, width=16)
'1111111111110010'

When using signed integer types, the result is the bit-wise NOT of
the unsigned type, interpreted as a signed integer:

>>> np.invert(np.array([13], dtype=np.int8))
array([-14], dtype=int8)
>>> np.binary_repr(-14, width=8)
'11110010'

Booleans are accepted as well:

>>> np.invert(np.array([True, False]))
array([False,  True])

The ``~`` operator can be used as a shorthand for ``np.invert`` on
ndarrays.

>>> x1 = np.array([True, False])
>>> ~x1
array([False,  True])
"""
bitwise_or: numpy.ufunc  # value = <ufunc 'bitwise_or'>
"""
bitwise_or(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the bit-wise OR of two arrays element-wise.

Computes the bit-wise OR of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``|``.

Parameters
----------
x1, x2 : array_like
    Only integer and boolean types are handled.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_or
bitwise_and
bitwise_xor
binary_repr :
    Return the binary representation of the input number as a string.

Examples
--------
>>> import numpy as np

The number 13 has the binary representation ``00001101``. Likewise,
16 is represented by ``00010000``.  The bit-wise OR of 13 and 16 is
then ``00011101``, or 29:

>>> np.bitwise_or(13, 16)
29
>>> np.binary_repr(29)
'11101'

>>> np.bitwise_or(32, 2)
34
>>> np.bitwise_or([33, 4], 1)
array([33,  5])
>>> np.bitwise_or([33, 4], [1, 2])
array([33,  6])

>>> np.bitwise_or(np.array([2, 5, 255]), np.array([4, 4, 4]))
array([  6,   5, 255])
>>> np.array([2, 5, 255]) | np.array([4, 4, 4])
array([  6,   5, 255])
>>> np.bitwise_or(np.array([2, 5, 255, 2147483647], dtype=np.int32),
...               np.array([4, 4, 4, 2147483647], dtype=np.int32))
array([         6,          5,        255, 2147483647], dtype=int32)
>>> np.bitwise_or([True, True], [False, True])
array([ True,  True])

The ``|`` operator can be used as a shorthand for ``np.bitwise_or`` on
ndarrays.

>>> x1 = np.array([2, 5, 255])
>>> x2 = np.array([4, 4, 4])
>>> x1 | x2
array([  6,   5, 255])
"""
bitwise_right_shift: numpy.ufunc  # value = <ufunc 'right_shift'>
"""
right_shift(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Shift the bits of an integer to the right.

Bits are shifted to the right `x2`.  Because the internal
representation of numbers is in binary format, this operation is
equivalent to dividing `x1` by ``2**x2``.

Parameters
----------
x1 : array_like, int
    Input values.
x2 : array_like, int
    Number of bits to remove at the right of `x1`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray, int
    Return `x1` with bits shifted `x2` times to the right.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
left_shift : Shift the bits of an integer to the left.
binary_repr : Return the binary representation of the input number
    as a string.

Examples
--------
>>> import numpy as np
>>> np.binary_repr(10)
'1010'
>>> np.right_shift(10, 1)
5
>>> np.binary_repr(5)
'101'

>>> np.right_shift(10, [1,2,3])
array([5, 2, 1])

The ``>>`` operator can be used as a shorthand for ``np.right_shift`` on
ndarrays.

>>> x1 = 10
>>> x2 = np.array([1,2,3])
>>> x1 >> x2
array([5, 2, 1])
"""
bitwise_xor: numpy.ufunc  # value = <ufunc 'bitwise_xor'>
"""
bitwise_xor(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the bit-wise XOR of two arrays element-wise.

Computes the bit-wise XOR of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``^``.

Parameters
----------
x1, x2 : array_like
    Only integer and boolean types are handled.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_xor
bitwise_and
bitwise_or
binary_repr :
    Return the binary representation of the input number as a string.

Examples
--------
>>> import numpy as np

The number 13 is represented by ``00001101``. Likewise, 17 is
represented by ``00010001``.  The bit-wise XOR of 13 and 17 is
therefore ``00011100``, or 28:

>>> np.bitwise_xor(13, 17)
28
>>> np.binary_repr(28)
'11100'

>>> np.bitwise_xor(31, 5)
26
>>> np.bitwise_xor([31,3], 5)
array([26,  6])

>>> np.bitwise_xor([31,3], [5,6])
array([26,  5])
>>> np.bitwise_xor([True, True], [False, True])
array([ True, False])

The ``^`` operator can be used as a shorthand for ``np.bitwise_xor`` on
ndarrays.

>>> x1 = np.array([True, True])
>>> x2 = np.array([False, True])
>>> x1 ^ x2
array([ True, False])
"""
c_: numpy.lib._index_tricks_impl.CClass  # value = <numpy.lib._index_tricks_impl.CClass object>
cbrt: numpy.ufunc  # value = <ufunc 'cbrt'>
"""
cbrt(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the cube-root of an array, element-wise.

Parameters
----------
x : array_like
    The values whose cube-roots are required.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    An array of the same shape as `x`, containing the
    cube root of each element in `x`.
    If `out` was provided, `y` is a reference to it.
    This is a scalar if `x` is a scalar.


Examples
--------
>>> import numpy as np
>>> np.cbrt([1,8,27])
array([ 1.,  2.,  3.])
"""
ceil: numpy.ufunc  # value = <ufunc 'ceil'>
r"""
ceil(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the ceiling of the input, element-wise.

The ceil of the scalar `x` is the smallest integer `i`, such that
``i >= x``.  It is often denoted as :math:`\lceil x \rceil`.

Parameters
----------
x : array_like
    Input data.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The ceiling of each element in `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
floor, trunc, rint, fix

Examples
--------
>>> import numpy as np

>>> a = np.array([-1.7, -1.5, -0.2, 0.2, 1.5, 1.7, 2.0])
>>> np.ceil(a)
array([-1., -1., -0.,  1.,  2.,  2.,  2.])
"""
conj: numpy.ufunc  # value = <ufunc 'conjugate'>
"""
conjugate(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the complex conjugate, element-wise.

The complex conjugate of a complex number is obtained by changing the
sign of its imaginary part.

Parameters
----------
x : array_like
    Input value.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The complex conjugate of `x`, with same dtype as `y`.
    This is a scalar if `x` is a scalar.

Notes
-----
`conj` is an alias for `conjugate`:

>>> np.conj is np.conjugate
True

Examples
--------
>>> import numpy as np
>>> np.conjugate(1+2j)
(1-2j)

>>> x = np.eye(2) + 1j * np.eye(2)
>>> np.conjugate(x)
array([[ 1.-1.j,  0.-0.j],
       [ 0.-0.j,  1.-1.j]])
"""
conjugate: numpy.ufunc  # value = <ufunc 'conjugate'>
"""
conjugate(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the complex conjugate, element-wise.

The complex conjugate of a complex number is obtained by changing the
sign of its imaginary part.

Parameters
----------
x : array_like
    Input value.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The complex conjugate of `x`, with same dtype as `y`.
    This is a scalar if `x` is a scalar.

Notes
-----
`conj` is an alias for `conjugate`:

>>> np.conj is np.conjugate
True

Examples
--------
>>> import numpy as np
>>> np.conjugate(1+2j)
(1-2j)

>>> x = np.eye(2) + 1j * np.eye(2)
>>> np.conjugate(x)
array([[ 1.-1.j,  0.-0.j],
       [ 0.-0.j,  1.-1.j]])
"""
copysign: numpy.ufunc  # value = <ufunc 'copysign'>
"""
copysign(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Change the sign of x1 to that of x2, element-wise.

If `x2` is a scalar, its sign will be copied to all elements of `x1`.

Parameters
----------
x1 : array_like
    Values to change the sign of.
x2 : array_like
    The sign of `x2` is copied to `x1`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    The values of `x1` with the sign of `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
sign
signbit

Examples
--------
>>> import numpy as np
>>> np.copysign(1.3, -1)
-1.3
>>> 1/np.copysign(0, 1)
inf
>>> 1/np.copysign(0, -1)
-inf

>>> np.copysign([-1, 0, 1], -1.1)
array([-1., -0., -1.])
>>> np.copysign([-1, 0, 1], np.arange(3)-1)
array([-1.,  0.,  1.])
"""
cos: numpy.ufunc  # value = <ufunc 'cos'>
"""
cos(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Cosine element-wise.

Parameters
----------
x : array_like
    Input array in radians.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding cosine values.
    This is a scalar if `x` is a scalar.

Notes
-----
If `out` is provided, the function writes the result into it,
and returns a reference to `out`.  (See Examples)

References
----------
M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions.
New York, NY: Dover, 1972.

Examples
--------
>>> import numpy as np
>>> np.cos(np.array([0, np.pi/2, np.pi]))
array([  1.00000000e+00,   6.12303177e-17,  -1.00000000e+00])
>>>
>>> # Example of providing the optional output parameter
>>> out1 = np.array([0], dtype='d')
>>> out2 = np.cos([0.1], out1)
>>> out2 is out1
True
>>>
>>> # Example of ValueError due to provision of shape mis-matched `out`
>>> np.cos(np.zeros((3,3)),np.zeros((2,2)))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: operands could not be broadcast together with shapes (3,3) (2,2)
"""
cosh: numpy.ufunc  # value = <ufunc 'cosh'>
"""
cosh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Hyperbolic cosine, element-wise.

Equivalent to ``1/2 * (np.exp(x) + np.exp(-x))`` and ``np.cos(1j*x)``.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array of same shape as `x`.
    This is a scalar if `x` is a scalar.

Examples
--------
>>> import numpy as np
>>> np.cosh(0)
1.0

The hyperbolic cosine describes the shape of a hanging cable:

>>> import matplotlib.pyplot as plt
>>> x = np.linspace(-4, 4, 1000)
>>> plt.plot(x, np.cosh(x))
>>> plt.show()
"""
deg2rad: numpy.ufunc  # value = <ufunc 'deg2rad'>
"""
deg2rad(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Convert angles from degrees to radians.

Parameters
----------
x : array_like
    Angles in degrees.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding angle in radians.
    This is a scalar if `x` is a scalar.

See Also
--------
rad2deg : Convert angles from radians to degrees.
unwrap : Remove large jumps in angle by wrapping.

Notes
-----
``deg2rad(x)`` is ``x * pi / 180``.

Examples
--------
>>> import numpy as np
>>> np.deg2rad(180)
3.1415926535897931
"""
degrees: numpy.ufunc  # value = <ufunc 'degrees'>
"""
degrees(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Convert angles from radians to degrees.

Parameters
----------
x : array_like
    Input array in radians.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray of floats
    The corresponding degree values; if `out` was supplied this is a
    reference to it.
    This is a scalar if `x` is a scalar.

See Also
--------
rad2deg : equivalent function

Examples
--------
Convert a radian array to degrees

>>> import numpy as np
>>> rad = np.arange(12.)*np.pi/6
>>> np.degrees(rad)
array([   0.,   30.,   60.,   90.,  120.,  150.,  180.,  210.,  240.,
        270.,  300.,  330.])

>>> out = np.zeros((rad.shape))
>>> r = np.degrees(rad, out)
>>> np.all(r == out)
True
"""
divide: numpy.ufunc  # value = <ufunc 'divide'>
"""
divide(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Divide arguments element-wise.

Parameters
----------
x1 : array_like
    Dividend array.
x2 : array_like
    Divisor array.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The quotient ``x1/x2``, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
seterr : Set whether to raise or warn on overflow, underflow and
         division by zero.

Notes
-----
Equivalent to ``x1`` / ``x2`` in terms of array-broadcasting.

The ``true_divide(x1, x2)`` function is an alias for
``divide(x1, x2)``.

Examples
--------
>>> import numpy as np
>>> np.divide(2.0, 4.0)
0.5
>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> np.divide(x1, x2)
array([[nan, 1. , 1. ],
       [inf, 4. , 2.5],
       [inf, 7. , 4. ]])

The ``/`` operator can be used as a shorthand for ``np.divide`` on
ndarrays.

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = 2 * np.ones(3)
>>> x1 / x2
array([[0. , 0.5, 1. ],
       [1.5, 2. , 2.5],
       [3. , 3.5, 4. ]])
"""
divmod: numpy.ufunc  # value = <ufunc 'divmod'>
"""
divmod(x1, x2[, out1, out2], / [, out=(None, None)], *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return element-wise quotient and remainder simultaneously.

``np.divmod(x, y)`` is equivalent to ``(x // y, x % y)``, but faster
because it avoids redundant work. It is used to implement the Python
built-in function ``divmod`` on NumPy arrays.

Parameters
----------
x1 : array_like
    Dividend array.
x2 : array_like
    Divisor array.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out1 : ndarray
    Element-wise quotient resulting from floor division.
    This is a scalar if both `x1` and `x2` are scalars.
out2 : ndarray
    Element-wise remainder from floor division.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
floor_divide : Equivalent to Python's ``//`` operator.
remainder : Equivalent to Python's ``%`` operator.
modf : Equivalent to ``divmod(x, 1)`` for positive ``x`` with the return
       values switched.

Examples
--------
>>> import numpy as np
>>> np.divmod(np.arange(5), 3)
(array([0, 0, 0, 1, 1]), array([0, 1, 2, 0, 1]))

The `divmod` function can be used as a shorthand for ``np.divmod`` on
ndarrays.

>>> x = np.arange(5)
>>> divmod(x, 3)
(array([0, 0, 0, 1, 1]), array([0, 1, 2, 0, 1]))
"""
e: float = 2.718281828459045
equal: numpy.ufunc  # value = <ufunc 'equal'>
"""
equal(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return (x1 == x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
not_equal, greater_equal, less_equal, greater, less

Examples
--------
>>> import numpy as np
>>> np.equal([0, 1, 3], np.arange(3))
array([ True,  True, False])

What is compared are values, not types. So an int (1) and an array of
length one can evaluate as True:

>>> np.equal(1, np.ones(1))
array([ True])

The ``==`` operator can be used as a shorthand for ``np.equal`` on
ndarrays.

>>> a = np.array([2, 4, 6])
>>> b = np.array([2, 4, 2])
>>> a == b
array([ True,  True, False])
"""
euler_gamma: float = 0.5772156649015329
exp: numpy.ufunc  # value = <ufunc 'exp'>
r"""
exp(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Calculate the exponential of all elements in the input array.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise exponential of `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
expm1 : Calculate ``exp(x) - 1`` for all elements in the array.
exp2  : Calculate ``2**x`` for all elements in the array.

Notes
-----
The irrational number ``e`` is also known as Euler's number.  It is
approximately 2.718281, and is the base of the natural logarithm,
``ln`` (this means that, if :math:`x = \ln y = \log_e y`,
then :math:`e^x = y`. For real input, ``exp(x)`` is always positive.

For complex arguments, ``x = a + ib``, we can write
:math:`e^x = e^a e^{ib}`.  The first term, :math:`e^a`, is already
known (it is the real argument, described above).  The second term,
:math:`e^{ib}`, is :math:`\cos b + i \sin b`, a function with
magnitude 1 and a periodic phase.

References
----------
.. [1] Wikipedia, "Exponential function",
       https://en.wikipedia.org/wiki/Exponential_function
.. [2] M. Abramovitz and I. A. Stegun, "Handbook of Mathematical Functions
       with Formulas, Graphs, and Mathematical Tables," Dover, 1964, p. 69,
       https://personal.math.ubc.ca/~cbm/aands/page_69.htm

Examples
--------
Plot the magnitude and phase of ``exp(x)`` in the complex plane:

>>> import numpy as np

>>> import matplotlib.pyplot as plt
>>> import numpy as np

>>> x = np.linspace(-2*np.pi, 2*np.pi, 100)
>>> xx = x + 1j * x[:, np.newaxis] # a + ib over complex plane
>>> out = np.exp(xx)

>>> plt.subplot(121)
>>> plt.imshow(np.abs(out),
...            extent=[-2*np.pi, 2*np.pi, -2*np.pi, 2*np.pi], cmap='gray')
>>> plt.title('Magnitude of exp(x)')

>>> plt.subplot(122)
>>> plt.imshow(np.angle(out),
...            extent=[-2*np.pi, 2*np.pi, -2*np.pi, 2*np.pi], cmap='hsv')
>>> plt.title('Phase (angle) of exp(x)')
>>> plt.show()
"""
exp2: numpy.ufunc  # value = <ufunc 'exp2'>
"""
exp2(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Calculate `2**p` for all `p` in the input array.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Element-wise 2 to the power `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
power

Examples
--------
>>> import numpy as np
>>> np.exp2([2, 3])
array([ 4.,  8.])
"""
expm1: numpy.ufunc  # value = <ufunc 'expm1'>
"""
expm1(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Calculate ``exp(x) - 1`` for all elements in the array.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Element-wise exponential minus one: ``out = exp(x) - 1``.
    This is a scalar if `x` is a scalar.

See Also
--------
log1p : ``log(1 + x)``, the inverse of expm1.


Notes
-----
This function provides greater precision than ``exp(x) - 1``
for small values of ``x``.

Examples
--------

The true value of ``exp(1e-10) - 1`` is ``1.00000000005e-10`` to
about 32 significant digits. This example shows the superiority of
expm1 in this case.

>>> import numpy as np
>>> np.expm1(1e-10)
1.00000000005e-10
>>> np.exp(1e-10) - 1
1.000000082740371e-10
"""
fabs: numpy.ufunc  # value = <ufunc 'fabs'>
"""
fabs(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the absolute values element-wise.

This function returns the absolute values (positive magnitude) of the
data in `x`. Complex values are not handled, use `absolute` to find the
absolute values of complex data.

Parameters
----------
x : array_like
    The array of numbers for which the absolute values are required. If
    `x` is a scalar, the result `y` will also be a scalar.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The absolute values of `x`, the returned values are always floats.
    This is a scalar if `x` is a scalar.

See Also
--------
absolute : Absolute values including `complex` types.

Examples
--------
>>> import numpy as np
>>> np.fabs(-1)
1.0
>>> np.fabs([-1.2, 1.2])
array([ 1.2,  1.2])
"""
float_power: numpy.ufunc  # value = <ufunc 'float_power'>
"""
float_power(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

First array elements raised to powers from second array, element-wise.

Raise each base in `x1` to the positionally-corresponding power in `x2`.
`x1` and `x2` must be broadcastable to the same shape. This differs from
the power function in that integers, float16, and float32  are promoted to
floats with a minimum precision of float64 so that the result is always
inexact.  The intent is that the function will return a usable result for
negative powers and seldom overflow for positive powers.

Negative values raised to a non-integral value will return ``nan``.
To get complex results, cast the input to complex, or specify the
``dtype`` to be ``complex`` (see the example below).

Parameters
----------
x1 : array_like
    The bases.
x2 : array_like
    The exponents.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The bases in `x1` raised to the exponents in `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
power : power function that preserves type

Examples
--------
>>> import numpy as np

Cube each element in a list.

>>> x1 = range(6)
>>> x1
[0, 1, 2, 3, 4, 5]
>>> np.float_power(x1, 3)
array([   0.,    1.,    8.,   27.,   64.,  125.])

Raise the bases to different exponents.

>>> x2 = [1.0, 2.0, 3.0, 3.0, 2.0, 1.0]
>>> np.float_power(x1, x2)
array([  0.,   1.,   8.,  27.,  16.,   5.])

The effect of broadcasting.

>>> x2 = np.array([[1, 2, 3, 3, 2, 1], [1, 2, 3, 3, 2, 1]])
>>> x2
array([[1, 2, 3, 3, 2, 1],
       [1, 2, 3, 3, 2, 1]])
>>> np.float_power(x1, x2)
array([[  0.,   1.,   8.,  27.,  16.,   5.],
       [  0.,   1.,   8.,  27.,  16.,   5.]])

Negative values raised to a non-integral value will result in ``nan``
(and a warning will be generated).

>>> x3 = np.array([-1, -4])
>>> with np.errstate(invalid='ignore'):
...     p = np.float_power(x3, 1.5)
...
>>> p
array([nan, nan])

To get complex results, give the argument ``dtype=complex``.

>>> np.float_power(x3, 1.5, dtype=complex)
array([-1.83697020e-16-1.j, -1.46957616e-15-8.j])
"""
floor: numpy.ufunc  # value = <ufunc 'floor'>
r"""
floor(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the floor of the input, element-wise.

The floor of the scalar `x` is the largest integer `i`, such that
`i <= x`.  It is often denoted as :math:`\lfloor x \rfloor`.

Parameters
----------
x : array_like
    Input data.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The floor of each element in `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
ceil, trunc, rint, fix

Notes
-----
Some spreadsheet programs calculate the "floor-towards-zero", where
``floor(-2.5) == -2``.  NumPy instead uses the definition of
`floor` where `floor(-2.5) == -3`. The "floor-towards-zero"
function is called ``fix`` in NumPy.

Examples
--------
>>> import numpy as np
>>> a = np.array([-1.7, -1.5, -0.2, 0.2, 1.5, 1.7, 2.0])
>>> np.floor(a)
array([-2., -2., -1.,  0.,  1.,  1.,  2.])
"""
floor_divide: numpy.ufunc  # value = <ufunc 'floor_divide'>
"""
floor_divide(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the largest integer smaller or equal to the division of the inputs.
It is equivalent to the Python ``//`` operator and pairs with the
Python ``%`` (`remainder`), function so that ``a = a % b + b * (a // b)``
up to roundoff.

Parameters
----------
x1 : array_like
    Numerator.
x2 : array_like
    Denominator.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    y = floor(`x1`/`x2`)
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
remainder : Remainder complementary to floor_divide.
divmod : Simultaneous floor division and remainder.
divide : Standard division.
floor : Round a number to the nearest integer toward minus infinity.
ceil : Round a number to the nearest integer toward infinity.

Examples
--------
>>> import numpy as np
>>> np.floor_divide(7,3)
2
>>> np.floor_divide([1., 2., 3., 4.], 2.5)
array([ 0.,  0.,  1.,  1.])

The ``//`` operator can be used as a shorthand for ``np.floor_divide``
on ndarrays.

>>> x1 = np.array([1., 2., 3., 4.])
>>> x1 // 2.5
array([0., 0., 1., 1.])
"""
fmax: numpy.ufunc  # value = <ufunc 'fmax'>
"""
fmax(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise maximum of array elements.

Compare two arrays and return a new array containing the element-wise
maxima. If one of the elements being compared is a NaN, then the
non-nan element is returned. If both elements are NaNs then the first
is returned.  The latter distinction is important for complex NaNs,
which are defined as at least one of the real or imaginary parts being
a NaN. The net effect is that NaNs are ignored when possible.

Parameters
----------
x1, x2 : array_like
    The arrays holding the elements to be compared.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The maximum of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
fmin :
    Element-wise minimum of two arrays, ignores NaNs.
maximum :
    Element-wise maximum of two arrays, propagates NaNs.
amax :
    The maximum value of an array along a given axis, propagates NaNs.
nanmax :
    The maximum value of an array along a given axis, ignores NaNs.

minimum, amin, nanmin

Notes
-----
The fmax is equivalent to ``np.where(x1 >= x2, x1, x2)`` when neither
x1 nor x2 are NaNs, but it is faster and does proper broadcasting.

Examples
--------
>>> import numpy as np
>>> np.fmax([2, 3, 4], [1, 5, 2])
array([ 2,  5,  4])

>>> np.fmax(np.eye(2), [0.5, 2])
array([[ 1. ,  2. ],
       [ 0.5,  2. ]])

>>> np.fmax([np.nan, 0, np.nan],[0, np.nan, np.nan])
array([ 0.,  0., nan])
"""
fmin: numpy.ufunc  # value = <ufunc 'fmin'>
"""
fmin(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise minimum of array elements.

Compare two arrays and return a new array containing the element-wise
minima. If one of the elements being compared is a NaN, then the
non-nan element is returned. If both elements are NaNs then the first
is returned.  The latter distinction is important for complex NaNs,
which are defined as at least one of the real or imaginary parts being
a NaN. The net effect is that NaNs are ignored when possible.

Parameters
----------
x1, x2 : array_like
    The arrays holding the elements to be compared.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The minimum of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
fmax :
    Element-wise maximum of two arrays, ignores NaNs.
minimum :
    Element-wise minimum of two arrays, propagates NaNs.
amin :
    The minimum value of an array along a given axis, propagates NaNs.
nanmin :
    The minimum value of an array along a given axis, ignores NaNs.

maximum, amax, nanmax

Notes
-----
The fmin is equivalent to ``np.where(x1 <= x2, x1, x2)`` when neither
x1 nor x2 are NaNs, but it is faster and does proper broadcasting.

Examples
--------
>>> import numpy as np
>>> np.fmin([2, 3, 4], [1, 5, 2])
array([1, 3, 2])

>>> np.fmin(np.eye(2), [0.5, 2])
array([[ 0.5,  0. ],
       [ 0. ,  1. ]])

>>> np.fmin([np.nan, 0, np.nan],[0, np.nan, np.nan])
array([ 0.,  0., nan])
"""
fmod: numpy.ufunc  # value = <ufunc 'fmod'>
"""
fmod(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns the element-wise remainder of division.

This is the NumPy implementation of the C library function fmod, the
remainder has the same sign as the dividend `x1`. It is equivalent to
the Matlab(TM) ``rem`` function and should not be confused with the
Python modulus operator ``x1 % x2``.

Parameters
----------
x1 : array_like
    Dividend.
x2 : array_like
    Divisor.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : array_like
    The remainder of the division of `x1` by `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
remainder : Equivalent to the Python ``%`` operator.
divide

Notes
-----
The result of the modulo operation for negative dividend and divisors
is bound by conventions. For `fmod`, the sign of result is the sign of
the dividend, while for `remainder` the sign of the result is the sign
of the divisor. The `fmod` function is equivalent to the Matlab(TM)
``rem`` function.

Examples
--------
>>> import numpy as np
>>> np.fmod([-3, -2, -1, 1, 2, 3], 2)
array([-1,  0, -1,  1,  0,  1])
>>> np.remainder([-3, -2, -1, 1, 2, 3], 2)
array([1, 0, 1, 1, 0, 1])

>>> np.fmod([5, 3], [2, 2.])
array([ 1.,  1.])
>>> a = np.arange(-3, 3).reshape(3, 2)
>>> a
array([[-3, -2],
       [-1,  0],
       [ 1,  2]])
>>> np.fmod(a, [2,2])
array([[-1,  0],
       [-1,  0],
       [ 1,  0]])
"""
frexp: numpy.ufunc  # value = <ufunc 'frexp'>
"""
frexp(x[, out1, out2], / [, out=(None, None)], *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Decompose the elements of x into mantissa and twos exponent.

Returns (`mantissa`, `exponent`), where ``x = mantissa * 2**exponent``.
The mantissa lies in the open interval(-1, 1), while the twos
exponent is a signed integer.

Parameters
----------
x : array_like
    Array of numbers to be decomposed.
out1 : ndarray, optional
    Output array for the mantissa. Must have the same shape as `x`.
out2 : ndarray, optional
    Output array for the exponent. Must have the same shape as `x`.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
mantissa : ndarray
    Floating values between -1 and 1.
    This is a scalar if `x` is a scalar.
exponent : ndarray
    Integer exponents of 2.
    This is a scalar if `x` is a scalar.

See Also
--------
ldexp : Compute ``y = x1 * 2**x2``, the inverse of `frexp`.

Notes
-----
Complex dtypes are not supported, they will raise a TypeError.

Examples
--------
>>> import numpy as np
>>> x = np.arange(9)
>>> y1, y2 = np.frexp(x)
>>> y1
array([ 0.   ,  0.5  ,  0.5  ,  0.75 ,  0.5  ,  0.625,  0.75 ,  0.875,
        0.5  ])
>>> y2
array([0, 1, 2, 2, 3, 3, 3, 3, 4], dtype=int32)
>>> y1 * 2**y2
array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.])
"""
gcd: numpy.ufunc  # value = <ufunc 'gcd'>
"""
gcd(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns the greatest common divisor of ``|x1|`` and ``|x2|``

Parameters
----------
x1, x2 : array_like, int
    Arrays of values.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).

Returns
-------
y : ndarray or scalar
    The greatest common divisor of the absolute value of the inputs
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
lcm : The lowest common multiple

Examples
--------
>>> import numpy as np
>>> np.gcd(12, 20)
4
>>> np.gcd.reduce([15, 25, 35])
5
>>> np.gcd(np.arange(6), 20)
array([20,  1,  2,  1,  4,  5])
"""
greater: numpy.ufunc  # value = <ufunc 'greater'>
"""
greater(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the truth value of (x1 > x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.


See Also
--------
greater_equal, less, less_equal, equal, not_equal

Examples
--------
>>> import numpy as np
>>> np.greater([4,2],[2,2])
array([ True, False])

The ``>`` operator can be used as a shorthand for ``np.greater`` on
ndarrays.

>>> a = np.array([4, 2])
>>> b = np.array([2, 2])
>>> a > b
array([ True, False])
"""
greater_equal: numpy.ufunc  # value = <ufunc 'greater_equal'>
"""
greater_equal(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the truth value of (x1 >= x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : bool or ndarray of bool
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
greater, less, less_equal, equal, not_equal

Examples
--------
>>> import numpy as np
>>> np.greater_equal([4, 2, 1], [2, 2, 2])
array([ True, True, False])

The ``>=`` operator can be used as a shorthand for ``np.greater_equal``
on ndarrays.

>>> a = np.array([4, 2, 1])
>>> b = np.array([2, 2, 2])
>>> a >= b
array([ True,  True, False])
"""
heaviside: numpy.ufunc  # value = <ufunc 'heaviside'>
"""
heaviside(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the Heaviside step function.

The Heaviside step function [1]_ is defined as::

                          0   if x1 < 0
    heaviside(x1, x2) =  x2   if x1 == 0
                          1   if x1 > 0

where `x2` is often taken to be 0.5, but 0 and 1 are also sometimes used.

Parameters
----------
x1 : array_like
    Input values.
x2 : array_like
    The value of the function when x1 is 0.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    The output array, element-wise Heaviside step function of `x1`.
    This is a scalar if both `x1` and `x2` are scalars.

References
----------
.. [1] Wikipedia, "Heaviside step function",
       https://en.wikipedia.org/wiki/Heaviside_step_function

Examples
--------
>>> import numpy as np
>>> np.heaviside([-1.5, 0, 2.0], 0.5)
array([ 0. ,  0.5,  1. ])
>>> np.heaviside([-1.5, 0, 2.0], 1)
array([ 0.,  1.,  1.])
"""
hypot: numpy.ufunc  # value = <ufunc 'hypot'>
"""
hypot(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Given the "legs" of a right triangle, return its hypotenuse.

Equivalent to ``sqrt(x1**2 + x2**2)``, element-wise.  If `x1` or
`x2` is scalar_like (i.e., unambiguously cast-able to a scalar type),
it is broadcast for use with each element of the other argument.
(See Examples)

Parameters
----------
x1, x2 : array_like
    Leg of the triangle(s).
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
z : ndarray
    The hypotenuse of the triangle(s).
    This is a scalar if both `x1` and `x2` are scalars.

Examples
--------
>>> import numpy as np
>>> np.hypot(3*np.ones((3, 3)), 4*np.ones((3, 3)))
array([[ 5.,  5.,  5.],
       [ 5.,  5.,  5.],
       [ 5.,  5.,  5.]])

Example showing broadcast of scalar_like argument:

>>> np.hypot(3*np.ones((3, 3)), [4])
array([[ 5.,  5.,  5.],
       [ 5.,  5.,  5.],
       [ 5.,  5.,  5.]])
"""
index_exp: numpy.lib._index_tricks_impl.IndexExpression  # value = <numpy.lib._index_tricks_impl.IndexExpression object>
inf: float  # value = inf
invert: numpy.ufunc  # value = <ufunc 'invert'>
"""
invert(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute bit-wise inversion, or bit-wise NOT, element-wise.

Computes the bit-wise NOT of the underlying binary representation of
the integers in the input arrays. This ufunc implements the C/Python
operator ``~``.

For signed integer inputs, the bit-wise NOT of the absolute value is
returned. In a two's-complement system, this operation effectively flips
all the bits, resulting in a representation that corresponds to the
negative of the input plus one. This is the most common method of
representing signed integers on computers [1]_. A N-bit two's-complement
system can represent every integer in the range :math:`-2^{N-1}` to
:math:`+2^{N-1}-1`.

Parameters
----------
x : array_like
    Only integer and boolean types are handled.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Result.
    This is a scalar if `x` is a scalar.

See Also
--------
bitwise_and, bitwise_or, bitwise_xor
logical_not
binary_repr :
    Return the binary representation of the input number as a string.

Notes
-----
``numpy.bitwise_not`` is an alias for `invert`:

>>> np.bitwise_not is np.invert
True

References
----------
.. [1] Wikipedia, "Two's complement",
    https://en.wikipedia.org/wiki/Two's_complement

Examples
--------
>>> import numpy as np

We've seen that 13 is represented by ``00001101``.
The invert or bit-wise NOT of 13 is then:

>>> x = np.invert(np.array(13, dtype=np.uint8))
>>> x
np.uint8(242)
>>> np.binary_repr(x, width=8)
'11110010'

The result depends on the bit-width:

>>> x = np.invert(np.array(13, dtype=np.uint16))
>>> x
np.uint16(65522)
>>> np.binary_repr(x, width=16)
'1111111111110010'

When using signed integer types, the result is the bit-wise NOT of
the unsigned type, interpreted as a signed integer:

>>> np.invert(np.array([13], dtype=np.int8))
array([-14], dtype=int8)
>>> np.binary_repr(-14, width=8)
'11110010'

Booleans are accepted as well:

>>> np.invert(np.array([True, False]))
array([False,  True])

The ``~`` operator can be used as a shorthand for ``np.invert`` on
ndarrays.

>>> x1 = np.array([True, False])
>>> ~x1
array([False,  True])
"""
isfinite: numpy.ufunc  # value = <ufunc 'isfinite'>
"""
isfinite(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Test element-wise for finiteness (not infinity and not Not a Number).

The result is returned as a boolean array.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray, bool
    True where ``x`` is not positive infinity, negative infinity,
    or NaN; false otherwise.
    This is a scalar if `x` is a scalar.

See Also
--------
isinf, isneginf, isposinf, isnan

Notes
-----
Not a Number, positive infinity and negative infinity are considered
to be non-finite.

NumPy uses the IEEE Standard for Binary Floating-Point for Arithmetic
(IEEE 754). This means that Not a Number is not equivalent to infinity.
Also that positive infinity is not equivalent to negative infinity. But
infinity is equivalent to positive infinity.  Errors result if the
second argument is also supplied when `x` is a scalar input, or if
first and second arguments have different shapes.

Examples
--------
>>> import numpy as np
>>> np.isfinite(1)
True
>>> np.isfinite(0)
True
>>> np.isfinite(np.nan)
False
>>> np.isfinite(np.inf)
False
>>> np.isfinite(-np.inf)
False
>>> np.isfinite([np.log(-1.),1.,np.log(0)])
array([False,  True, False])

>>> x = np.array([-np.inf, 0., np.inf])
>>> y = np.array([2, 2, 2])
>>> np.isfinite(x, y)
array([0, 1, 0])
>>> y
array([0, 1, 0])
"""
isinf: numpy.ufunc  # value = <ufunc 'isinf'>
"""
isinf(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Test element-wise for positive or negative infinity.

Returns a boolean array of the same shape as `x`, True where ``x ==
+/-inf``, otherwise False.

Parameters
----------
x : array_like
    Input values
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : bool (scalar) or boolean ndarray
    True where ``x`` is positive or negative infinity, false otherwise.
    This is a scalar if `x` is a scalar.

See Also
--------
isneginf, isposinf, isnan, isfinite

Notes
-----
NumPy uses the IEEE Standard for Binary Floating-Point for Arithmetic
(IEEE 754).

Errors result if the second argument is supplied when the first
argument is a scalar, or if the first and second arguments have
different shapes.

Examples
--------
>>> import numpy as np
>>> np.isinf(np.inf)
True
>>> np.isinf(np.nan)
False
>>> np.isinf(-np.inf)
True
>>> np.isinf([np.inf, -np.inf, 1.0, np.nan])
array([ True,  True, False, False])

>>> x = np.array([-np.inf, 0., np.inf])
>>> y = np.array([2, 2, 2])
>>> np.isinf(x, y)
array([1, 0, 1])
>>> y
array([1, 0, 1])
"""
isnan: numpy.ufunc  # value = <ufunc 'isnan'>
"""
isnan(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Test element-wise for NaN and return result as a boolean array.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or bool
    True where ``x`` is NaN, false otherwise.
    This is a scalar if `x` is a scalar.

See Also
--------
isinf, isneginf, isposinf, isfinite, isnat

Notes
-----
NumPy uses the IEEE Standard for Binary Floating-Point for Arithmetic
(IEEE 754). This means that Not a Number is not equivalent to infinity.

Examples
--------
>>> import numpy as np
>>> np.isnan(np.nan)
True
>>> np.isnan(np.inf)
False
>>> np.isnan([np.log(-1.),1.,np.log(0)])
array([ True, False, False])
"""
isnat: numpy.ufunc  # value = <ufunc 'isnat'>
"""
isnat(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Test element-wise for NaT (not a time) and return result as a boolean array.

Parameters
----------
x : array_like
    Input array with datetime or timedelta data type.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or bool
    True where ``x`` is NaT, false otherwise.
    This is a scalar if `x` is a scalar.

See Also
--------
isnan, isinf, isneginf, isposinf, isfinite

Examples
--------
>>> import numpy as np
>>> np.isnat(np.datetime64("NaT"))
True
>>> np.isnat(np.datetime64("2016-01-01"))
False
>>> np.isnat(np.array(["NaT", "2016-01-01"], dtype="datetime64[ns]"))
array([ True, False])
"""
lcm: numpy.ufunc  # value = <ufunc 'lcm'>
"""
lcm(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns the lowest common multiple of ``|x1|`` and ``|x2|``

Parameters
----------
x1, x2 : array_like, int
    Arrays of values.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).

Returns
-------
y : ndarray or scalar
    The lowest common multiple of the absolute value of the inputs
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
gcd : The greatest common divisor

Examples
--------
>>> import numpy as np
>>> np.lcm(12, 20)
60
>>> np.lcm.reduce([3, 12, 20])
60
>>> np.lcm.reduce([40, 12, 20])
120
>>> np.lcm(np.arange(6), 20)
array([ 0, 20, 20, 60, 20, 20])
"""
ldexp: numpy.ufunc  # value = <ufunc 'ldexp'>
"""
ldexp(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns x1 * 2**x2, element-wise.

The mantissas `x1` and twos exponents `x2` are used to construct
floating point numbers ``x1 * 2**x2``.

Parameters
----------
x1 : array_like
    Array of multipliers.
x2 : array_like, int
    Array of twos exponents.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The result of ``x1 * 2**x2``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
frexp : Return (y1, y2) from ``x = y1 * 2**y2``, inverse to `ldexp`.

Notes
-----
Complex dtypes are not supported, they will raise a TypeError.

`ldexp` is useful as the inverse of `frexp`, if used by itself it is
more clear to simply use the expression ``x1 * 2**x2``.

Examples
--------
>>> import numpy as np
>>> np.ldexp(5, np.arange(4))
array([ 5., 10., 20., 40.], dtype=float16)

>>> x = np.arange(6)
>>> np.ldexp(*np.frexp(x))
array([ 0.,  1.,  2.,  3.,  4.,  5.])
"""
left_shift: numpy.ufunc  # value = <ufunc 'left_shift'>
"""
left_shift(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Shift the bits of an integer to the left.

Bits are shifted to the left by appending `x2` 0s at the right of `x1`.
Since the internal representation of numbers is in binary format, this
operation is equivalent to multiplying `x1` by ``2**x2``.

Parameters
----------
x1 : array_like of integer type
    Input values.
x2 : array_like of integer type
    Number of zeros to append to `x1`. Has to be non-negative.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : array of integer type
    Return `x1` with bits shifted `x2` times to the left.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
right_shift : Shift the bits of an integer to the right.
binary_repr : Return the binary representation of the input number
    as a string.

Examples
--------
>>> import numpy as np
>>> np.binary_repr(5)
'101'
>>> np.left_shift(5, 2)
20
>>> np.binary_repr(20)
'10100'

>>> np.left_shift(5, [1,2,3])
array([10, 20, 40])

Note that the dtype of the second argument may change the dtype of the
result and can lead to unexpected results in some cases (see
:ref:`Casting Rules <ufuncs.casting>`):

>>> a = np.left_shift(np.uint8(255), np.int64(1))  # Expect 254
>>> print(a, type(a)) # Unexpected result due to upcasting
510 <class 'numpy.int64'>
>>> b = np.left_shift(np.uint8(255), np.uint8(1))
>>> print(b, type(b))
254 <class 'numpy.uint8'>

The ``<<`` operator can be used as a shorthand for ``np.left_shift`` on
ndarrays.

>>> x1 = 5
>>> x2 = np.array([1, 2, 3])
>>> x1 << x2
array([10, 20, 40])
"""
less: numpy.ufunc  # value = <ufunc 'less'>
"""
less(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the truth value of (x1 < x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
greater, less_equal, greater_equal, equal, not_equal

Examples
--------
>>> import numpy as np
>>> np.less([1, 2], [2, 2])
array([ True, False])

The ``<`` operator can be used as a shorthand for ``np.less`` on ndarrays.

>>> a = np.array([1, 2])
>>> b = np.array([2, 2])
>>> a < b
array([ True, False])
"""
less_equal: numpy.ufunc  # value = <ufunc 'less_equal'>
"""
less_equal(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the truth value of (x1 <= x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
greater, less, greater_equal, equal, not_equal

Examples
--------
>>> import numpy as np
>>> np.less_equal([4, 2, 1], [2, 2, 2])
array([False,  True,  True])

The ``<=`` operator can be used as a shorthand for ``np.less_equal`` on
ndarrays.

>>> a = np.array([4, 2, 1])
>>> b = np.array([2, 2, 2])
>>> a <= b
array([False,  True,  True])
"""
little_endian: bool = True
log: numpy.ufunc  # value = <ufunc 'log'>
"""
log(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Natural logarithm, element-wise.

The natural logarithm `log` is the inverse of the exponential function,
so that `log(exp(x)) = x`. The natural logarithm is logarithm in base
`e`.

Parameters
----------
x : array_like
    Input value.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The natural logarithm of `x`, element-wise.
    This is a scalar if `x` is a scalar.

See Also
--------
log10, log2, log1p, emath.log

Notes
-----
Logarithm is a multivalued function: for each `x` there is an infinite
number of `z` such that `exp(z) = x`. The convention is to return the
`z` whose imaginary part lies in `(-pi, pi]`.

For real-valued input data types, `log` always returns real output. For
each value that cannot be expressed as a real number or infinity, it
yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `log` is a complex analytical function that
has a branch cut `[-inf, 0]` and is continuous from above on it. `log`
handles the floating-point negative zero as an infinitesimal negative
number, conforming to the C99 standard.

In the cases where the input has a negative real part and a very small
negative complex part (approaching 0), the result is so close to `-pi`
that it evaluates to exactly `-pi`.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 67.
       https://personal.math.ubc.ca/~cbm/aands/page_67.htm
.. [2] Wikipedia, "Logarithm". https://en.wikipedia.org/wiki/Logarithm

Examples
--------
>>> import numpy as np
>>> np.log([1, np.e, np.e**2, 0])
array([  0.,   1.,   2., -inf])
"""
log10: numpy.ufunc  # value = <ufunc 'log10'>
"""
log10(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the base 10 logarithm of the input array, element-wise.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The logarithm to the base 10 of `x`, element-wise. NaNs are
    returned where x is negative.
    This is a scalar if `x` is a scalar.

See Also
--------
emath.log10

Notes
-----
Logarithm is a multivalued function: for each `x` there is an infinite
number of `z` such that `10**z = x`. The convention is to return the
`z` whose imaginary part lies in `(-pi, pi]`.

For real-valued input data types, `log10` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `log10` is a complex analytical function that
has a branch cut `[-inf, 0]` and is continuous from above on it.
`log10` handles the floating-point negative zero as an infinitesimal
negative number, conforming to the C99 standard.

In the cases where the input has a negative real part and a very small
negative complex part (approaching 0), the result is so close to `-pi`
that it evaluates to exactly `-pi`.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 67.
       https://personal.math.ubc.ca/~cbm/aands/page_67.htm
.. [2] Wikipedia, "Logarithm". https://en.wikipedia.org/wiki/Logarithm

Examples
--------
>>> import numpy as np
>>> np.log10([1e-15, -3.])
array([-15.,  nan])
"""
log1p: numpy.ufunc  # value = <ufunc 'log1p'>
"""
log1p(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the natural logarithm of one plus the input array, element-wise.

Calculates ``log(1 + x)``.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    Natural logarithm of `1 + x`, element-wise.
    This is a scalar if `x` is a scalar.

See Also
--------
expm1 : ``exp(x) - 1``, the inverse of `log1p`.

Notes
-----
For real-valued input, `log1p` is accurate also for `x` so small
that `1 + x == 1` in floating-point accuracy.

Logarithm is a multivalued function: for each `x` there is an infinite
number of `z` such that `exp(z) = 1 + x`. The convention is to return
the `z` whose imaginary part lies in `[-pi, pi]`.

For real-valued input data types, `log1p` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `log1p` is a complex analytical function that
has a branch cut `[-inf, -1]` and is continuous from above on it.
`log1p` handles the floating-point negative zero as an infinitesimal
negative number, conforming to the C99 standard.

References
----------
.. [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
       10th printing, 1964, pp. 67.
       https://personal.math.ubc.ca/~cbm/aands/page_67.htm
.. [2] Wikipedia, "Logarithm". https://en.wikipedia.org/wiki/Logarithm

Examples
--------
>>> import numpy as np
>>> np.log1p(1e-99)
1e-99
>>> np.log(1 + 1e-99)
0.0
"""
log2: numpy.ufunc  # value = <ufunc 'log2'>
"""
log2(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Base-2 logarithm of `x`.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    Base-2 logarithm of `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
log, log10, log1p, emath.log2

Notes
-----
Logarithm is a multivalued function: for each `x` there is an infinite
number of `z` such that `2**z = x`. The convention is to return the `z`
whose imaginary part lies in `(-pi, pi]`.

For real-valued input data types, `log2` always returns real output.
For each value that cannot be expressed as a real number or infinity,
it yields ``nan`` and sets the `invalid` floating point error flag.

For complex-valued input, `log2` is a complex analytical function that
has a branch cut `[-inf, 0]` and is continuous from above on it. `log2`
handles the floating-point negative zero as an infinitesimal negative
number, conforming to the C99 standard.

In the cases where the input has a negative real part and a very small
negative complex part (approaching 0), the result is so close to `-pi`
that it evaluates to exactly `-pi`.

Examples
--------
>>> import numpy as np
>>> x = np.array([0, 1, 2, 2**4])
>>> np.log2(x)
array([-inf,   0.,   1.,   4.])

>>> xi = np.array([0+1.j, 1, 2+0.j, 4.j])
>>> np.log2(xi)
array([ 0.+2.26618007j,  0.+0.j        ,  1.+0.j        ,  2.+2.26618007j])
"""
logaddexp: numpy.ufunc  # value = <ufunc 'logaddexp'>
"""
logaddexp(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Logarithm of the sum of exponentiations of the inputs.

Calculates ``log(exp(x1) + exp(x2))``. This function is useful in
statistics where the calculated probabilities of events may be so small
as to exceed the range of normal floating point numbers.  In such cases
the logarithm of the calculated probability is stored. This function
allows adding probabilities stored in such a fashion.

Parameters
----------
x1, x2 : array_like
    Input values.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
result : ndarray
    Logarithm of ``exp(x1) + exp(x2)``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logaddexp2: Logarithm of the sum of exponentiations of inputs in base 2.

Examples
--------
>>> import numpy as np
>>> prob1 = np.log(1e-50)
>>> prob2 = np.log(2.5e-50)
>>> prob12 = np.logaddexp(prob1, prob2)
>>> prob12
-113.87649168120691
>>> np.exp(prob12)
3.5000000000000057e-50
"""
logaddexp2: numpy.ufunc  # value = <ufunc 'logaddexp2'>
"""
logaddexp2(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Logarithm of the sum of exponentiations of the inputs in base-2.

Calculates ``log2(2**x1 + 2**x2)``. This function is useful in machine
learning when the calculated probabilities of events may be so small as
to exceed the range of normal floating point numbers.  In such cases
the base-2 logarithm of the calculated probability can be used instead.
This function allows adding probabilities stored in such a fashion.

Parameters
----------
x1, x2 : array_like
    Input values.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
result : ndarray
    Base-2 logarithm of ``2**x1 + 2**x2``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logaddexp: Logarithm of the sum of exponentiations of the inputs.

Examples
--------
>>> import numpy as np
>>> prob1 = np.log2(1e-50)
>>> prob2 = np.log2(2.5e-50)
>>> prob12 = np.logaddexp2(prob1, prob2)
>>> prob1, prob2, prob12
(-166.09640474436813, -164.77447664948076, -164.28904982231052)
>>> 2**prob12
3.4999999999999914e-50
"""
logical_and: numpy.ufunc  # value = <ufunc 'logical_and'>
"""
logical_and(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the truth value of x1 AND x2 element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or bool
    Boolean result of the logical AND operation applied to the elements
    of `x1` and `x2`; the shape is determined by broadcasting.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_or, logical_not, logical_xor
bitwise_and

Examples
--------
>>> import numpy as np
>>> np.logical_and(True, False)
False
>>> np.logical_and([True, False], [False, False])
array([False, False])

>>> x = np.arange(5)
>>> np.logical_and(x>1, x<4)
array([False, False,  True,  True, False])


The ``&`` operator can be used as a shorthand for ``np.logical_and`` on
boolean ndarrays.

>>> a = np.array([True, False])
>>> b = np.array([False, False])
>>> a & b
array([False, False])
"""
logical_not: numpy.ufunc  # value = <ufunc 'logical_not'>
"""
logical_not(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the truth value of NOT x element-wise.

Parameters
----------
x : array_like
    Logical NOT is applied to the elements of `x`.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : bool or ndarray of bool
    Boolean result with the same shape as `x` of the NOT operation
    on elements of `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
logical_and, logical_or, logical_xor

Examples
--------
>>> import numpy as np
>>> np.logical_not(3)
False
>>> np.logical_not([True, False, 0, 1])
array([False,  True,  True, False])

>>> x = np.arange(5)
>>> np.logical_not(x<3)
array([False, False, False,  True,  True])
"""
logical_or: numpy.ufunc  # value = <ufunc 'logical_or'>
"""
logical_or(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the truth value of x1 OR x2 element-wise.

Parameters
----------
x1, x2 : array_like
    Logical OR is applied to the elements of `x1` and `x2`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or bool
    Boolean result of the logical OR operation applied to the elements
    of `x1` and `x2`; the shape is determined by broadcasting.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_and, logical_not, logical_xor
bitwise_or

Examples
--------
>>> import numpy as np
>>> np.logical_or(True, False)
True
>>> np.logical_or([True, False], [False, False])
array([ True, False])

>>> x = np.arange(5)
>>> np.logical_or(x < 1, x > 3)
array([ True, False, False, False,  True])

The ``|`` operator can be used as a shorthand for ``np.logical_or`` on
boolean ndarrays.

>>> a = np.array([True, False])
>>> b = np.array([False, False])
>>> a | b
array([ True, False])
"""
logical_xor: numpy.ufunc  # value = <ufunc 'logical_xor'>
"""
logical_xor(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute the truth value of x1 XOR x2, element-wise.

Parameters
----------
x1, x2 : array_like
    Logical XOR is applied to the elements of `x1` and `x2`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : bool or ndarray of bool
    Boolean result of the logical XOR operation applied to the elements
    of `x1` and `x2`; the shape is determined by broadcasting.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
logical_and, logical_or, logical_not, bitwise_xor

Examples
--------
>>> import numpy as np
>>> np.logical_xor(True, False)
True
>>> np.logical_xor([True, True, False, False], [True, False, True, False])
array([False,  True,  True, False])

>>> x = np.arange(5)
>>> np.logical_xor(x < 1, x > 3)
array([ True, False, False, False,  True])

Simple example showing support of broadcasting

>>> np.logical_xor(0, np.eye(2))
array([[ True, False],
       [False,  True]])
"""
matmul: numpy.ufunc  # value = <ufunc 'matmul'>
"""
matmul(x1, x2, /, out=None, *, casting='same_kind', order='K', dtype=None, subok=True[, signature, axes, axis])

Matrix product of two arrays.

Parameters
----------
x1, x2 : array_like
    Input arrays, scalars not allowed.
out : ndarray, optional
    A location into which the result is stored. If provided, it must have
    a shape that matches the signature `(n,k),(k,m)->(n,m)`. If not
    provided or None, a freshly-allocated array is returned.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The matrix product of the inputs.
    This is a scalar only when both x1, x2 are 1-d vectors.

Raises
------
ValueError
    If the last dimension of `x1` is not the same size as
    the second-to-last dimension of `x2`.

    If a scalar value is passed in.

See Also
--------
vecdot : Complex-conjugating dot product for stacks of vectors.
matvec : Matrix-vector product for stacks of matrices and vectors.
vecmat : Vector-matrix product for stacks of vectors and matrices.
tensordot : Sum products over arbitrary axes.
einsum : Einstein summation convention.
dot : alternative matrix product with different broadcasting rules.

Notes
-----
The behavior depends on the arguments in the following way.

- If both arguments are 2-D they are multiplied like conventional
  matrices.
- If either argument is N-D, N > 2, it is treated as a stack of
  matrices residing in the last two indexes and broadcast accordingly.
- If the first argument is 1-D, it is promoted to a matrix by
  prepending a 1 to its dimensions. After matrix multiplication
  the prepended 1 is removed. (For stacks of vectors, use ``vecmat``.)
- If the second argument is 1-D, it is promoted to a matrix by
  appending a 1 to its dimensions. After matrix multiplication
  the appended 1 is removed. (For stacks of vectors, use ``matvec``.)

``matmul`` differs from ``dot`` in two important ways:

- Multiplication by scalars is not allowed, use ``*`` instead.
- Stacks of matrices are broadcast together as if the matrices
  were elements, respecting the signature ``(n,k),(k,m)->(n,m)``:

  >>> a = np.ones([9, 5, 7, 4])
  >>> c = np.ones([9, 5, 4, 3])
  >>> np.dot(a, c).shape
  (9, 5, 7, 9, 5, 3)
  >>> np.matmul(a, c).shape
  (9, 5, 7, 3)
  >>> # n is 7, k is 4, m is 3

The matmul function implements the semantics of the ``@`` operator
defined in :pep:`465`.

It uses an optimized BLAS library when possible (see `numpy.linalg`).

Examples
--------
For 2-D arrays it is the matrix product:

>>> import numpy as np
>>> a = np.array([[1, 0],
...               [0, 1]])
>>> b = np.array([[4, 1],
...               [2, 2]])
>>> np.matmul(a, b)
array([[4, 1],
       [2, 2]])

For 2-D mixed with 1-D, the result is the usual.

>>> a = np.array([[1, 0],
...               [0, 1]])
>>> b = np.array([1, 2])
>>> np.matmul(a, b)
array([1, 2])
>>> np.matmul(b, a)
array([1, 2])


Broadcasting is conventional for stacks of arrays

>>> a = np.arange(2 * 2 * 4).reshape((2, 2, 4))
>>> b = np.arange(2 * 2 * 4).reshape((2, 4, 2))
>>> np.matmul(a,b).shape
(2, 2, 2)
>>> np.matmul(a, b)[0, 1, 1]
98
>>> sum(a[0, 1, :] * b[0 , :, 1])
98

Vector, vector returns the scalar inner product, but neither argument
is complex-conjugated:

>>> np.matmul([2j, 3j], [2j, 3j])
(-13+0j)

Scalar multiplication raises an error.

>>> np.matmul([1,2], 3)
Traceback (most recent call last):
...
ValueError: matmul: Input operand 1 does not have enough dimensions ...

The ``@`` operator can be used as a shorthand for ``np.matmul`` on
ndarrays.

>>> x1 = np.array([2j, 3j])
>>> x2 = np.array([2j, 3j])
>>> x1 @ x2
(-13+0j)
"""
matvec: numpy.ufunc  # value = <ufunc 'matvec'>
r"""
matvec(x1, x2, /, out=None, *, casting='same_kind', order='K', dtype=None, subok=True[, signature, axes, axis])

Matrix-vector dot product of two arrays.

Given a matrix (or stack of matrices) :math:`\mathbf{A}` in ``x1`` and
a vector (or stack of vectors) :math:`\mathbf{v}` in ``x2``, the
matrix-vector product is defined as:

.. math::
   \mathbf{A} \cdot \mathbf{v} = \sum_{j=0}^{n-1} A_{ij} v_j

where the sum is over the last dimensions in ``x1`` and ``x2``
(unless ``axes`` is specified).  (For a matrix-vector product with the
vector conjugated, use ``np.vecmat(x2, x1.mT)``.)

.. versionadded:: 2.2.0

Parameters
----------
x1, x2 : array_like
    Input arrays, scalars not allowed.
out : ndarray, optional
    A location into which the result is stored. If provided, it must have
    the broadcasted shape of ``x1`` and ``x2`` with the summation axis
    removed. If not provided or None, a freshly-allocated array is used.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The matrix-vector product of the inputs.

Raises
------
ValueError
    If the last dimensions of ``x1`` and ``x2`` are not the same size.

    If a scalar value is passed in.

See Also
--------
vecdot : Vector-vector product.
vecmat : Vector-matrix product.
matmul : Matrix-matrix product.
einsum : Einstein summation convention.

Examples
--------
Rotate a set of vectors from Y to X along Z.

>>> a = np.array([[0., 1., 0.],
...               [-1., 0., 0.],
...               [0., 0., 1.]])
>>> v = np.array([[1., 0., 0.],
...               [0., 1., 0.],
...               [0., 0., 1.],
...               [0., 6., 8.]])
>>> np.matvec(a, v)
array([[ 0., -1.,  0.],
       [ 1.,  0.,  0.],
       [ 0.,  0.,  1.],
       [ 6.,  0.,  8.]])
"""
maximum: numpy.ufunc  # value = <ufunc 'maximum'>
"""
maximum(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise maximum of array elements.

Compare two arrays and return a new array containing the element-wise
maxima. If one of the elements being compared is a NaN, then that
element is returned. If both elements are NaNs then the first is
returned. The latter distinction is important for complex NaNs, which
are defined as at least one of the real or imaginary parts being a NaN.
The net effect is that NaNs are propagated.

Parameters
----------
x1, x2 : array_like
    The arrays holding the elements to be compared.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The maximum of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
minimum :
    Element-wise minimum of two arrays, propagates NaNs.
fmax :
    Element-wise maximum of two arrays, ignores NaNs.
amax :
    The maximum value of an array along a given axis, propagates NaNs.
nanmax :
    The maximum value of an array along a given axis, ignores NaNs.

fmin, amin, nanmin

Notes
-----
The maximum is equivalent to ``np.where(x1 >= x2, x1, x2)`` when
neither x1 nor x2 are nans, but it is faster and does proper
broadcasting.

Examples
--------
>>> import numpy as np
>>> np.maximum([2, 3, 4], [1, 5, 2])
array([2, 5, 4])

>>> np.maximum(np.eye(2), [0.5, 2]) # broadcasting
array([[ 1. ,  2. ],
       [ 0.5,  2. ]])

>>> np.maximum([np.nan, 0, np.nan], [0, np.nan, np.nan])
array([nan, nan, nan])
>>> np.maximum(np.inf, 1)
inf
"""
mgrid: numpy.lib._index_tricks_impl.MGridClass  # value = <numpy.lib._index_tricks_impl.MGridClass object>
minimum: numpy.ufunc  # value = <ufunc 'minimum'>
"""
minimum(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Element-wise minimum of array elements.

Compare two arrays and return a new array containing the element-wise
minima. If one of the elements being compared is a NaN, then that
element is returned. If both elements are NaNs then the first is
returned. The latter distinction is important for complex NaNs, which
are defined as at least one of the real or imaginary parts being a NaN.
The net effect is that NaNs are propagated.

Parameters
----------
x1, x2 : array_like
    The arrays holding the elements to be compared.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The minimum of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
maximum :
    Element-wise maximum of two arrays, propagates NaNs.
fmin :
    Element-wise minimum of two arrays, ignores NaNs.
amin :
    The minimum value of an array along a given axis, propagates NaNs.
nanmin :
    The minimum value of an array along a given axis, ignores NaNs.

fmax, amax, nanmax

Notes
-----
The minimum is equivalent to ``np.where(x1 <= x2, x1, x2)`` when
neither x1 nor x2 are NaNs, but it is faster and does proper
broadcasting.

Examples
--------
>>> import numpy as np
>>> np.minimum([2, 3, 4], [1, 5, 2])
array([1, 3, 2])

>>> np.minimum(np.eye(2), [0.5, 2]) # broadcasting
array([[ 0.5,  0. ],
       [ 0. ,  1. ]])

>>> np.minimum([np.nan, 0, np.nan],[0, np.nan, np.nan])
array([nan, nan, nan])
>>> np.minimum(-np.inf, 1)
-inf
"""
mod: numpy.ufunc  # value = <ufunc 'remainder'>
"""
remainder(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns the element-wise remainder of division.

Computes the remainder complementary to the `floor_divide` function.  It is
equivalent to the Python modulus operator ``x1 % x2`` and has the same sign
as the divisor `x2`. The MATLAB function equivalent to ``np.remainder``
is ``mod``.

.. warning::

    This should not be confused with:

    * Python's `math.remainder` and C's ``remainder``, which
      compute the IEEE remainder, which are the complement to
      ``round(x1 / x2)``.
    * The MATLAB ``rem`` function and or the C ``%`` operator which is the
      complement to ``int(x1 / x2)``.

Parameters
----------
x1 : array_like
    Dividend array.
x2 : array_like
    Divisor array.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The element-wise remainder of the quotient ``floor_divide(x1, x2)``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
floor_divide : Equivalent of Python ``//`` operator.
divmod : Simultaneous floor division and remainder.
fmod : Equivalent of the MATLAB ``rem`` function.
divide, floor

Notes
-----
Returns 0 when `x2` is 0 and both `x1` and `x2` are (arrays of)
integers.
``mod`` is an alias of ``remainder``.

Examples
--------
>>> import numpy as np
>>> np.remainder([4, 7], [2, 3])
array([0, 1])
>>> np.remainder(np.arange(7), 5)
array([0, 1, 2, 3, 4, 0, 1])

The ``%`` operator can be used as a shorthand for ``np.remainder`` on
ndarrays.

>>> x1 = np.arange(7)
>>> x1 % 5
array([0, 1, 2, 3, 4, 0, 1])
"""
modf: numpy.ufunc  # value = <ufunc 'modf'>
"""
modf(x[, out1, out2], / [, out=(None, None)], *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the fractional and integral parts of an array, element-wise.

The fractional and integral parts are negative if the given number is
negative.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y1 : ndarray
    Fractional part of `x`.
    This is a scalar if `x` is a scalar.
y2 : ndarray
    Integral part of `x`.
    This is a scalar if `x` is a scalar.

Notes
-----
For integer input the return values are floats.

See Also
--------
divmod : ``divmod(x, 1)`` is equivalent to ``modf`` with the return values
         switched, except it always has a positive remainder.

Examples
--------
>>> import numpy as np
>>> np.modf([0, 3.5])
(array([ 0. ,  0.5]), array([ 0.,  3.]))
>>> np.modf(-0.5)
(-0.5, -0)
"""
multiply: numpy.ufunc  # value = <ufunc 'multiply'>
"""
multiply(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Multiply arguments element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays to be multiplied.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The product of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

Notes
-----
Equivalent to `x1` * `x2` in terms of array broadcasting.

Examples
--------
>>> import numpy as np
>>> np.multiply(2.0, 4.0)
8.0

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> np.multiply(x1, x2)
array([[  0.,   1.,   4.],
       [  0.,   4.,  10.],
       [  0.,   7.,  16.]])

The ``*`` operator can be used as a shorthand for ``np.multiply`` on
ndarrays.

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> x1 * x2
array([[  0.,   1.,   4.],
       [  0.,   4.,  10.],
       [  0.,   7.,  16.]])
"""
nan: float  # value = nan
negative: numpy.ufunc  # value = <ufunc 'negative'>
"""
negative(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Numerical negative, element-wise.

Parameters
----------
x : array_like or scalar
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    Returned array or scalar: `y = -x`.
    This is a scalar if `x` is a scalar.

Examples
--------
>>> import numpy as np
>>> np.negative([1.,-1.])
array([-1.,  1.])

The unary ``-`` operator can be used as a shorthand for ``np.negative`` on
ndarrays.

>>> x1 = np.array(([1., -1.]))
>>> -x1
array([-1.,  1.])
"""
newaxis = None
nextafter: numpy.ufunc  # value = <ufunc 'nextafter'>
"""
nextafter(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the next floating-point value after x1 towards x2, element-wise.

Parameters
----------
x1 : array_like
    Values to find the next representable value of.
x2 : array_like
    The direction where to look for the next representable value of `x1`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    The next representable values of `x1` in the direction of `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

Examples
--------
>>> import numpy as np
>>> eps = np.finfo(np.float64).eps
>>> np.nextafter(1, 2) == eps + 1
True
>>> np.nextafter([1, 2], [2, 1]) == [eps + 1, 2 - eps]
array([ True,  True])
"""
not_equal: numpy.ufunc  # value = <ufunc 'not_equal'>
"""
not_equal(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return (x1 != x2) element-wise.

Parameters
----------
x1, x2 : array_like
    Input arrays.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array, element-wise comparison of `x1` and `x2`.
    Typically of type bool, unless ``dtype=object`` is passed.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
equal, greater, greater_equal, less, less_equal

Examples
--------
>>> import numpy as np
>>> np.not_equal([1.,2.], [1., 3.])
array([False,  True])
>>> np.not_equal([1, 2], [[1, 3],[1, 4]])
array([[False,  True],
       [False,  True]])

The ``!=`` operator can be used as a shorthand for ``np.not_equal`` on
ndarrays.

>>> a = np.array([1., 2.])
>>> b = np.array([1., 3.])
>>> a != b
array([False,  True])
"""
ogrid: numpy.lib._index_tricks_impl.OGridClass  # value = <numpy.lib._index_tricks_impl.OGridClass object>
pi: float = 3.141592653589793
positive: numpy.ufunc  # value = <ufunc 'positive'>
"""
positive(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Numerical positive, element-wise.

Parameters
----------
x : array_like or scalar
    Input array.

Returns
-------
y : ndarray or scalar
    Returned array or scalar: `y = +x`.
    This is a scalar if `x` is a scalar.

Notes
-----
Equivalent to `x.copy()`, but only defined for types that support
arithmetic.

Examples
--------
>>> import numpy as np

>>> x1 = np.array(([1., -1.]))
>>> np.positive(x1)
array([ 1., -1.])

The unary ``+`` operator can be used as a shorthand for ``np.positive`` on
ndarrays.

>>> x1 = np.array(([1., -1.]))
>>> +x1
array([ 1., -1.])
"""
pow: numpy.ufunc  # value = <ufunc 'power'>
"""
power(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

First array elements raised to powers from second array, element-wise.

Raise each base in `x1` to the positionally-corresponding power in
`x2`.  `x1` and `x2` must be broadcastable to the same shape.

An integer type raised to a negative integer power will raise a
``ValueError``.

Negative values raised to a non-integral value will return ``nan``.
To get complex results, cast the input to complex, or specify the
``dtype`` to be ``complex`` (see the example below).

Parameters
----------
x1 : array_like
    The bases.
x2 : array_like
    The exponents.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The bases in `x1` raised to the exponents in `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
float_power : power function that promotes integers to float

Examples
--------
>>> import numpy as np

Cube each element in an array.

>>> x1 = np.arange(6)
>>> x1
[0, 1, 2, 3, 4, 5]
>>> np.power(x1, 3)
array([  0,   1,   8,  27,  64, 125])

Raise the bases to different exponents.

>>> x2 = [1.0, 2.0, 3.0, 3.0, 2.0, 1.0]
>>> np.power(x1, x2)
array([  0.,   1.,   8.,  27.,  16.,   5.])

The effect of broadcasting.

>>> x2 = np.array([[1, 2, 3, 3, 2, 1], [1, 2, 3, 3, 2, 1]])
>>> x2
array([[1, 2, 3, 3, 2, 1],
       [1, 2, 3, 3, 2, 1]])
>>> np.power(x1, x2)
array([[ 0,  1,  8, 27, 16,  5],
       [ 0,  1,  8, 27, 16,  5]])

The ``**`` operator can be used as a shorthand for ``np.power`` on
ndarrays.

>>> x2 = np.array([1, 2, 3, 3, 2, 1])
>>> x1 = np.arange(6)
>>> x1 ** x2
array([ 0,  1,  8, 27, 16,  5])

Negative values raised to a non-integral value will result in ``nan``
(and a warning will be generated).

>>> x3 = np.array([-1.0, -4.0])
>>> with np.errstate(invalid='ignore'):
...     p = np.power(x3, 1.5)
...
>>> p
array([nan, nan])

To get complex results, give the argument ``dtype=complex``.

>>> np.power(x3, 1.5, dtype=complex)
array([-1.83697020e-16-1.j, -1.46957616e-15-8.j])
"""
power: numpy.ufunc  # value = <ufunc 'power'>
"""
power(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

First array elements raised to powers from second array, element-wise.

Raise each base in `x1` to the positionally-corresponding power in
`x2`.  `x1` and `x2` must be broadcastable to the same shape.

An integer type raised to a negative integer power will raise a
``ValueError``.

Negative values raised to a non-integral value will return ``nan``.
To get complex results, cast the input to complex, or specify the
``dtype`` to be ``complex`` (see the example below).

Parameters
----------
x1 : array_like
    The bases.
x2 : array_like
    The exponents.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The bases in `x1` raised to the exponents in `x2`.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
float_power : power function that promotes integers to float

Examples
--------
>>> import numpy as np

Cube each element in an array.

>>> x1 = np.arange(6)
>>> x1
[0, 1, 2, 3, 4, 5]
>>> np.power(x1, 3)
array([  0,   1,   8,  27,  64, 125])

Raise the bases to different exponents.

>>> x2 = [1.0, 2.0, 3.0, 3.0, 2.0, 1.0]
>>> np.power(x1, x2)
array([  0.,   1.,   8.,  27.,  16.,   5.])

The effect of broadcasting.

>>> x2 = np.array([[1, 2, 3, 3, 2, 1], [1, 2, 3, 3, 2, 1]])
>>> x2
array([[1, 2, 3, 3, 2, 1],
       [1, 2, 3, 3, 2, 1]])
>>> np.power(x1, x2)
array([[ 0,  1,  8, 27, 16,  5],
       [ 0,  1,  8, 27, 16,  5]])

The ``**`` operator can be used as a shorthand for ``np.power`` on
ndarrays.

>>> x2 = np.array([1, 2, 3, 3, 2, 1])
>>> x1 = np.arange(6)
>>> x1 ** x2
array([ 0,  1,  8, 27, 16,  5])

Negative values raised to a non-integral value will result in ``nan``
(and a warning will be generated).

>>> x3 = np.array([-1.0, -4.0])
>>> with np.errstate(invalid='ignore'):
...     p = np.power(x3, 1.5)
...
>>> p
array([nan, nan])

To get complex results, give the argument ``dtype=complex``.

>>> np.power(x3, 1.5, dtype=complex)
array([-1.83697020e-16-1.j, -1.46957616e-15-8.j])
"""
r_: numpy.lib._index_tricks_impl.RClass  # value = <numpy.lib._index_tricks_impl.RClass object>
rad2deg: numpy.ufunc  # value = <ufunc 'rad2deg'>
"""
rad2deg(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Convert angles from radians to degrees.

Parameters
----------
x : array_like
    Angle in radians.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding angle in degrees.
    This is a scalar if `x` is a scalar.

See Also
--------
deg2rad : Convert angles from degrees to radians.
unwrap : Remove large jumps in angle by wrapping.

Notes
-----
rad2deg(x) is ``180 * x / pi``.

Examples
--------
>>> import numpy as np
>>> np.rad2deg(np.pi/2)
90.0
"""
radians: numpy.ufunc  # value = <ufunc 'radians'>
"""
radians(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Convert angles from degrees to radians.

Parameters
----------
x : array_like
    Input array in degrees.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding radian values.
    This is a scalar if `x` is a scalar.

See Also
--------
deg2rad : equivalent function

Examples
--------
>>> import numpy as np

Convert a degree array to radians

>>> deg = np.arange(12.) * 30.
>>> np.radians(deg)
array([ 0.        ,  0.52359878,  1.04719755,  1.57079633,  2.0943951 ,
        2.61799388,  3.14159265,  3.66519143,  4.1887902 ,  4.71238898,
        5.23598776,  5.75958653])

>>> out = np.zeros((deg.shape))
>>> ret = np.radians(deg, out)
>>> ret is out
True
"""
reciprocal: numpy.ufunc  # value = <ufunc 'reciprocal'>
"""
reciprocal(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the reciprocal of the argument, element-wise.

Calculates ``1/x``.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    Return array.
    This is a scalar if `x` is a scalar.

Notes
-----
.. note::
    This function is not designed to work with integers.

For integer arguments with absolute value larger than 1 the result is
always zero because of the way Python handles integer division.  For
integer zero the result is an overflow.

Examples
--------
>>> import numpy as np
>>> np.reciprocal(2.)
0.5
>>> np.reciprocal([1, 2., 3.33])
array([ 1.       ,  0.5      ,  0.3003003])
"""
remainder: numpy.ufunc  # value = <ufunc 'remainder'>
"""
remainder(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns the element-wise remainder of division.

Computes the remainder complementary to the `floor_divide` function.  It is
equivalent to the Python modulus operator ``x1 % x2`` and has the same sign
as the divisor `x2`. The MATLAB function equivalent to ``np.remainder``
is ``mod``.

.. warning::

    This should not be confused with:

    * Python's `math.remainder` and C's ``remainder``, which
      compute the IEEE remainder, which are the complement to
      ``round(x1 / x2)``.
    * The MATLAB ``rem`` function and or the C ``%`` operator which is the
      complement to ``int(x1 / x2)``.

Parameters
----------
x1 : array_like
    Dividend array.
x2 : array_like
    Divisor array.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The element-wise remainder of the quotient ``floor_divide(x1, x2)``.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
floor_divide : Equivalent of Python ``//`` operator.
divmod : Simultaneous floor division and remainder.
fmod : Equivalent of the MATLAB ``rem`` function.
divide, floor

Notes
-----
Returns 0 when `x2` is 0 and both `x1` and `x2` are (arrays of)
integers.
``mod`` is an alias of ``remainder``.

Examples
--------
>>> import numpy as np
>>> np.remainder([4, 7], [2, 3])
array([0, 1])
>>> np.remainder(np.arange(7), 5)
array([0, 1, 2, 3, 4, 0, 1])

The ``%`` operator can be used as a shorthand for ``np.remainder`` on
ndarrays.

>>> x1 = np.arange(7)
>>> x1 % 5
array([0, 1, 2, 3, 4, 0, 1])
"""
right_shift: numpy.ufunc  # value = <ufunc 'right_shift'>
"""
right_shift(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Shift the bits of an integer to the right.

Bits are shifted to the right `x2`.  Because the internal
representation of numbers is in binary format, this operation is
equivalent to dividing `x1` by ``2**x2``.

Parameters
----------
x1 : array_like, int
    Input values.
x2 : array_like, int
    Number of bits to remove at the right of `x1`.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray, int
    Return `x1` with bits shifted `x2` times to the right.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
left_shift : Shift the bits of an integer to the left.
binary_repr : Return the binary representation of the input number
    as a string.

Examples
--------
>>> import numpy as np
>>> np.binary_repr(10)
'1010'
>>> np.right_shift(10, 1)
5
>>> np.binary_repr(5)
'101'

>>> np.right_shift(10, [1,2,3])
array([5, 2, 1])

The ``>>`` operator can be used as a shorthand for ``np.right_shift`` on
ndarrays.

>>> x1 = 10
>>> x2 = np.array([1,2,3])
>>> x1 >> x2
array([5, 2, 1])
"""
rint: numpy.ufunc  # value = <ufunc 'rint'>
"""
rint(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Round elements of the array to the nearest integer.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Output array is same shape and type as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
fix, ceil, floor, trunc

Notes
-----
For values exactly halfway between rounded decimal values, NumPy
rounds to the nearest even value. Thus 1.5 and 2.5 round to 2.0,
-0.5 and 0.5 round to 0.0, etc.

Examples
--------
>>> import numpy as np
>>> a = np.array([-1.7, -1.5, -0.2, 0.2, 1.5, 1.7, 2.0])
>>> np.rint(a)
array([-2., -2., -0.,  0.,  2.,  2.,  2.])
"""
s_: numpy.lib._index_tricks_impl.IndexExpression  # value = <numpy.lib._index_tricks_impl.IndexExpression object>
sctypeDict: dict = {'bool': numpy.bool, 'float16': numpy.float16, 'float32': numpy.float32, 'float64': numpy.float64, 'longdouble': numpy.longdouble, 'complex64': numpy.complex64, 'complex128': numpy.complex128, 'clongdouble': numpy.clongdouble, 'bytes_': numpy.bytes_, 'str_': numpy.str_, 'void': numpy.void, 'object_': numpy.object_, 'datetime64': numpy.datetime64, 'timedelta64': numpy.timedelta64, 'int8': numpy.int8, 'byte': numpy.int8, 'uint8': numpy.uint8, 'ubyte': numpy.uint8, 'int16': numpy.int16, 'short': numpy.int16, 'uint16': numpy.uint16, 'ushort': numpy.uint16, 'int32': numpy.int32, 'intc': numpy.int32, 'uint32': numpy.uint32, 'uintc': numpy.uint32, 'int64': numpy.int64, 'long': numpy.int64, 'uint64': numpy.uint64, 'ulong': numpy.uint64, 'longlong': numpy.longlong, 'ulonglong': numpy.ulonglong, 'intp': numpy.int64, 'uintp': numpy.uint64, 'double': numpy.float64, 'cdouble': numpy.complex128, 'single': numpy.float32, 'csingle': numpy.complex64, 'half': numpy.float16, 'bool_': numpy.bool, 'int_': numpy.int64, 'uint': numpy.uint64, 'float': numpy.float64, 'complex': numpy.complex128, 'object': numpy.object_, 'bytes': numpy.bytes_, 'a': numpy.bytes_, 'int': numpy.int64, 'str': numpy.str_, 'unicode': numpy.str_}
sign: numpy.ufunc  # value = <ufunc 'sign'>
r"""
sign(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns an element-wise indication of the sign of a number.

The `sign` function returns ``-1 if x < 0, 0 if x==0, 1 if x > 0``.  nan
is returned for nan inputs.

For complex inputs, the `sign` function returns ``x / abs(x)``, the
generalization of the above (and ``0 if x==0``).

.. versionchanged:: 2.0.0
    Definition of complex sign changed to follow the Array API standard.

Parameters
----------
x : array_like
    Input values.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The sign of `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
signbit
copysign

Notes
-----
There is more than one definition of sign in common use for complex
numbers.  The definition used here, :math:`x/|x|`, is the more common
and useful one, but is different from the one used in numpy prior to
version 2.0, :math:`x/\sqrt{x*x}`, which is equivalent to
``sign(x.real) + 0j if x.real != 0 else sign(x.imag) + 0j``.

Examples
--------
>>> import numpy as np
>>> np.sign([-5., 4.5])
array([-1.,  1.])
>>> np.sign(0)
0
>>> np.sign([3-4j, 8j])
array([0.6-0.8j, 0. +1.j ])
"""
signbit: numpy.ufunc  # value = <ufunc 'signbit'>
"""
signbit(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Returns element-wise True where signbit is set (less than zero).

Parameters
----------
x : array_like
    The input value(s).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
result : ndarray of bool
    Output array, or reference to `out` if that was supplied.
    This is a scalar if `x` is a scalar.

See Also
--------
sign
copysign

Examples
--------
>>> import numpy as np
>>> np.signbit(-1.2)
True
>>> np.signbit(np.array([1, -2.3, 2.1]))
array([False,  True, False])
"""
sin: numpy.ufunc  # value = <ufunc 'sin'>
r"""
sin(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Trigonometric sine, element-wise.

Parameters
----------
x : array_like
    Angle, in radians (:math:`2 \pi` rad equals 360 degrees).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : array_like
    The sine of each element of x.
    This is a scalar if `x` is a scalar.

See Also
--------
arcsin, sinh, cos

Notes
-----
The sine is one of the fundamental functions of trigonometry (the
mathematical study of triangles).  Consider a circle of radius 1
centered on the origin.  A ray comes in from the :math:`+x` axis, makes
an angle at the origin (measured counter-clockwise from that axis), and
departs from the origin.  The :math:`y` coordinate of the outgoing
ray's intersection with the unit circle is the sine of that angle.  It
ranges from -1 for :math:`x=3\pi / 2` to +1 for :math:`\pi / 2.`  The
function has zeroes where the angle is a multiple of :math:`\pi`.
Sines of angles between :math:`\pi` and :math:`2\pi` are negative.
The numerous properties of the sine and related functions are included
in any standard trigonometry text.

Examples
--------
>>> import numpy as np

Print sine of one angle:

>>> np.sin(np.pi/2.)
1.0

Print sines of an array of angles given in degrees:

>>> np.sin(np.array((0., 30., 45., 60., 90.)) * np.pi / 180. )
array([ 0.        ,  0.5       ,  0.70710678,  0.8660254 ,  1.        ])

Plot the sine function:

>>> import matplotlib.pylab as plt
>>> x = np.linspace(-np.pi, np.pi, 201)
>>> plt.plot(x, np.sin(x))
>>> plt.xlabel('Angle [rad]')
>>> plt.ylabel('sin(x)')
>>> plt.axis('tight')
>>> plt.show()
"""
sinh: numpy.ufunc  # value = <ufunc 'sinh'>
"""
sinh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Hyperbolic sine, element-wise.

Equivalent to ``1/2 * (np.exp(x) - np.exp(-x))`` or
``-1j * np.sin(1j*x)``.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding hyperbolic sine values.
    This is a scalar if `x` is a scalar.

Notes
-----
If `out` is provided, the function writes the result into it,
and returns a reference to `out`.  (See Examples)

References
----------
M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions.
New York, NY: Dover, 1972, pg. 83.

Examples
--------
>>> import numpy as np
>>> np.sinh(0)
0.0
>>> np.sinh(np.pi*1j/2)
1j
>>> np.sinh(np.pi*1j) # (exact value is 0)
1.2246063538223773e-016j
>>> # Discrepancy due to vagaries of floating point arithmetic.

>>> # Example of providing the optional output parameter
>>> out1 = np.array([0], dtype='d')
>>> out2 = np.sinh([0.1], out1)
>>> out2 is out1
True

>>> # Example of ValueError due to provision of shape mis-matched `out`
>>> np.sinh(np.zeros((3,3)),np.zeros((2,2)))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: operands could not be broadcast together with shapes (3,3) (2,2)
"""
spacing: numpy.ufunc  # value = <ufunc 'spacing'>
"""
spacing(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the distance between x and the nearest adjacent number.

Parameters
----------
x : array_like
    Values to find the spacing of.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    The spacing of values of `x`.
    This is a scalar if `x` is a scalar.

Notes
-----
It can be considered as a generalization of EPS:
``spacing(np.float64(1)) == np.finfo(np.float64).eps``, and there
should not be any representable number between ``x + spacing(x)`` and
x for any finite x.

Spacing of +- inf and NaN is NaN.

Examples
--------
>>> import numpy as np
>>> np.spacing(1) == np.finfo(np.float64).eps
True
"""
sqrt: numpy.ufunc  # value = <ufunc 'sqrt'>
"""
sqrt(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the non-negative square-root of an array, element-wise.

Parameters
----------
x : array_like
    The values whose square-roots are required.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    An array of the same shape as `x`, containing the positive
    square-root of each element in `x`.  If any element in `x` is
    complex, a complex array is returned (and the square-roots of
    negative reals are calculated).  If all of the elements in `x`
    are real, so is `y`, with negative elements returning ``nan``.
    If `out` was provided, `y` is a reference to it.
    This is a scalar if `x` is a scalar.

See Also
--------
emath.sqrt
    A version which returns complex numbers when given negative reals.
    Note that 0.0 and -0.0 are handled differently for complex inputs.

Notes
-----
*sqrt* has--consistent with common convention--as its branch cut the
real "interval" [`-inf`, 0), and is continuous from above on it.
A branch cut is a curve in the complex plane across which a given
complex function fails to be continuous.

Examples
--------
>>> import numpy as np
>>> np.sqrt([1,4,9])
array([ 1.,  2.,  3.])

>>> np.sqrt([4, -1, -3+4J])
array([ 2.+0.j,  0.+1.j,  1.+2.j])

>>> np.sqrt([4, -1, np.inf])
array([ 2., nan, inf])
"""
square: numpy.ufunc  # value = <ufunc 'square'>
"""
square(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the element-wise square of the input.

Parameters
----------
x : array_like
    Input data.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
out : ndarray or scalar
    Element-wise `x*x`, of the same shape and dtype as `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
numpy.linalg.matrix_power
sqrt
power

Examples
--------
>>> import numpy as np
>>> np.square([-1j, 1])
array([-1.-0.j,  1.+0.j])
"""
subtract: numpy.ufunc  # value = <ufunc 'subtract'>
"""
subtract(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Subtract arguments, element-wise.

Parameters
----------
x1, x2 : array_like
    The arrays to be subtracted from each other.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The difference of `x1` and `x2`, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

Notes
-----
Equivalent to ``x1 - x2`` in terms of array broadcasting.

Examples
--------
>>> import numpy as np
>>> np.subtract(1.0, 4.0)
-3.0

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> np.subtract(x1, x2)
array([[ 0.,  0.,  0.],
       [ 3.,  3.,  3.],
       [ 6.,  6.,  6.]])

The ``-`` operator can be used as a shorthand for ``np.subtract`` on
ndarrays.

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> x1 - x2
array([[0., 0., 0.],
       [3., 3., 3.],
       [6., 6., 6.]])
"""
tan: numpy.ufunc  # value = <ufunc 'tan'>
"""
tan(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute tangent element-wise.

Equivalent to ``np.sin(x)/np.cos(x)`` element-wise.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding tangent values.
    This is a scalar if `x` is a scalar.

Notes
-----
If `out` is provided, the function writes the result into it,
and returns a reference to `out`.  (See Examples)

References
----------
M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions.
New York, NY: Dover, 1972.

Examples
--------
>>> import numpy as np
>>> from math import pi
>>> np.tan(np.array([-pi,pi/2,pi]))
array([  1.22460635e-16,   1.63317787e+16,  -1.22460635e-16])
>>>
>>> # Example of providing the optional output parameter illustrating
>>> # that what is returned is a reference to said parameter
>>> out1 = np.array([0], dtype='d')
>>> out2 = np.cos([0.1], out1)
>>> out2 is out1
True
>>>
>>> # Example of ValueError due to provision of shape mis-matched `out`
>>> np.cos(np.zeros((3,3)),np.zeros((2,2)))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: operands could not be broadcast together with shapes (3,3) (2,2)
"""
tanh: numpy.ufunc  # value = <ufunc 'tanh'>
"""
tanh(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Compute hyperbolic tangent element-wise.

Equivalent to ``np.sinh(x)/np.cosh(x)`` or ``-1j * np.tan(1j*x)``.

Parameters
----------
x : array_like
    Input array.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The corresponding hyperbolic tangent values.
    This is a scalar if `x` is a scalar.

Notes
-----
If `out` is provided, the function writes the result into it,
and returns a reference to `out`.  (See Examples)

References
----------
.. [1] M. Abramowitz and I. A. Stegun, Handbook of Mathematical Functions.
       New York, NY: Dover, 1972, pg. 83.
       https://personal.math.ubc.ca/~cbm/aands/page_83.htm

.. [2] Wikipedia, "Hyperbolic function",
       https://en.wikipedia.org/wiki/Hyperbolic_function

Examples
--------
>>> import numpy as np
>>> np.tanh((0, np.pi*1j, np.pi*1j/2))
array([ 0. +0.00000000e+00j,  0. -1.22460635e-16j,  0. +1.63317787e+16j])

>>> # Example of providing the optional output parameter illustrating
>>> # that what is returned is a reference to said parameter
>>> out1 = np.array([0], dtype='d')
>>> out2 = np.tanh([0.1], out1)
>>> out2 is out1
True

>>> # Example of ValueError due to provision of shape mis-matched `out`
>>> np.tanh(np.zeros((3,3)),np.zeros((2,2)))
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: operands could not be broadcast together with shapes (3,3) (2,2)
"""
test: numpy._pytesttester.PytestTester  # value = <numpy._pytesttester.PytestTester object>
true_divide: numpy.ufunc  # value = <ufunc 'divide'>
"""
divide(x1, x2, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Divide arguments element-wise.

Parameters
----------
x1 : array_like
    Dividend array.
x2 : array_like
    Divisor array.
    If ``x1.shape != x2.shape``, they must be broadcastable to a common
    shape (which becomes the shape of the output).
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The quotient ``x1/x2``, element-wise.
    This is a scalar if both `x1` and `x2` are scalars.

See Also
--------
seterr : Set whether to raise or warn on overflow, underflow and
         division by zero.

Notes
-----
Equivalent to ``x1`` / ``x2`` in terms of array-broadcasting.

The ``true_divide(x1, x2)`` function is an alias for
``divide(x1, x2)``.

Examples
--------
>>> import numpy as np
>>> np.divide(2.0, 4.0)
0.5
>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = np.arange(3.0)
>>> np.divide(x1, x2)
array([[nan, 1. , 1. ],
       [inf, 4. , 2.5],
       [inf, 7. , 4. ]])

The ``/`` operator can be used as a shorthand for ``np.divide`` on
ndarrays.

>>> x1 = np.arange(9.0).reshape((3, 3))
>>> x2 = 2 * np.ones(3)
>>> x1 / x2
array([[0. , 0.5, 1. ],
       [1.5, 2. , 2.5],
       [3. , 3.5, 4. ]])
"""
trunc: numpy.ufunc  # value = <ufunc 'trunc'>
"""
trunc(x, /, out=None, *, where=True, casting='same_kind', order='K', dtype=None, subok=True[, signature])

Return the truncated value of the input, element-wise.

The truncated value of the scalar `x` is the nearest integer `i` which
is closer to zero than `x` is. In short, the fractional part of the
signed number `x` is discarded.

Parameters
----------
x : array_like
    Input data.
out : ndarray, None, or tuple of ndarray and None, optional
    A location into which the result is stored. If provided, it must have
    a shape that the inputs broadcast to. If not provided or None,
    a freshly-allocated array is returned. A tuple (possible only as a
    keyword argument) must have length equal to the number of outputs.
where : array_like, optional
    This condition is broadcast over the input. At locations where the
    condition is True, the `out` array will be set to the ufunc result.
    Elsewhere, the `out` array will retain its original value.
    Note that if an uninitialized `out` array is created via the default
    ``out=None``, locations within it where the condition is False will
    remain uninitialized.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray or scalar
    The truncated value of each element in `x`.
    This is a scalar if `x` is a scalar.

See Also
--------
ceil, floor, rint, fix

Examples
--------
>>> import numpy as np
>>> a = np.array([-1.7, -1.5, -0.2, 0.2, 1.5, 1.7, 2.0])
>>> np.trunc(a)
array([-1., -1., -0.,  0.,  1.,  1.,  2.])
"""
typecodes: dict = {'Character': 'c', 'Integer': 'bhilqnp', 'UnsignedInteger': 'BHILQNP', 'Float': 'efdg', 'Complex': 'FDG', 'AllInteger': 'bBhHiIlLqQnNpP', 'AllFloat': 'efdgFDG', 'Datetime': 'Mm', 'All': '?bhilqnpBHILQNPefdgFDGSUVOMm'}
vecdot: numpy.ufunc  # value = <ufunc 'vecdot'>
r"""
vecdot(x1, x2, /, out=None, *, casting='same_kind', order='K', dtype=None, subok=True[, signature, axes, axis])

Vector dot product of two arrays.

Let :math:`\mathbf{a}` be a vector in `x1` and :math:`\mathbf{b}` be
a corresponding vector in `x2`. The dot product is defined as:

.. math::
   \mathbf{a} \cdot \mathbf{b} = \sum_{i=0}^{n-1} \overline{a_i}b_i

where the sum is over the last dimension (unless `axis` is specified) and
where :math:`\overline{a_i}` denotes the complex conjugate if :math:`a_i`
is complex and the identity otherwise.

.. versionadded:: 2.0.0

Parameters
----------
x1, x2 : array_like
    Input arrays, scalars not allowed.
out : ndarray, optional
    A location into which the result is stored. If provided, it must have
    the broadcasted shape of `x1` and `x2` with the last axis removed.
    If not provided or None, a freshly-allocated array is used.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The vector dot product of the inputs.
    This is a scalar only when both x1, x2 are 1-d vectors.

Raises
------
ValueError
    If the last dimension of `x1` is not the same size as
    the last dimension of `x2`.

    If a scalar value is passed in.

See Also
--------
vdot : same but flattens arguments first
matmul : Matrix-matrix product.
vecmat : Vector-matrix product.
matvec : Matrix-vector product.
einsum : Einstein summation convention.

Examples
--------
>>> import numpy as np

Get the projected size along a given normal for an array of vectors.

>>> v = np.array([[0., 5., 0.], [0., 0., 10.], [0., 6., 8.]])
>>> n = np.array([0., 0.6, 0.8])
>>> np.vecdot(v, n)
array([ 3.,  8., 10.])
"""
vecmat: numpy.ufunc  # value = <ufunc 'vecmat'>
r"""
vecmat(x1, x2, /, out=None, *, casting='same_kind', order='K', dtype=None, subok=True[, signature, axes, axis])

Vector-matrix dot product of two arrays.

Given a vector (or stack of vector) :math:`\mathbf{v}` in ``x1`` and
a matrix (or stack of matrices) :math:`\mathbf{A}` in ``x2``, the
vector-matrix product is defined as:

.. math::
   \mathbf{v} \cdot \mathbf{A} = \sum_{i=0}^{n-1} \overline{v_i}A_{ij}

where the sum is over the last dimension of ``x1`` and the one-but-last
dimensions in ``x2`` (unless `axes` is specified) and where
:math:`\overline{v_i}` denotes the complex conjugate if :math:`v`
is complex and the identity otherwise. (For a non-conjugated vector-matrix
product, use ``np.matvec(x2.mT, x1)``.)

.. versionadded:: 2.2.0

Parameters
----------
x1, x2 : array_like
    Input arrays, scalars not allowed.
out : ndarray, optional
    A location into which the result is stored. If provided, it must have
    the broadcasted shape of ``x1`` and ``x2`` with the summation axis
    removed. If not provided or None, a freshly-allocated array is used.
**kwargs
    For other keyword-only arguments, see the
    :ref:`ufunc docs <ufuncs.kwargs>`.

Returns
-------
y : ndarray
    The vector-matrix product of the inputs.

Raises
------
ValueError
    If the last dimensions of ``x1`` and the one-but-last dimension of
    ``x2`` are not the same size.

    If a scalar value is passed in.

See Also
--------
vecdot : Vector-vector product.
matvec : Matrix-vector product.
matmul : Matrix-matrix product.
einsum : Einstein summation convention.

Examples
--------
Project a vector along X and Y.

>>> v = np.array([0., 4., 2.])
>>> a = np.array([[1., 0., 0.],
...               [0., 1., 0.],
...               [0., 0., 0.]])
>>> np.vecmat(v, a)
array([ 0.,  4., 0.])
"""
