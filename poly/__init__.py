"""Polynomial Interface

Provides access to dense and sparse polynomial representations through a single 
`Poly` class.

"""

from . import sparse
from . import dense
from . import utils

from .sparse import *
from .dense import *
from .utils import *




class Poly:
  def __init__(
      self,
      coeffs: list[int] | dict[int, int] | str | PolyDense | PolySparse,
      poly_ring_mod = None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    self._is_sparse = False
    self._poly_sparse = None
    self._poly_dense = None
    # Check coeffs
    if coeffs.__class__ == PolyDense:
      self._is_sparse = False
      self._poly_dense = coeffs
      self._poly_sparse = None
    elif coeffs.__class__ == PolySparse:
      self._is_sparse = True
      self._poly_sparse = coeffs
      self._poly_dense = None
    elif coeffs.__class__ == dict:
      self._is_sparse = True
      self._poly_sparse = PolySparse(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == list:
      self._is_sparse = False
      self._poly_dense = PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == str:
      # For now return a dense representation.
      try:
        coeffs_dense = poly_parse_string(coeffs, var="x")
      except SyntaxError:
        raise SyntaxError("Invalid polynomial string")
      self._is_sparse = False
      self._poly_dense = PolyDense(
          coeffs_dense,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == int:
      self._is_sparse = False
      self._poly_dense = PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    else:
      raise ValueError("coeffs must be a list, dictionary, or string")

  def _poly(self):
    if self._is_sparse:
      return self._poly_sparse
    else:
      return self._poly_dense

  def __str__(self):
    return self._poly().__str__()

  def __add__(self, other):
    return Poly(self._poly() + other._poly())

  def __sub__(self, other):
    return Poly(self._poly() - other._poly())

  def __mul__(self, other):
    return Poly(self._poly() * other._poly())

  def __pow__(self, n):
    return generalized_fast_pow(self, n, Poly.__mul__)

  def __floordiv__(self, other):
    q, r = self.__divmod__(other)
    return q

  def __mod__(self, other):
    q, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
    q, r = self._poly().__divmod__(other._poly())
    return Poly(q), Poly(r)

  def __deepcopy__(self):
    return Poly(self._poly())

  def is_sparse(self):
    return self._is_sparse

  def poly_ring_mod(self):
    return self._poly()._poly_ring_mod
