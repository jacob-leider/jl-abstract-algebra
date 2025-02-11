"""Polynomial Interface

Provides access to dense and sparse polynomial representations through a single 
`Poly` class.

"""

from .sparse import *
from .dense import *
from .utils import *

class Poly:
  def __init__(
      self,
      coeffs: list[int] | dict[int, int] | str | PolyDense | PolySparse | int,
      poly_ring_mod = None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    # Check coeffs
    if isinstance(coeffs, PolyDense):
      self._rep = coeffs
    elif isinstance(coeffs, PolySparse):
      self._rep = coeffs
    elif isinstance(coeffs, dict):
      self._rep = PolySparse(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif isinstance(coeffs, list):
      self._rep = PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif isinstance(coeffs, str):
      # For now return a dense representation. usually th
      try:
        coeffs_dense = poly_parse_string(coeffs, var="x")
      except SyntaxError:
        raise SyntaxError("Invalid polynomial string")
      self._rep = PolyDense(
          coeffs_dense,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif isinstance(coeffs, int):
      self._rep = PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    else:
      raise ValueError("coeffs must be a list, dictionary, or string")


  def __str__(self, **kwargs):
    return poly_string(self._poly()._coeffs, **kwargs)

  def __add__(self, other):
    if isinstance(other, Poly):
      return Poly(self._poly() + other._poly())
    else:
      return Poly(self._poly() + other)

  def __sub__(self, other):
    if isinstance(other, Poly):
      return Poly(self._poly() - other._poly())
    else:
      return Poly(self._poly() - other)

  def __rmul__(self, other):
    if isinstance(other, Poly):
      return Poly(self._poly() * other._poly())
    else:
      return Poly(self._poly() * other)

  def __mul__(self, other):
    if isinstance(other, Poly):
      return Poly(self._poly() * other._poly())
    else:
      return Poly(self._poly() * other)

  def __pow__(self, n):
    if isinstance(n, int):
      return generalized_fast_pow(self, n, Poly.__mul__)
    else:
      return NotImplemented

  def __floordiv__(self, other):
    q, _ = self.__divmod__(other)
    return q

  def __mod__(self, other):
    _, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
    if isinstance(other, Poly):
      q, r = self._poly().__divmod__(other._poly())
    else:
      q, r = self._poly().__divmod__(other)
    return Poly(q), Poly(r)


  def __deepcopy__(self):
    return Poly(self._poly())

  def _poly(self) -> PolySparse | PolyDense:
    if not (isinstance(self._rep, PolyDense) or isinstance(self._rep, PolySparse)):
      # This code should never be reached.
      raise AttributeError(f"Poly object {self} contains neither a sparse or a dense representation.")
    else:
      return self._rep

  def is_sparse(self):
    return isinstance(self._rep, PolySparse)

  def poly_ring_mod(self):
    return self._poly().poly_ring_mod()

  def coeff_field_order(self):
    return self._poly().coeff_field_order()
