"""Polynomial Interface

Provides access to dense and sparse polynomial representations through a single 
`Poly` class.

"""

from .sparse import *
from .dense import *
from .utils import *

# May need to make coeff_field_order a kwarg, so we know if user intended to 
# switch the field when passing a Poly to Poly.

# Conventions
#   1. If the representation is ambiguous, always revert to dense.
#   2. if coeffs specifies a representation (i.e. is a list, dict, PolyDense 
#      or PolySparse), use this representation for the Poly instance.


def make_sparse(poly, **kwargs) -> dict | None:
    if isinstance(poly, Poly) or isinstance(poly, PolyDense) or isinstance(poly, PolySparse):
        poly = poly.coeffs()

    if isinstance(poly, list):
        return poly_dense_to_sparse(poly)
    elif isinstance(poly, dict):
        return poly
    elif isinstance(poly, str):
        var = "x"
        if "var" in kwargs:
            var = kwargs["var"]
        try:
          return poly_dense_to_sparse(poly_parse_string(poly, var=var))
        except SyntaxError:
          raise SyntaxError("Invalid polynomial string")
    elif isinstance(poly, int):
        return {0: poly}
    else:
        return None

def make_dense(poly, **kwargs) -> list | None:
    if isinstance(poly, Poly) or isinstance(poly, PolyDense) or isinstance(poly, PolySparse):
        poly = poly.coeffs()

    if isinstance(poly, list):
        return poly
    elif isinstance(poly, dict):
        return poly_sparse_to_dense(poly)
    elif isinstance(poly, str):
        var = "x"
        if "var" in kwargs:
            var = kwargs["var"]
        try:
          return poly_parse_string(poly, var=var)
        except SyntaxError:
          raise SyntaxError("Invalid polynomial string")
    elif isinstance(poly, int):
        return [poly]
    else:
        return None


def determine_rep(coeffs, **kwargs):
    rep = "dense"

    # Second highest priority.
    if isinstance(coeffs, PolyDense) or isinstance(coeffs, list):
        rep = "dense"
    elif isinstance(coeffs, PolySparse) or isinstance(coeffs, dict):
        rep = "sparse"

    # Highest priority.
    if "rep" in kwargs:
        if kwargs["rep"] == "dense":
            pass
        elif kwargs["rep"] == "sparse":
            pass
        else:
            raise ValueError("keyword arg rep must be \"dense\" or \"sparse\"")

    return rep


class Poly:
  def __init__(
      self,
      coeffs: "list | dict | str | Poly | PolyDense | PolySparse | int",
      poly_ring_mod: "list | dict | str | 'Poly' | PolyDense | PolySparse | int | None"=None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    rep = determine_rep(coeffs, **kwargs)

    if rep == "dense":
        dense_coeffs = make_dense(coeffs, **kwargs)
        if dense_coeffs is None:
            raise ValueError(f"bad coeffs (type = {type(coeffs)})")
        else:
            self._rep = PolyDense(
                    dense_coeffs,
                    poly_ring_mod=make_dense(poly_ring_mod, **kwargs),
                    coeff_field_order=coeff_field_order)
    else: # rep == "sparse"
        sparse_coeffs = make_sparse(coeffs, **kwargs)
        if sparse_coeffs is None:
            raise ValueError(f"bad coeffs (type = {type(coeffs)})")
        else:
            self._rep = PolySparse(
                    sparse_coeffs,
                    poly_ring_mod=make_sparse(poly_ring_mod, **kwargs),
                    coeff_field_order=coeff_field_order)

  def __str__(self, **kwargs):
    return poly_string(self._poly()._coeffs, **kwargs)

  def __add__(self, other):
    if isinstance(other, Poly):
      return self._copy_ring(self._poly() + other._poly())
    else:
      return self._copy_ring(self._poly() + other)

  def __sub__(self, other):
    if isinstance(other, Poly):
      return self._copy_ring(self._poly() - other._poly())
    else:
      return self._copy_ring(self._poly() - other)

  def __rmul__(self, other):
    if isinstance(other, Poly):
      return self._copy_ring(self._poly() * other._poly())
    else:
      return self._copy_ring(self._poly() * other)

  def __mul__(self, other):
    if isinstance(other, Poly):
      return self._copy_ring(self._poly() * other._poly())
    else:
      return self._copy_ring(self._poly() * other)

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
    return self._copy_ring(q), self._copy_ring(r)


  def __deepcopy__(self):
    return self._copy_ring(self._poly())

  def _poly(self) -> PolySparse | PolyDense:
    if not (isinstance(self._rep, PolyDense) or isinstance(self._rep, PolySparse)):
      # This code should never be reached.
      raise AttributeError(f"Poly object {self} contains neither a sparse or a dense representation.")
    else:
      return self._rep

  def _copy_ring(self, coeffs):
    """
    Create a new instance copying all attributes from this instance except 
    _coeffs.
    """
    return Poly(coeffs, poly_ring_mod=self.poly_ring_mod(), 
                coeff_field_order=self.coeff_field_order())

  def is_sparse(self):
    return isinstance(self._rep, PolySparse)

  def poly_ring_mod(self):
    return self._poly().poly_ring_mod()

  def coeff_field_order(self):
    return self._poly().coeff_field_order()

  def coeffs(self) -> list | dict:
    return self._poly().coeffs()
