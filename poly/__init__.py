"""Polynomial Interface

Provides access to dense and sparse polynomial representations through a single 
`Poly` class.

"""

from . import sparse
from . import dense

def poly_sparse_to_dense(sparse: dict[int, int]) -> list[int]:
  """
  Converts a sparse polynomial to a dense polynomial.

  Args:
    sparse (dict[int, int]): The sparse polynomial.

  Returns:
    dense (list[int]) The dense polynomial.
  """
  dense = []
  for i in range(max(sparse.keys()) + 1):
    if i in sparse:
      dense.append(sparse[i])
    else:
      dense.append(0)
  return dense


def poly_dense_to_sparse(dense: list[int]) -> dict[int, int]:
  """
  Converts a dense polynomial to a sparse polynomial.

  Args:
    dense (list[int]): The dense polynomial.

  Returns:
    sparse (dict[int, int]) The sparse polynomial.
    """
  sparse = {}
  for deg, coeff in enumerate(dense):
    if coeff != 0:
      sparse[deg] = coeff
  return sparse


def poly_string(poly: list[int] | dict[int, int], space_after_coeff=False) -> str:
  """
  Returns a string representation of a polynomial.

  Args:
    poly (list[int]): The polynomial to print.
  
  Returns:
    poly_str (str): The string representation of the polynomial.
  """
  if poly.__class__ == dict:
    mono_str_vec = []
    for deg, coeff in sorted(poly.items(), reverse=True):
      # Coefficients cannot be zero.
      if space_after_coeff:
        mono_str_vec.append(f"{coeff} x^{deg}")
      else:
        mono_str_vec.append(f"{coeff}x^{deg}")

    poly_str = " + ".join(mono_str_vec)
    poly_str = poly_str.replace("x^0", "")
    poly_str = poly_str.replace("1x", "x")
  else:
    if len(poly) == 0:
      return "0"

    mono_str_vec = []
    for i, c in enumerate(poly):
      if c != 0:
        if space_after_coeff:
          mono_str_vec.append(f"{c} x^{i}")
        else:
          mono_str_vec.append(f"{c}x^{i}")

    poly_str = " + ".join(reversed(mono_str_vec))
    poly_str = poly_str.replace("x^0", "")
    poly_str = poly_str.replace("1x", "x")
  
  return poly_str


def poly_parse_string(poly_str: str, var: str) -> list[int]:
  """
  Parses a polynomial string.

  Args:
    poly_str (str): The polynomial string to parse.

  Returns:
    poly (list[int]) The parsed polynomial.
  """
  poly_str = poly_str.replace(" ", "")
  poly_str = poly_str.replace("^", "")
  poly_str = poly_str.replace("-", "+-")

  coeffs_map = {}
  for mono_str in poly_str.split("+"):
    if var in mono_str:
      coeff, degree = mono_str.split(var)
      if coeff == '':
        # Degree 1 may be implicit.
        coeff = 1
      else:
        try:
          coeff = int(coeff)
        except:
          raise ValueError(f"Invalid coefficient: {coeff}")
      if degree == '':
        # Identity coefficient may be implicit.
        degree = 1
      else:
        try:
          degree = int(degree)
        except:
          raise ValueError(f"Invalid degree: {degree}")
      if degree in coeffs_map:
        coeffs_map[degree] += coeff
      else:
        coeffs_map[degree] = coeff
    else:
      # Degree 0 may be implicit.
      if len(mono_str) > 0:
        try:
          coeff = eval(mono_str)
          if 0 in coeffs_map:
            coeffs_map[0] += coeff
          else:
            coeffs_map[0] = coeff
        except:
          raise ValueError(f"Invalid coefficient: {mono_str}")
        
  return poly_sparse_to_dense(coeffs_map)


class Poly:
  def __init__(
      self,
      coeffs: list[int] | dict[int, int],
      poly_ring_mod: list[int] | dict[int, int]=None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    self.is_sparse_ = False
    self.poly_sparse_ = None
    self.poly_dense_ = None

    if coeffs.__class__ == dict[int, int]:
      self.is_sparse_ = True
      self.poly_sparse_ = sparse.PolySparse(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == list[int]:
      self.is_sparse_ = False
      self.poly_dense_ = dense.PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    else:
      raise ValueError("coeffs must be a list or dictionary")

  def __str__(self):
    return poly_string(self.coeffs_)

  def __add__(self, other):
    if self.is_sparse_:
      return self.poly_sparse_ + other
    else:
      return self.poly_dense_ + other

  def __sub__(self, other):
    if self.is_sparse_:
      return self.poly_sparse_ - other
    else:
      return self.poly_dense_ - other

  def __mul__(self, other):
    if self.is_sparse_:
      return self.poly_sparse_ * other
    else:
      return self.poly_dense_ * other

  def __pow__(self, n):
    if self.is_sparse_:
      return self.poly_sparse_ ** n
    else:
      return self.poly_dense_ ** n

  def __floordiv__(self, other):
    q, r = self.__divmod__(other)
    return q

  def __mod__(self, other):
    q, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
    if self.is_sparse_:
      q, r = self.poly_sparse_.__divmod__(other)
    else:
      q, r = self.poly_dense_.__divmod__(other)
    return q, r

  def is_sparse(self):
    return self.is_sparse_
