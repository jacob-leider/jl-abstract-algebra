from sparse import *
from dense import *


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


def poly_string(poly: (list[int]|dict[int, int]), space_after_coeff=False) -> str:
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
      coeffs: list[int],
      poly_ring_mod=None,
      coeff_field_order=None,
      *args,
      **kwargs):
    self.coeffs_ = clear_leading_zeros(coeffs)
    self.coeff_field_order_ = coeff_field_order
    self.poly_ring_mod_ = poly_ring_mod

  def __str__(self):
    return poly_string(self.coeffs_)

  def __add__(self, other):
    return Poly(
        poly_add(
            self.coeffs_,
            other.coeffs_,
            coeff_field_order=self.coeff_field_order_),
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)

  def __sub__(self, other):
    return Poly(
        poly_subtract(
            self.coeffs_,
            other.coeffs_,
            coeff_field_order=self.coeff_field_order_),
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)

  def __mul__(self, other):
    prm_coeffs = None
    if self.poly_ring_mod_ != None:
      prm_coeffs = self.poly_ring_mod_.coeffs_

    return Poly(
        poly_multiply(
            self.coeffs_,
            other.coeffs_,
            poly_ring_mod=prm_coeffs,
            coeff_field_order=self.coeff_field_order_),
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)

  def __pow__(self, n):
    prm_coeffs = None
    if self.poly_ring_mod_ != None:
      prm_coeffs = self.poly_ring_mod_.coeffs_

    return Poly(
        poly_fast_pow(
            self.coeffs_,
            n,
            poly_ring_mod=prm_coeffs,
            coeff_field_order=self.coeff_field_order_),
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)

  def __floordiv__(self, other):
    q, r = self.__divmod__(other)
    return q

  def __mod__(self, other):
    q, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
    q, r = poly_quotient_remainder(
        self.coeffs_,
        other.coeffs_,
        coeff_field_order=self.coeff_field_order_)
    q = Poly(
        q,
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)
    r = Poly(
        r,
        poly_ring_mod=self.poly_ring_mod_,
        coeff_field_order=self.coeff_field_order_)
    return q, r

  def sparsity(self):
    """
    Returns the sparsity of the polynomial's list representation.
    """
    return self.coeffs_.count(0) / len(self.coeffs_)
