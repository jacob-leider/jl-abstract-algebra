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


def generalized_fast_pow(a, n: int, func: callable, *args, **kwargs):
  res = a.__deepcopy__() # Start out with first power since we don't know a^0
  n -= 1

  while n > 0:
    if n % 2 == 1:
      res = func(res, a, *args, **kwargs)
      n = n - 1 # For mathematical clarity.

    a = func(a, a, *args, **kwargs)
    n = n / 2

  return res

class Poly:
  def __init__(
      self,
      coeffs: list[int] | dict[int, int] | str | PolyDense | PolySparse,
      poly_ring_mod = None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    """
    A representation of an element of (possible the quotient of) a polynomial ring where 
    implementation details (e.g. sparse vs. dense) are abstracted away.s
    """
    self.is_sparse_ = False
    self.poly_sparse_ = None
    self.poly_dense_ = None
    # Check coeffs    
    if coeffs.__class__ == PolyDense:
      self.is_sparse_ = False
      self.poly_dense_ = coeffs
      self.poly_sparse_ = None
    elif coeffs.__class__ == PolySparse:
      self.is_sparse_ = True
      self.poly_sparse_ = coeffs
      self.poly_dense_ = None
    elif coeffs.__class__ == dict:
      self.is_sparse_ = True
      self.poly_sparse_ = PolySparse(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == list:
      self.is_sparse_ = False
      self.poly_dense_ = PolyDense(
          coeffs,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    elif coeffs.__class__ == str:
      # For now return a dense representation.
      try:
        coeffs_dense = poly_parse_string(coeffs)
      except SyntaxError:
        raise SyntaxError("Invalid polynomial string")
      self.is_sparse_ = False
      self.poly_dense_ = PolyDense(
          coeffs_dense,
          poly_ring_mod=poly_ring_mod,
          coeff_field_order=coeff_field_order)
    else:
      raise ValueError("coeffs must be a list, dictionary, or string")
    

  def __str__(self):
    
    if self.is_sparse():
      return self.poly_sparse_.__str__()
    else:
      return self.poly_dense_.__str__()

  def __add__(self, other):
    if self.is_sparse_:
      return Poly(self.poly_sparse_ + other.poly_sparse_)
    else:
      return Poly(self.poly_dense_ + other.poly_dense_)

  def __sub__(self, other):
    if self.is_sparse_:
      return Poly(self.poly_sparse_ - other.poly_sparse_)
    else:
      return Poly(self.poly_dense_ - other.poly_dense_)

  def __mul__(self, other):
    if self.is_sparse_:
      return Poly(self.poly_sparse_ * other.poly_sparse_)
    else:
      return Poly(self.poly_dense_ * other.poly_dense_)

  def __pow__(self, n):
    return generalized_fast_pow(self, n, Poly.__mul__)

  def __floordiv__(self, other):
    q, r = self.__divmod__(other)
    return Poly(q)

  def __mod__(self, other):
    q, r = self.__divmod__(other)
    return Poly(r)

  def __divmod__(self, other):
    if self.is_sparse_:
      q, r = self.poly_sparse_.__divmod__(other.poly_sparse_)
    else:
      q, r = self.poly_dense_.__divmod__(other.poly_dense_)
    return Poly(q), Poly(r)

  def __deepcopy__(self):
    if self.is_sparse_:
      return Poly(self.poly_sparse_.coeffs_,
                  self.poly_sparse_.poly_ring_mod_,
                  self.poly_sparse_.coeff_field_order_)
    else:
      return Poly(self.poly_dense_.coeffs_,
                  self.poly_dense_.poly_ring_mod_,
                  self.poly_dense_.coeff_field_order_)

  def is_sparse(self):
    return self.is_sparse_
