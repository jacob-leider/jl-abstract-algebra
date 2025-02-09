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


def poly_string(
    poly: (list[int]|dict[int, int]),
    space_after_coeff: bool=False,
    order_increasing=True,
    implicit_powers=True,
    implicit_coeffs=True, **kwargs) -> str:
  """
  Returns a string representation of a polynomial.

  Args:
    poly (list[int]): The polynomial to print.

  Returns:
    poly_str (str): The string representation of the polynomial.
  """
  # Check param mode.
  mode = "txt"
  if mode in kwargs:
    mode = kwargs["mode"]
  # Convert polynomial to a list of (deg, coeff) pairs.
  deg_coeff_list = None
  if poly.__class__ == dict:
    deg_coeff_list = list(sorted(poly.items()))
  elif poly.__class__ == list:
    deg_coeff_list = list(enumerate(poly))
  else:
    raise ValueError(f"Invalid polynomial type: {poly.__class__}")
  if not order_increasing:
    deg_coeff_list.reverse()
  # Build a list of monomial strings.
  mono_str_vec = []
  for deg, coeff in deg_coeff_list:
    if coeff != 0:
      # Coefficient.
      if coeff == 1:
        if implicit_coeffs & deg > 0: # if deg == 0, coeff must be explicit.
          mono_str = ""
        else:
          mono_str = "1"
      elif coeff == -1:
        if implicit_coeffs:
          mono_str = "-"
        else:
          mono_str = "-1"
      else:
        mono_str = f"{coeff}"
      # Space.
      if space_after_coeff:
        mono_str += " "
      # Variable (x^).
      if deg == 0:
        if implicit_powers:
          mono_str_vec.append(mono_str)
          continue
        else:
          mono_str += "x^"
      elif deg == 1:
        if implicit_powers:
          mono_str += "x"
          mono_str_vec.append(mono_str)
          continue
        else:
          mono_str += "x^"
      else:
        mono_str += "x^"
      # Exponent.
      if mode == "latex":
        mono_str += "{" + str(deg) + "}"
      else:
        mono_str += str(deg)
      # Done.
      mono_str_vec.append(mono_str)

  # Join and clean.
  poly_str = " + ".join(mono_str_vec)
  if implicit_powers:
    poly_str = poly_str.replace("x^0", "")
    poly_str = poly_str.replace("x^1", "x")
  if implicit_coeffs:
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
        coeffs_dense = poly_parse_string(coeffs)
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
