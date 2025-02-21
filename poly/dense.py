from poly.utils import *

def poly_dense_add(
    a: list,
    b: list,
    coeff_field_order: int | None) -> list:
  """
  Adds two polynomials.

  Args:
    poly_ring_mod (list): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list): The first polynomial.
    b (list): The second polynomial.

  Returns:
    c (list) The sum of the two polynomials.
  """
  c = [0] * max(len(a), len(b))

  for i in range(len(a)):
    c[i] += a[i]
    if coeff_field_order != None:
      c[i] %= coeff_field_order

  for i in range(len(b)):
    c[i] += b[i]
    if coeff_field_order != None:
      c[i] %= coeff_field_order

  c = clear_leading_zeros(c)

  return c


def poly_dense_subtract(
    a: list,
    b: list,
    coeff_field_order: int | None) -> list:
  """
  Subtracts two polynomials.

  Args:
    a (list): The first polynomial.
    b (list): The second polynomial.
    poly_ring_mod (list): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    c (list) The difference `a` - `b` of the two polynomials.
  """
  neg_b = poly_dense_scale(-1, b, coeff_field_order)
  return poly_dense_add(a, neg_b, coeff_field_order)


def poly_dense_scale(
    s: int,
    a: list,
    coeff_field_order: int | None) -> list:
  """
  Scales a polynomial.

  Args:
    s (int): The scalar.
    a (list): The polynomial to scale.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    a (list) The scaled polynomial.
  """
  scaled_a = a.copy()

  for i in range(len(scaled_a)):
    scaled_a[i] *= s
    if coeff_field_order != None:
      scaled_a[i] %= coeff_field_order

  return scaled_a


def degree_n_dense_monomial(n: int):
  """
  Returns the coefficient list of a monic monomial of degree n.

  Args:
    n (int): The degree of the monomial.

  Returns:
    [0] * n + [1] The monomial of degree n.
  """
  return [0] * n + [1]


def poly_dense_quotient_remainder(
    a: list,
    b: list,
    coeff_field_order: int | None) -> tuple[list, list]:
  """
  Reduces a polynomial a(x) to the lowest degree polynomial b(x) such that
  a(x) = b(x) + g(x)poly_ring_mod(x) for some g(x).

  Args:
    a (list): The polynomial to reduce.
    b (list): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    q (list) The quotient.
    r (list) The remainder.
  """

  if not isinstance(coeff_field_order, int):
    raise ValueError("Currently only supports prime fields (coeff_field_order) must be an int")

  r = a.copy()
  q = []
  n = len(b)
  N = len(r)

  # Division-by-zero-error.
  if (len(b) == 0):
    raise ValueError("dividend (`b`) is zero")

  # Check for ill-formatted polynomials.
  if (a[-1] == 0):
    raise ValueError("`a` has a leading zero coefficient")
  if (b[-1] == 0):
    raise ValueError("`b` has a leading zero coefficient")

  b_leading_coeff_inv = pow(b[-1], -1, coeff_field_order)
  monic_b = poly_dense_scale(b_leading_coeff_inv, b, coeff_field_order)

  for i in range(N - 1, n - 2, -1):
    if i >= len(r):
      q.append(0)
    else:
      # Add a term to the quotient.
      r_leading_coeff = r[i]
      q.append((r_leading_coeff * b_leading_coeff_inv) % coeff_field_order)
      # Subtract off the leading term from the remainder.
      to_subtract = poly_dense_scale(
          r_leading_coeff, 
          monic_b, 
          coeff_field_order)
      to_subtract = poly_dense_multiply(degree_n_dense_monomial(
          len(r) - n), 
          to_subtract, 
          None,
          coeff_field_order)
      r = poly_dense_subtract(r, to_subtract, coeff_field_order)

  # Quotient terms currently in order of degree, decreasing.
  q.reverse()

  return q, r


def poly_dense_multiply(
    a: list,
    b: list,
    poly_ring_mod: list | None,
    coeff_field_order: int | None) -> list:
  """
  Multiplies two polynomials.

  Args:
    poly_ring_mod (list): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list): The first polynomial.
    b (list): The second polynomial.

  Returns:
    c (list) The product of the two polynomials.
  """
  c = [0] * (len(a) + len(b))

  # Convolve.
  for i, x in enumerate(a):
    for j, y in enumerate(b):
      c[i + j] += x * y
      if coeff_field_order != None:
        c[i + j] %= coeff_field_order

  c = clear_leading_zeros(c)

  # Reduce.
  if poly_ring_mod != None:
    _, c = poly_dense_quotient_remainder(c, poly_ring_mod, coeff_field_order)

  return c


def poly_dense_fast_pow(
    a: list,
    n: int,
    poly_ring_mod: list | None,
    coeff_field_order: int | None,
    mod_is_irreducible: bool=False) -> list:
  """
  Raises a polynomial to a power.

  Args:
    poly_ring_mod (list): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list): The polynomial to raise to a power.x
    n (int): The power to raise the polynomial to.

  Returns:
    a (list) The polynomial raised to the power.
  """
  res = [1]
  factor = a.copy()

  # Error handling.
  if isinstance(poly_ring_mod, list):
    if n < 0:
      if coeff_field_order != None:
        if mod_is_irreducible:
          poly_field_order = pow(coeff_field_order, len(poly_ring_mod) + 1)
          n = n % poly_field_order
    else:
      if coeff_field_order != None:
        poly_field_order = pow(coeff_field_order, len(poly_ring_mod) + 1)
        n = n % poly_field_order
  elif poly_ring_mod == None:
    if n < 0:
      raise ValueError("Cannot invert an element whose domain isn't a field (poly_ring_mod == None)")
    if mod_is_irreducible:
      raise ValueError("Cannot specify an irreducible modulus when no modulus provided")
  else:
      raise ValueError("poly_ring_mod must be a list")

  # Algorithm.
  while n > 0:
    if n % 2 == 1:
      res = poly_dense_multiply(res, factor, poly_ring_mod, coeff_field_order)
      n = n - 1 # For mathematical clarity.

    factor = poly_dense_multiply(factor, factor, poly_ring_mod, coeff_field_order)
    n = n // 2

  return res


class PolyDense:
  def __init__(
      self,
      coeffs: list,
      poly_ring_mod: list | None=None,
      coeff_field_order: int | None=None,
      *args,
      **kwargs):
    # check coeffs
    self._coeffs = coeffs
    self._poly_ring_mod = poly_ring_mod
    self._coeff_field_order = coeff_field_order

  def __str__(self):
    return poly_string(self._coeffs)

  def __add__(self, other):
    if isinstance(other, PolyDense):
      return PolyDense(
         poly_dense_add(
              self._coeffs,
              other._coeffs,
              coeff_field_order=self._coeff_field_order),
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
    else:
      return NotImplemented

  def __sub__(self, other):
    if isinstance(other, PolyDense):
      return PolyDense(
          poly_dense_subtract(
              self._coeffs,
              other._coeffs,
              coeff_field_order=self._coeff_field_order),
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
    else:
      return NotImplemented

  def __rmul__(self, other):
    if isinstance(other, int):
      return PolyDense(
          poly_dense_scale(
              other,
              self._coeffs,
              coeff_field_order=self._coeff_field_order))
    elif isinstance(other, PolyDense):
      return other.__mul__(self)
    else:
      return NotImplemented

  def __mul__(self, other):
    if isinstance(other, int):
      return PolyDense(
          poly_dense_scale(
              other, 
              self._coeffs, 
              coeff_field_order=self._coeff_field_order),
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
    elif isinstance(other, PolyDense):
      return PolyDense(
          poly_dense_multiply(
              self._coeffs,
              other._coeffs,
              poly_ring_mod=self._poly_ring_mod,
              coeff_field_order=self._coeff_field_order),
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
    else:
      return NotImplemented

  def __pow__(self, n):
    return PolyDense(
        poly_dense_fast_pow(
            self._coeffs,
            n,
            poly_ring_mod=self._poly_ring_mod,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __floordiv__(self, other):
    q, _ = self.__divmod__(other)
    return q

  def __mod__(self, other):
    _, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
    # TODO: Scalar division makes sense. Just invert, scalar multiply, and 
    # return a remainder of 0.
    if isinstance(other, PolyDense):
      q, r = poly_dense_quotient_remainder(
          self._coeffs,
          other._coeffs,
          coeff_field_order=self._coeff_field_order)
      q = PolyDense(
          q,
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
      r = PolyDense(
          r,
          poly_ring_mod=self._poly_ring_mod,
          coeff_field_order=self._coeff_field_order)
      return q, r
    else:
      return NotImplemented

  def __copy__(self):
    return PolyDense(
        self._coeffs.copy(),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  # Alias for container class
  def _poly(self) -> 'PolyDense':
    return self

  def coeffs(self) -> list:
      return self._coeffs

  def coeff_field_order(self) -> int | None:
    return self._coeff_field_order

  def poly_ring_mod(self) -> list | None:
    return self._poly_ring_mod

