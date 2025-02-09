def poly_dense_add(
    a: list[int],
    b: list[int],
    coeff_field_order: int) -> list[int]:
  """
  Adds two polynomials.

  Args:
    poly_ring_mod (list[int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list[int]): The first polynomial.
    b (list[int]): The second polynomial.

  Returns:
    c (list[int]) The sum of the two polynomials.
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
    a: list[int],
    b: list[int],
    coeff_field_order: int) -> list[int]:
  """
  Subtracts two polynomials.

  Args:
    a (list[int]): The first polynomial.
    b (list[int]): The second polynomial.
    poly_ring_mod (list[int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    c (list[int]) The difference `a` - `b` of the two polynomials.
  """
  neg_b = poly_dense_scale(-1, b, coeff_field_order)
  return poly_dense_add(a, neg_b, coeff_field_order)


def poly_dense_scale(
    s: int,
    a: list[int],
    coeff_field_order: int) -> list[int]:
  """
  Scales a polynomial.

  Args:
    s (int): The scalar.
    a (list[int]): The polynomial to scale.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    a (list[int]) The scaled polynomial.
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
    a: list[int],
    b: list[int],
    coeff_field_order: int) -> list[int]:
  """
  Reduces a polynomial a(x) to the lowest degree polynomial b(x) such that
  a(x) = b(x) + g(x)poly_ring_mod(x) for some g(x).

  Args:
    a (list[int]): The polynomial to reduce.
    b (list[int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    q (list[int]) The quotient.
    r (list[int]) The remainder.
  """
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
      to_subtract = poly_dense_scale(r_leading_coeff,
                               monic_b,
                               coeff_field_order)
      to_subtract = poly_dense_multiply(degree_n_dense_monomial(len(r) - n),
                                  to_subtract,
                                  None,
                                  coeff_field_order)
      r = poly_dense_subtract(r, to_subtract, coeff_field_order)

  # Quotient terms currently in order of degree, decreasing.
  q.reverse()

  return q, r


def poly_dense_multiply(
    a: list[int],
    b: list[int],
    poly_ring_mod: list[int],
    coeff_field_order: int) -> list[int]:
  """
  Multiplies two polynomials.

  Args:
    poly_ring_mod (list[int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list[int]): The first polynomial.
    b (list[int]): The second polynomial.

  Returns:
    c (list[int]) The product of the two polynomials.
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
    a: list[int],
    n: int,
    poly_ring_mod: list[int],
    coeff_field_order: int,
    mod_is_irreducible: bool=False) -> list[int]:
  """
  Raises a polynomial to a power.

  Args:
    poly_ring_mod (list[int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.
    a (list[int]): The polynomial to raise to a power.x
    n (int): The power to raise the polynomial to.

  Returns:
    a (list[int]) The polynomial raised to the power.
  """
  res = [1]
  factor = a.copy()

  # Error handling.
  if n < 0:
    if poly_ring_mod == None:
      raise ValueError("Cannot invert an element whose domain isn't a field (poly_ring_mod == None)")
    else:
      poly_field_order = pow(coeff_field_order, len(poly_ring_mod) + 1)
      n = n % poly_field_order

  if mod_is_irreducible:
    if poly_ring_mod == None:
      raise ValueError("Cannot specify an irreducible modulus when no modulus provided")
    n %= coeff_field_order ** (len(poly_ring_mod) - 1)

  # Algorithm.
  while n > 0:
    if n % 2 == 1:
      res = poly_dense_multiply(res, factor, poly_ring_mod, coeff_field_order)
      n = n - 1 # For mathematical clarity.

    factor = poly_dense_multiply(factor, factor, poly_ring_mod, coeff_field_order)
    n = n / 2

  return res


class PolyDense:
  def __init__(
      self,
      coeffs: list[int],
      poly_ring_mod=None,
      coeff_field_order=None,
      *args,
      **kwargs):
    # check coeffs
    if coeffs == None:
      self._coeffs = []
    elif coeffs.__class__ == dict:
      self._coeffs = poly_sparse_to_dense(coeffs)
    elif coeffs.__class__ == list:
      self._coeffs = coeffs
    else:
      raise ValueError(f"coeffs must be a list or dictionary. recieved an instance of {coeffs.__class__}")
    # check poly_ring_mod
    if poly_ring_mod != None:
      if poly_ring_mod.__class__ == dict:
        self._poly_ring_mod = poly_sparse_to_dense(poly_ring_mod)
      elif poly_ring_mod.__class__ == list:
        self._poly_ring_mod = poly_ring_mod
      else:
        raise ValueError(f"poly_ring_mod must be a list or dictionary. recieved an instance of {poly_ring_mod.__class__}")
    else:
      self._poly_ring_mod = None
    # check coeff_field_order
    if coeff_field_order == None:
      self._coeff_field_order = None
    elif coeff_field_order.__class__ == int:
      self._coeff_field_order = coeff_field_order
    else:
      raise ValueError(f"coeff_field_order must be a list or dictionary. recieved an instance of {coeff_field_order.__class__}")

  def __str__(self):
    return poly_string(self._coeffs)

  def __add__(self, other):
    return PolyDense(
        poly_dense_add(
            self._coeffs,
            other._coeffs,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __sub__(self, other):
    return PolyDense(
        poly_dense_subtract(
            self._coeffs,
            other._coeffs,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __mul__(self, other):
    return PolyDense(
        poly_dense_multiply(
            self._coeffs,
            other._coeffs,
            poly_ring_mod=self._poly_ring_mod,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

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
    q, r = self.__divmod__(other)
    return q

  def __mod__(self, other):
    q, r = self.__divmod__(other)
    return r

  def __divmod__(self, other):
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

  def __copy__(self):
    return PolyDense(
        self._coeffs.copy(),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  # Alias for container class
  def _poly(self):
    return self
