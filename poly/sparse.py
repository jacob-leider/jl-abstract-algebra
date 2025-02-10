from poly.utils import *

def poly_sparse_scale(
    s: int,
    a: dict[int, int],
    coeff_field_order: None | int) -> dict[int, int]:
  """
  Scales a sparse polynomial.

  Args:
    s (int): The scalar.
    a (dict[int, int]): The sparse polynomial to scale.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    scaled_a (dict[int, int]) The scaled sparse polynomial.
  """
  scaled_a = {}
  for deg, coeff in a.items():
    if coeff_field_order == None:
      scaled_a[deg] = s * coeff
    else:
      scaled_a[deg] = (s * coeff) % coeff_field_order

  return scaled_a


def poly_sparse_add(
    a: dict[int, int],
    b: dict[int, int],
    coeff_field_order: int | None) -> dict[int, int]:
  """
  Adds two sparse polynomials.

  Args:
    a (dict[int, int]): The first sparse polynomial.
    b (dict[int, int]): The second sparse polynomial.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    c (dict[int, int]) The sum of the two sparse polynomials.
  """
  c = a.copy()

  for deg, coeff in b.items():
    if deg in c:
      c[deg] += coeff
    else:
      c[deg] = coeff
    if coeff_field_order != None:
      c[deg] %= coeff_field_order

    if c[deg] == 0:
      del c[deg]

  return c


def poly_sparse_subtract(
    a: dict[int, int],
    b: dict[int, int],
    coeff_field_order: int | None) -> dict[int, int]:
  """
  Subtracts two sparse polynomials.

  Args:
    a (dict[int, int]): The first polynomial.
    b (dict[int, int]): The second polynomial.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    c (dict[int, int]) The sum a - b.
  """
  neg_b = poly_sparse_scale(-1, b, coeff_field_order)
  return poly_sparse_add(a, neg_b, coeff_field_order)


def poly_sparse_multiply(
    a: dict[int, int],
    b: dict[int, int],
    poly_ring_mod: None | dict[int, int],
    coeff_field_order: None | int) -> dict[int, int]:
  c = {}

  for a_deg, a_coeff in a.items():
    for b_deg, b_coeff in b.items():
      if a_deg + b_deg in c:
        c[a_deg + b_deg] += a_coeff * b_coeff
      else:
        c[a_deg + b_deg] = a_coeff * b_coeff

      if coeff_field_order != None:
        c[a_deg + b_deg] %= coeff_field_order

      if c[a_deg + b_deg] == 0:
        del c[a_deg + b_deg]

  if poly_ring_mod != None:
    _, c = poly_sparse_quotient_remainder(c, poly_ring_mod, coeff_field_order)

  return c


def degree_n_sparse_monomial(n: int):
  """
  Returns the coefficient list of a monic monomial of degree n.

  Args:
    n (int): The degree of the monomial.

  Returns:
    The monomial of degree n.
  """
  return {n : 1}

def poly_sparse_quotient_remainder(
    a: dict[int, int],
    b: dict[int, int],
    coeff_field_order: None | int) -> tuple[dict[int, int], dict[int, int]]:
  """
  Reduces a (sparse) polynomial a(x) to the lowest degree (sparse) polynomial
  b(x) such that a(x) = b(x) + g(x)poly_ring_mod(x) for some g(x).

  Args:
    a (dict[int, int]): The polynomial to reduce.
    b (dict[int, int]): Mod for the polynomial ring.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    q (dict[int, int]) The remainder.
    r (dict[int, int]) The quotient.
  """
  
  # FIXME: Support any ring.
  if not isinstance(coeff_field_order, int):
    raise ValueError("Currently only accepts prime field coefficients.")

  # nz ~ nonzero.
  r_nz = sorted(a.keys(), reverse = True)
  b_nz = sorted(b.keys(), reverse = True)

  q_map = {}
  r_map = a.copy()

  deg_g = b_nz[0]
  leading_g_coeff_inv = pow(b[deg_g], -1, coeff_field_order)
  monic_b_map = poly_sparse_scale(
      leading_g_coeff_inv, # FIXME: This doesn't work with polynomial rings over general coefficient rings.
      b.copy(),
      coeff_field_order)

  for deg_r in r_nz:
    # Euclidean function-check.
    if deg_r < deg_g:
      break
    # Already zero-check.
    if deg_r not in r_map:
      continue

    # Subtract off a term.
    leading_r_coeff = r_map[deg_r]
    to_subtract = poly_sparse_multiply(
        degree_n_sparse_monomial(deg_r - deg_g),
        monic_b_map,
        None,
        coeff_field_order)
    to_subtract = poly_sparse_scale(
            leading_r_coeff,
        to_subtract,
        coeff_field_order)
    r_map = poly_sparse_subtract(r_map, to_subtract, coeff_field_order)

    # Append the appropriate q term.
    if deg_r - deg_g in q_map:
      q_map[deg_r - deg_g] += leading_r_coeff
    else:
      q_map[deg_r - deg_g] = leading_r_coeff
    if coeff_field_order != None:
      q_map[deg_r - deg_g] %= coeff_field_order

  # Ensure all zeros are implicit.
  for deg, coeff in r_map.copy().items():
    if coeff == 0:
      del r_map[deg]
  for deg, coeff in q_map.copy().items():
    if coeff == 0:
      del q_map[deg]

  return q_map, r_map


def poly_sparse_fast_pow(
    a: dict[int, int],
    n: int, poly_ring_mod: None | dict[int, int],
    coeff_field_order: int | None,
    mod_is_irreducible: bool=False) -> dict[int, int]:
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
  res = {0: 1}
  factor = a.copy()

  # Error handling.
  if poly_ring_mod == None:
      if n < 0:
        raise ValueError("Cannot invert an element whose domain isn't a field (poly_ring_mod == None)")
      if mod_is_irreducible:
        raise ValueError("Cannot specify an irreducible modulus when no modulus provided")
  else:
    if len(poly_ring_mod) == 0:
      raise ValueError("poly_ring_mod cannot be an empty dict.")
    # Reduce the power modulo the order of the multiplicative group of 
    # coeff_field_order[x] / (poly_ring_mod).
    poly_ring_mod_deg = max(poly_ring_mod.keys())
    if isinstance(coeff_field_order, int):
      poly_field_order = pow(coeff_field_order, poly_ring_mod_deg)
      if n < 0:
        n = n % poly_field_order
      if mod_is_irreducible:
        n %= poly_field_order

  # Algorithm.
  while n > 0:
    if n % 2 == 1:
      res = poly_sparse_multiply(
          res,
          factor,
          poly_ring_mod,
          coeff_field_order)
      n = n - 1 # For mathematical clarity.

    factor = poly_sparse_multiply(
        factor,
        factor,
        poly_ring_mod,
        coeff_field_order)
    n = n // 2

  return res


class PolySparse:
  def __init__(
      self,
      coeffs: dict[int, int] | int,
      poly_ring_mod: None | dict[int, int]=None,
      coeff_field_order=None,
      *args,
      **kwargs):
    """
    Sparse representation of an element of a polynomial ring.

    Args:
      coeffs (list[int] | dict[int, int]): The coefficients of the polynomial.
      poly_ring_mod (list[int] | dict[int,int]): Mod for the polynomial ring.
      coeff_field_order (int): mod for coefficient field.
    """
    # check coeffs
    if isinstance(coeffs, dict):
      self._coeffs = coeffs
    elif isinstance(coeffs, int):
      self._coeffs = {0:coeffs}
    else:
      raise ValueError(f"coeffs must be a dictionary. recieved an instance of {coeffs.__class__}")
    # check poly_ring_mod
    if poly_ring_mod != None:
      if isinstance(poly_ring_mod, dict):
        self._poly_ring_mod = poly_ring_mod
      else:
        raise ValueError(f"poly_ring_mod must be a dictionary. recieved an instance of {poly_ring_mod.__class__}")
    else:
      self._poly_ring_mod = None
    # check coeff_field_order
    if coeff_field_order == None:
      self._coeff_field_order = None
    elif isinstance(coeff_field_order, int):
      self._coeff_field_order = coeff_field_order
    else:
      raise ValueError(f"coeff_field_order must be an int. recieved an instance of {coeffs.__class__}")

  def __str__(self):
    return poly_string(self._coeffs)

  def __add__(self, other):
    return PolySparse(
        poly_sparse_add(
            self._coeffs,
            other._coeffs,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __sub__(self, other):
    return PolySparse(
        poly_sparse_subtract(
            self._coeffs,
            other._coeffs,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __mul__(self, other):
    return PolySparse(
        poly_sparse_multiply(
            self._coeffs,
            other._coeffs,
            poly_ring_mod=self._poly_ring_mod,
            coeff_field_order=self._coeff_field_order),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  def __pow__(self, n):
    return PolySparse(
        poly_sparse_fast_pow(
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
    q, r = poly_sparse_quotient_remainder(
        self._coeffs,
        other._coeffs,
        coeff_field_order=self._coeff_field_order)
    q = PolySparse(
        q,
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)
    r = PolySparse(
        r,
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)
    return q, r

  def __copy__(self):
    return PolySparse(
        self._coeffs.copy(),
        poly_ring_mod=self._poly_ring_mod,
        coeff_field_order=self._coeff_field_order)

  # Alias for container class
  def _poly(self) -> 'PolySparse':
    return self

  def coeffs(self) -> dict[int, int]:
      return self._coeffs
