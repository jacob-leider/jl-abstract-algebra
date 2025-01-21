

def clear_leading_zeros(a: list[int]) -> list[int]:
  """
  Clears leading zeros from a list.

  Args:
    a (list[int]): The list to clear.

  Returns:
    a (list[int]) The list without leading zeros.
  """
  while len(a) > 0:
    if a[-1] == 0:
      a = a[:-1]
    else:
      break
  return a


def poly_add(a: list[int], b: list[int], coeff_field_order: int) -> list[int]:
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


def poly_subtract(a: list[int], b: list[int], coeff_field_order) -> list[int]:
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
  return poly_add(a, poly_scale(-1, b, coeff_field_order), coeff_field_order)


def poly_scale(s: int, a: list[int], coeff_field_order) -> list[int]:
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


def degree_n_monomial(n):
  """
  Returns the coefficient list of a monic monomial of degree n.

  Args:
    n (int): The degree of the monomial.

  Returns:
    [0] * n + [1] The monomial of degree n.
  """
  return [0] * n + [1]


def poly_quotient_remainder(a: list[int], b: list[int], coeff_field_order) -> list[int]:
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
  monic_b = poly_scale(b_leading_coeff_inv, b, coeff_field_order)

  for i in range(N - 1, n - 2, -1):
    if i >= len(r):
      q.append(0)
    else:
      # Add a term to the quotient.
      r_leading_coeff = r[i]
      q.append((r_leading_coeff * b_leading_coeff_inv) % coeff_field_order)
      # Subtract off the leading term from the remainder.
      to_subtract = poly_scale(r_leading_coeff,
                               monic_b,
                               coeff_field_order)
      to_subtract = poly_multiply(degree_n_monomial(len(r) - n),
                                  to_subtract,
                                  None,
                                  coeff_field_order)
      r = poly_subtract(r, to_subtract, coeff_field_order)

  # Quotient terms currently in order of degree, decreasing.
  q.reverse()

  return q, r


def poly_multiply(a: list[int], b: list[int], poly_ring_mod, coeff_field_order) -> list[int]:
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
    _, c = poly_quotient_remainder(c, poly_ring_mod, coeff_field_order)

  return c


def poly_fast_pow(a: list[int], n: int, poly_ring_mod, coeff_field_order, mod_is_irreducible=False) -> list[int]:
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
      res = poly_multiply(res, factor, poly_ring_mod, coeff_field_order)
      n = n - 1 # For mathematical clarity.

    factor = poly_multiply(factor, factor, poly_ring_mod, coeff_field_order)
    n = n / 2

  return res
