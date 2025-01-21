

def poly_sparse_scale(s: int, a: dict[int, int], coeff_field_order) -> dict[int, int]:
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


def poly_sparse_add(a: dict[int, int], b: dict[int, int], coeff_field_order: int) -> dict[int, int]:
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


def poly_sparse_subtract(a: dict[int, int], b: dict[int, int], coeff_field_order: int) -> dict[int, int]:
  """
  Subtracts two sparse polynomials.

  Args:
    a (dict[int, int]): The first polynomial.
    b (dict[int, int]): The second polynomial.
    coeff_field_order (int): mod for coefficient field.

  Returns:
    c (dict[int, int]) The sum a - b.
  """
  return poly_sparse_add(a, poly_sparse_scale(-1, b, coeff_field_order), coeff_field_order)


def poly_sparse_multiply(a: dict[int, int], b: dict[int, int], poly_ring_mod, coeff_field_order) -> dict[int, int]:
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
    pass

  return c


def poly_sparse_quotient_remainder(a: dict[int, int], b: dict[int, int], coeff_field_order: int) -> dict[int, int]:
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
  r_nz = sorted(a.keys(), reverse = True)
  b_nz = sorted(b.keys(), reverse = True)

  q_map = {}
  r_map = a.copy()
  
  deg_g = b_nz[0]
  leading_g_coeff_inv = pow(b[deg_g], -1, coeff_field_order)
  monic_b_map = poly_sparse_scale(leading_g_coeff_inv, b.copy(), coeff_field_order)

  ct = 10
  for deg_r in r_nz:
    if deg_r < deg_g:
      break
    
    ct -= 1
    if ct == 0:
      break
    
    if deg_r not in r_map:
      # Term vanished in a prior round.
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
  
  # Ensure proper form.
  for deg, coeff in r_map.copy().items():
    if coeff == 0:
      del r_map[deg]
  for deg, coeff in q_map.copy().items():
    if coeff == 0:
      del q_map[deg]

  return q_map, r_map


def degree_n_sparse_monomial(n):
  """
  Returns the coefficient list of a monic monomial of degree n.

  Args:
    n (int): The degree of the monomial.

  Returns:
    [0] * n + [1] The monomial of degree n.
  """
  m = {}
  m[n] = 1
  return m

def degree_n_monomial(n):
  """
  Returns the coefficient list of a monic monomial of degree n.

  Args:
    n (int): The degree of the monomial.

  Returns:
    [0] * n + [1] The monomial of degree n.
  """
  return [0] * n + [1]
