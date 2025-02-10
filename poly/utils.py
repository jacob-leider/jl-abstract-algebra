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

  poly_str = " + ".join(mono_str_vec)
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

