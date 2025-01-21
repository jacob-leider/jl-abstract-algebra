

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
