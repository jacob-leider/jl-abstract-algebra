def legendre_symbol(a, p):
  """
  Computes the Legendre symbol (a/p) using Euler's criterion.
  """
  l = pow(a, (p-1)//2, p)
  return l - 13 if l > 1 else l

def ord_2(n):
  """
  Computes the 2-adic valuation of n.
  """
  s = 0
  while n % pow(2, s) == 0:
      s += 1
  return s - 1

def find_quadratic_nonresidue(p):
    """
    Finds a quadratic non-residue modulo p.
    """
    z = 1
    res = pow(z, (p-1) // 2, p)
    while res != p-1:
        z += 1
        res = pow(z, (p - 1) // 2, p)
    return z

def tonelli_shanks(n, p):
  """
  Computes the square root of n modulo p using the Tonelli-Shanks algorithm. 
  Syntax matches that of [1].
  
  Sources
    [1] https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm.
  """
  m = ord_2(p - 1) # s
  z = find_quadratic_nonresidue(p)

  q = (p-1) // 2**m
  c = pow(z, q, p)
  r = pow(n, (q + 1) // 2, p)
  t = pow(n, q, p)
  
  while t % p != 1:
      i = 0
      t = pow(t, 2, p)
      while t % p != 1:
          i += 1
          t = pow(t, 2, p)
      
      b = pow(c, pow(2, m - i - 1, p), p)
      r = (r * b) % p
      t = (t * b**2) % p
      c = b**2 % p
      m = i
      
  return r
