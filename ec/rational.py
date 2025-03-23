import utils

def ec_get_rational_points(p, a, b) -> list:
  points: list[tuple[None, None] | tuple[int, int]] = [(None, None)]
  for x in range(p):
    ysquared = (pow(x, 3, p) + (a * x) % p + b) % p
    ls = utils.legendre_symbol(ysquared, p)
    if ls == 0:
      points.append((x, 0))
    elif ls == 1:
      y1 = utils.tonelli_shanks(ysquared, p)
      y2 = (y1 * (p - 1)) % p
      points.append((x, y1))
      points.append((x, y2))
  return points


def ec_count_rational_points(p, a, b) -> int:
  npoints: int = 1
  for x in range(p):
    ysquared = (pow(x, 3, p) + (a * x) % p + b) % p
    npoints += utils.legendre_symbol(ysquared, p) + 1
  return npoints

