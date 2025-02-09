# poly
Polynomial arithmetic, construction of polynomial rings and algebras, and various algorithms for algebraic structures.

## Documentation

`Poly` is the main class for polynomials. While instances of the `Poly` class maintain references to `PolySparse` and `PolyDense` instances, it is preferable to use the `Poly` class and its sparse/dense conversion utilities over `PolySparse` and `PolyDense`.

### Poly
A representation of a polynomial with integer, rational, or finite-field coefficients. Maintains a reference to a sparse representation (dict/map: see `PolySparse`) or a dense representation (list: See `PolyDense`). 

### PolySparse
A polynomial represented by a dict that maps degrees (ints) to coefficients.

### PolyDense
A polynomial represented by a list of coefficients where indices (starting from zero, increasing) represent degrees.
