#! /usr/bin/env python
#
# Provide some simple capabilities from number theory.
#
# Version of 2008.11.14.
#
# Written in 2005 and 2006 by Peter Pearson and placed in the public domain.
# Revision history:
#   2008.11.14: Use pow( base, exponent, modulus ) for modular_exp.
#               Make gcd and lcm accept arbitrarly many arguments.


class Error( Exception ):
  """Base class for exceptions in this module."""
  pass


class SquareRootError( Error ):
  pass


def polynomial_reduce_mod( poly, polymod, p ):
  """Reduce poly by polymod, integer arithmetic modulo p.

  Polynomials are represented as lists of coefficients
  of increasing powers of x."""

  # This module has been tested only by extensive use
  # in calculating modular square roots.

  # Just to make this easy, require a monic polynomial:
  assert polymod[-1] == 1

  assert len( polymod ) > 1

  while len( poly ) >= len( polymod ):
    if poly[-1] != 0:
      for i in range( 2, len( polymod ) + 1 ):
        poly[-i] = ( poly[-i] - poly[-1] * polymod[-i] ) % p
    poly = poly[0:-1]

  return poly



def polynomial_multiply_mod( m1, m2, polymod, p ):
  """Polynomial multiplication modulo a polynomial over ints mod p.

  Polynomials are represented as lists of coefficients
  of increasing powers of x."""

  # This is just a seat-of-the-pants implementation.

  # This module has been tested only by extensive use
  # in calculating modular square roots.

  # Initialize the product to zero:

  prod = ( len( m1 ) + len( m2 ) - 1 ) * [0]

  # Add together all the cross-terms:

  for i in range( len( m1 ) ):
    for j in range( len( m2 ) ):
      prod[i+j] = ( prod[i+j] + m1[i] * m2[j] ) % p

  return polynomial_reduce_mod( prod, polymod, p )


def polynomial_exp_mod( base, exponent, polymod, p ):
  """Polynomial exponentiation modulo a polynomial over ints mod p.

  Polynomials are represented as lists of coefficients
  of increasing powers of x."""

  # Based on the Handbook of Applied Cryptography, algorithm 2.227.

  # This module has been tested only by extensive use
  # in calculating modular square roots.

  assert exponent < p

  if exponent == 0: return [ 1 ]

  G = base
  k = exponent
  if k%2 == 1: s = G
  else:        s = [ 1 ]

  while k > 1:
    k = k // 2
    G = polynomial_multiply_mod( G, G, polymod, p )
    if k%2 == 1: s = polynomial_multiply_mod( G, s, polymod, p )

  return s



def jacobi( a, n ):
  """Jacobi symbol"""

  # Based on the Handbook of Applied Cryptography (HAC), algorithm 2.149.

  # This function has been tested by comparison with a small
  # table printed in HAC, and by extensive use in calculating
  # modular square roots.

  assert n >= 3
  assert n%2 == 1
  a = a % n
  if a == 0: return 0
  if a == 1: return 1
  a1, e = a, 0
  while a1%2 == 0:
    a1, e = a1//2, e+1
  if e%2 == 0 or n%8 == 1 or n%8 == 7: s = 1
  else: s = -1
  if a1 == 1: return s
  if n%4 == 3 and a1%4 == 3: s = -s
  return s * jacobi( n % a1, a1 )



def square_root_mod_prime( a, p ):
  """Modular square root of a, mod p, p prime."""

  # Based on the Handbook of Applied Cryptography, algorithms 3.34 to 3.39.

  # This module has been tested for all values in [0,p-1] for
  # every prime p from 3 to 1229.

  assert 0 <= a < p
  assert 1 < p

  if a == 0: return 0
  if p == 2: return a

  jac = jacobi( a, p )
  if jac == -1: raise SquareRootError( "%d has no square root modulo %d" \
                                       % ( a, p ) )

  if p % 4 == 3: return pow( a, (p+1)//4, p )

  if p % 8 == 5:
    d = pow( a, (p-1)//4, p )
    if d == 1: return pow( a, (p+3)//8, p )
    if d == p-1: return ( 2 * a * pow( 4*a, (p-5)//8, p ) ) % p
    raise RuntimeError("Shouldn't get here.")

  for b in range( 2, p ):
    if jacobi( b*b-4*a, p ) == -1:
      f = ( a, -b, 1 )
      ff = polynomial_exp_mod( ( 0, 1 ), (p+1)//2, f, p )
      assert ff[1] == 0
      return ff[0]
  raise RuntimeError("No b found.")

