# -*- coding: utf-8 -*-

from .field import FQ, int_types


class EIP1829Error(Exception):
    pass


def doubling(x1, y1, a, b):
    # From: http://www.secg.org/sec1-v2.pdf §2.2.1 "Elliptic Curves over Fp"
    # https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
    x3 = ((3*(x1**2)+a)**2)/((2*y1)**2)-x1-x1
    y3 = (2*x1+x1)*(3*(x1**2)+a)/(2*y1)-(3*(x1**2)+a)**3/(2*y1)**3-y1
    return (x3, y3)


def addition(x1, y1, x2, y2, a, b):
    # From: http://www.secg.org/sec1-v2.pdf §2.2.1 "Elliptic Curves over Fp"
    # https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
    if x1 is None and y1 is None:
        return (x2, y2)
    if x1 == x2 and y1 == y2:
        return doubling(x1, y1, a, b)

    x3 = (y2-y1)**2/(x2-x1)**2-x1-x2
    y3 = (2*x1+x2)*(y2-y1)/(x2-x1)-(y2-y1)**3/(x2-x1)**3-y1

    return (x3, y3)


def multiply(X, Y, s, alpha, beta):
    # double and add
    rX, rY = None, None
    aX, aY = X, Y
    while s > 0:
        b = s & 1
        if b == 1:
            rX, rY = addition(rX, rY, aX, aY, alpha, beta)
        aX, aY = addition(aX, aY, aX, aY, alpha, beta)
        s = s // 2
    return (rX, rY)


def is_on_curve(x, y, a, b):
    """
    y² = x³ + α ⋅ x + β  mod  m
    """
    yy = y * y
    rhs = x**3 + (a * x) + b
    return yy == rhs


def recover_y(x, a, b):
    yy = x**3 + (a * x) + b
    return yy.sqrt()


def coerce_integer(value):
    # From: http://www.secg.org/sec1-v2.pdf §2.3.8 "Octet-String-to-Integer Conversion"
    if not isinstance(value, int_types):
        raise EIP1829Error("Requires integer type")
    # TODO: conversion from bytes (big-endian)
    return value


def coerce_field_element(m, value):
    value = coerce_integer(value)
    if value < 0 or value >= m:
        raise EIP1829Error("Outside of field")
    return FQ(value, m)


def coerce_octet(value):
    if not isinstance(value, int_types):
        raise TypeError("Requires integer type")
    if (value & 0xFF) != value:
        raise EIP1829Error("Not an octet")
    return value


def coerce_Y(X, Y, a, b, p):
    # - 1. If M = 0, output P = O
    if Y == 0:
        Y = FQ(0, p)
        # XXX: is Infinity always (0,0) for Short Weierstrass curves?
        if X != 0:
            raise EIP1829Error("Must specify X coordinate as zero for infinity")
        return Y

    #  - 2.3. If Y = 2, set ỹ = 0, and if Y = 3, set ỹ = 1.
    #         Otherwise output “invalid” and stop.
    if Y in [2, 3]:
        ytilda = Y % 2

        #  - 2.4.1. If q = p is an odd prime, compute the field element α ≡ x^3 + ax + b (mod p)
        #           and compute a square root β of α modulo p.
        #           Output “invalid” and stop if there are no square roots of α modulo p
        Y = recover_y(X, a, b)

        # XXX: Different implementations may give different results, based on the sqrt method used
        # TODO: normalise the value of the result of the square root in `recover_y`?

        #           otherwise set y = β if β ≡ ỹ (mod 2), and set y = p − β if β !≡ ỹ (mod 2).
        if (int(Y) % 2) != ytilda:
            Y = -Y

        return Y

    raise EIP1829Error("Invalid value for Y coordinate signdness (%r)" % (Y,))


def eip1829(p, a, b, *args):
    """
    https://github.com/ethereum/EIPs/blob/af94e5498b9db0ee32f96bdba199a76b138906bc/EIPS/eip-1829.md

    C = s₀ ⋅ A₀ + s₁ ⋅ A₁ + ⋯ + s_n ⋅ A_n

    aka linear combination, inner product, multi-multiplication or even multi-exponentiation.

    (Cx, Cy) := ecmul(p, α, β, s0, Ax0, As0, s1, Ax1, As1, ...)
    """
    if len(args) % 3 != 0:
        raise EIP1829Error('Number of arguments must be divisible by 3')

    # field modulus for F_p
    p = coerce_integer(p)
    if p <= 2:
        raise EIP1829Error("Invalid Modulus")

    # Parameters `a` and `b` of the curve equation
    a = coerce_field_element(p, a)
    b = coerce_field_element(p, b)

    # Equation
    Cx, Cy = None, None

    # Iterate arguments, 3 elements at a time; (s, X_i, Y_i)
    it = iter(args)
    for s in it:
        X = next(it)
        Y = next(it)

        # From: http://www.secg.org/sec1-v2.pdf §2.3.4 "Octet-String-to-Elliptic-Curve-Point Conversion"
        #  - 2.2. Convert X to a field element xP of Fq using the conversion routine specified in Section 2.3.6.
        #         Output “invalid” and stop if the routine outputs “invalid”.
        X = coerce_field_element(p, X)

        # Then recover the Y coordinate from the X, where the input Y is a flag with special handling
        Y = coerce_Y(X, Y, a, b, p)

        # Multiply point by the scalar, and accumulate the result
        mX, mY = multiply(X, Y, s, a, b)
        Cx, Cy = addition(Cx, Cy, mX, mY, a, b)

    return (Cx, Cy)
