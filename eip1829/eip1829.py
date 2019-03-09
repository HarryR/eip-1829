from .field import FQ, int_types


def doubling(x1, y1, a, b):
    # https://www.hyperelliptic.org/EFD/g1p/auto-shortw.html
    x3 = ((3*(x1**2)+a)**2)/((2*y1)**2)-x1-x1
    y3 = (2*x1+x1)*(3*(x1**2)+a)/(2*y1)-(3*(x1**2)+a)**3/(2*y1)**3-y1
    return x3, y3


def addition(x1, y1, x2, y2, a, b):
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


def is_on_curve(x, y, alpha, beta):
    """
    y² = x³ + α ⋅ x + β  mod  m
    """
    yy = y * y
    rhs = x**3 + (alpha * x) + beta
    return yy == rhs


def recover_y(x, alpha, beta):
    yy = x**3 + (alpha * x) + beta
    return y.sqrt()


def is_valid_field_element(m, value):
    if not isinstance(value, int_types):
        raise TypeError("Invalid type")
    if value < 0 or value >= m:
        raise valueError("Outside of field")
    return True


def eip1829(m, alpha, beta, s, A):
    """
    https://github.com/ethereum/EIPs/blob/af94e5498b9db0ee32f96bdba199a76b138906bc/EIPS/eip-1829.md

    C = s₀ ⋅ A₀ + s₁ ⋅ A₁ + ⋯ + s_n ⋅ A_n

    aka linear combination, inner product, multi-multiplication or even multi-exponentiation.

    (Cx, Cy) := ecmul(m, α, β,  s0, Ax0, As0, s1, Ax1, As1, ...)
    """
    assert isinstance(m, int_types)
    assert m > 2
    assert is_valid_field_element(m, alpha)
    assert is_valid_field_element(m, beta)
    alpha = FQ(alpha, m)
    beta = FQ(beta, m)
    C = (None, None)
    for s_i, A_i in zip(s, A):
        X, Y = FQ(A_i[0], m), FQ(A_i[1], m)
        mX, mY = multiply(X, Y, s_i, alpha, beta)
        C = addition(C[0], C[1], mX, mY, alpha, beta)
    return C
