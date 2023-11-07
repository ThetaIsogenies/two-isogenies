# ============================================ #
#     Fast square root and quadratic roots     #
# ============================================ #


def canonical_root(a):
    """
    Very stupid and slow way, but it
    makes the sqrt match rust for all
    cases
    """
    a0, a1 = a.list()
    if a0.is_zero() and (int(a1) % 2) == 1:
        return -a
    if (int(a0) % 2) == 1:
        return -a
    return a


def invert_or_zero(x):
    """
    Used to match the rust, which returns the inverse of 0 as 0
    """
    if x == 0:
        return 0
    return 1 / x


def sqrt_Fp(x):
    """
    Faster computation of sqrt in Fp assuming p = 3 mod 4
    """
    p = x.parent().characteristic()
    exp = (p + 1) // 4

    r = x**exp
    if r * r != x:
        return 0

    # To ensure the result matches
    # Rust
    if int(r) % 2 != 0:
        return -r
    return r


def sqrt_Fp2(x, canonical=False):
    """
    Fast computation of square-roots in SageMath using that p = 3 mod 4
    """
    F = x.parent()
    x0, x1 = x.list()

    if x1 == 0:
        lx0 = x0.is_square()
        if lx0:
            y0 = sqrt_Fp(x0)
            root = F([y0, 0])
            if canonical:
                return canonical_root(root)
            return root

        else:
            y1 = sqrt_Fp(-x0)
            root = F([0, y1])
            if canonical:
                return canonical_root(root)
            return root

    delta = x0**2 + x1**2
    sqrt_delta = sqrt_Fp(delta)

    y02 = (x0 + sqrt_delta) / 2
    if not y02.is_square():
        y02 -= sqrt_delta

    y0 = sqrt_Fp(y02)
    y1 = x1 / (y0 + y0)
    root = F([y0, y1])

    if canonical:
        return canonical_root(root)
    return root
