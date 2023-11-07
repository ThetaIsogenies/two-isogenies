"""
Helper functions for the supersingular elliptic curve computations in FESTA
"""

# Sage Imports
from sage.all import ZZ

# Local imports
from utilities.order import has_order_D
from utilities.discrete_log import weil_pairing_pari
from utilities.fast_sqrt import sqrt_Fp2

# =========================================== #
#   Extract coefficent from Montgomery curve  #
# =========================================== #


def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError("The elliptic curve E is not in the Montgomery model.")
    return A


# =========================================== #
# Compute points of order D and Torsion Bases #
# =========================================== #


def random_point(E):
    """
    Returns a random point on the elliptic curve E
    assumed to be in Montgomery form with a base
    field which characteristic p = 3 mod 4
    """
    A = E.a_invariants()[1]
    if E.a_invariants() != (0, A, 0, 1, 0):
        raise ValueError(
            "Function `generate_point` assumes the curve E is in the Montgomery model"
        )

    # Try 1000 times then give up, just protection
    # for infinite loops
    F = E.base_ring()
    for _ in range(1000):
        x = F.random_element()
        y2 = x * (x**2 + A * x + 1)
        if y2.is_square():
            y = sqrt_Fp2(y2)
            return E(x, y)

    raise ValueError(
        "Generated 1000 points, something is probably going wrong somewhere."
    )


def generate_point(E, x_start=0):
    """
    Generate points on a curve E with x-coordinate
    i + x for x in Fp and i is the generator of Fp^2
    such that i^2 = -1.
    """
    F = E.base_ring()
    one = F.one()

    if x_start:
        x = x_start + one
    else:
        x = F.gen() + one

    A = E.a_invariants()[1]
    if E.a_invariants() != (0, A, 0, 1, 0):
        raise ValueError(
            "Function `generate_point` assumes the curve E is in the Montgomery model"
        )

    # Try 1000 times then give up, just protection
    # for infinite loops
    for _ in range(1000):
        y2 = x * (x**2 + A * x + 1)
        if y2.is_square():
            y = sqrt_Fp2(y2)
            yield E(x, y)
        x += one

    raise ValueError(
        "Generated 1000 points, something is probably going wrong somewhere."
    )


def generate_point_order_D(E, D, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            An integer D dividing (p +1)
    Output: A point P of order D.
    """
    p = E.base().characteristic()
    n = (p + 1) // D

    Ps = generate_point(E, x_start=x_start)
    for G in Ps:
        P = n * G

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            P._order = ZZ(D)
            yield P

    raise ValueError(f"Never found a point P of order D.")


def compute_point_order_D(E, D, x_start=0):
    """
    Wrapper function around a generator which returns the first
    point of order D
    """
    return generate_point_order_D(E, D, x_start=x_start).__next__()


def compute_linearly_independent_point_with_pairing(E, P, D, x_start=0):
    """
    Input:  An elliptic curve E / Fp2
            A point P ∈ E[D]
            An integer D dividing (p +1)
    Output: A point Q such that E[D] = <P, Q>
            The Weil pairing e(P,Q)
    """
    Qs = generate_point_order_D(E, D, x_start=x_start)
    for Q in Qs:
        # Make sure the point is linearly independent
        pair = weil_pairing_pari(P, Q, D)
        if has_order_D(pair, D, multiplicative=True):
            Q._order = ZZ(D)
            return Q, pair
    raise ValueError("Never found a point Q linearly independent to P")


def compute_linearly_independent_point(E, P, D, x_start=0):
    """
    Wrapper function around `compute_linearly_independent_point_with_pairing`
    which only returns a linearly independent point
    """
    Q, _ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=x_start)
    return Q


def torsion_basis_with_pairing(E, D):
    """
    Generate basis of E(Fp^2)[D] of supersingular curve

    While computing E[D] = <P, Q> we naturally compute the
    Weil pairing e(P,Q), which we also return as in some cases
    the Weil pairing is then used when solving the BiDLP
    """
    Fp2 = E.base()
    p = Fp2.characteristic()
    i = Fp2.gen()

    # Ensure D divides the curve's order
    if (p + 1) % D != 0:
        print(f"{ZZ(D).factor() = }")
        print(f"{ZZ(p+1).factor() = }")
        raise ValueError(f"D must divide the point's order")

    P = compute_point_order_D(E, D)
    x_start = P[0] + i
    Q, ePQ = compute_linearly_independent_point_with_pairing(E, P, D, x_start=x_start)

    return P, Q, ePQ


def torsion_basis(E, D):
    """
    Wrapper function around torsion_basis_with_pairing which only
    returns the torsion basis <P,Q> = E[D]
    """
    P, Q, _ = torsion_basis_with_pairing(E, D)
    return P, Q


# =============================================== #
#  Ensure Basis <P,Q> of E[2^k] has (0,0) under Q #
# =============================================== #


def fix_torsion_basis_renes(P, Q, k):
    """
    Set the torsion basis P,Q such that
    2^(k-1)Q = (0,0) to ensure that (0,0)
    is never a kernel of a two isogeny
    """
    cofactor = 2 ** (k - 1)

    R = cofactor * P
    if R[0] == 0:
        return Q, P
    R = cofactor * Q
    if R[0] == 0:
        return P, Q
    return P, P + Q
