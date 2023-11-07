from sage.all import EllipticCurve
from collections import namedtuple
from utilities.supersingular import montgomery_coefficient
from utilities.fast_sqrt import sqrt_Fp2

ThetaNullPoint = namedtuple("ThetaNullPoint_dim_1", "a b")


def theta_null_point_to_montgomery_curve(O0):
    """
    Given a level 2 theta null point (a:b), compute a Montgomery curve equation.
    We use the model where the 4-torsion point (1:0) above (a:-b) is sent to
    (1:1) in Montgomery coordinates.

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    a, b = O0

    aa = a**2
    bb = b**2

    T1 = aa + bb
    T2 = aa - bb

    # Montgomery coefficient
    A = -(T1**2 + T2**2) / (T1 * T2)

    # Construct curve
    F = a.parent()
    E = EllipticCurve(F, [0, A, 0, 1, 0])
    return E


def montgomery_curve_to_theta_null_point(E):
    """
    From an elliptic curve in Montgomery form, compute a theta null point
    (a:b).
    Let T1=(1:1) the canonical point of 4-torsion in Montgomery coordinates
    and T2 such that (T1,T2) forms a symplectic basis of E[4].
    There are 4 choices of T2, giving 4 different theta null points.

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    # Extract A from curve equation
    A = montgomery_coefficient(E)

    # alpha is a root of
    # x^2 + Ax + 1
    disc = A * A - 4
    assert disc.is_square()

    d = sqrt_Fp2(A * A - 4)
    alpha = (-A + d) / 2

    # The theta coordinates a^2 and b^2
    # are related to alpha
    aa = alpha + 1
    bb = alpha - 1

    aabb = aa * bb
    ab = sqrt_Fp2(aabb)

    # We aren't given (a,b) rational, but
    # (a/b) is rational so we use
    # (ab : b^2) as the theta null point
    O0 = ThetaNullPoint(ab, bb)
    return O0


def montgomery_torsion_to_theta_null_point(P):
    """
    From a four torsion point T2 on a Montgomery curve, such that
    (T1,T2) is a symplectic basis of E[4] and T1=(1:1),
    compute the corresponding theta null point.

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    r = P[0]
    s = P[2]
    O0 = ThetaNullPoint(r + s, r - s)
    return O0


def theta_point_to_montgomery_point(O0, O):
    """
    Given an elliptic curve in Montgomery form
       E : y^2 = x^3 + Ax^2 + x
    and the theta null point with coordinates
       O0 = (a, b)

    Converts a theta point O = (U : V) on θ_(a,b)
    to a point P = (X : Z) on E

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    a, b = O0
    U, V = O
    X = a * V + b * U
    Z = a * V - b * U

    return (X, Z)


def montgomery_point_to_theta_point(O0, P):
    """
    Given an elliptic curve in Montgomery form
       E : y^2 = x^3 + Ax^2 + x
    and the theta null point with coordinates
       (a, b)

    Converts a point P = (X : Z) on E to a point
    O = (U : V) on θ_(a,b)

    Algorithm from:
        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """
    if P.is_zero():
        return O0[:]

    a, b = O0
    X, Z = P[0], P[2]
    if X == Z and Z == 0:  # Do not forget the origin
        return (a, b)
    else:
        U = a * (X - Z)
        V = b * (X + Z)
        return (U, V)
