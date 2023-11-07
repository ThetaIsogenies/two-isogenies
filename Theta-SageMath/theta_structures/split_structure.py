from theta_structures.dimension_two import ThetaStructure, ThetaPoint
from theta_structures.dimension_one import (
    theta_null_point_to_montgomery_curve,
    theta_point_to_montgomery_point,
)
from theta_structures.couple_point import CouplePoint
from utilities.fast_sqrt import sqrt_Fp2


class SplitThetaStructure:
    """
    Given some ThetaStrcuture which is isomorphic to a product of elliptic
    curves E1 x E2, this class takes as input a ThetaStructure and returns a new
    ThetaStrcuture which has a compatible representation: the dual, squared
    coordinate is in the index position (11, 11). In this form, computing the
    dimension one theta points (a : b) from (X : Y : Z : W) is simply a case of
    selecting the splitting (X : Y) and (Y : W)

    Contains helper functions which return the elliptic products E1 and E2 and also
    allows the mapping of some ThetaPoint P on A to the Elliptic Product E1 x E2.
    """

    def __init__(self, T):
        if not isinstance(T, ThetaStructure):
            raise TypeError

        # Create dim 1 theta structures
        self.O1, self.O2 = self.split(T)

        # Compute Mont. curves
        self.E1 = theta_null_point_to_montgomery_curve(self.O1)
        self.E2 = theta_null_point_to_montgomery_curve(self.O2)

    def curves(self):
        """
        Returns the elliptic curves E1 and E2 in the Montgomery model using the
        formula of

            Models of Kummer lines and Galois representation,
            Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        """
        return self.E1, self.E2

    @staticmethod
    def split(P):
        """
        Assuming the zero index of the ThetaStructure is in the (11, 11)
        postion, we can always obtain the dim-1 theta structure in the following
        way

        P = (a : b : c : d)

        P1 = (a : b)
        P2 = (b : d)

        This is because the product structure of P1 and P2 is simply:

        P = (a * b : b * b : a * b : b * d)

        And we see up to an overall scale factor we recover the projective
        factors essentially for free
        """
        a, b, _, d = P.coords()

        P1 = (a, b)
        P2 = (b, d)

        return P1, P2

    @staticmethod
    def to_points(E, X, Z):
        """
        Given the (X : Z) point on the KummerLine of E compute the point
            ±P = (X : Y : Z) on the curve
        """
        if Z == 0:
            return E(0)

        x = X / Z

        A = E.a_invariants()[1]
        y2 = x * (x**2 + A * x + 1)
        y = sqrt_Fp2(y2)

        return E(x, y)

    def __call__(self, P, lift=True):
        """ """
        if not isinstance(P, ThetaPoint):
            raise TypeError

        # Dim 2 -> Dim 1 theta points
        P1, P2 = self.split(P)

        # Convert to Montgomery points
        Q1X, Q1Z = theta_point_to_montgomery_point(self.O1, P1)
        Q2X, Q2Z = theta_point_to_montgomery_point(self.O2, P2)

        if lift:
            # lift from the Kummer to the elliptic curve
            Q1 = self.to_points(self.E1, Q1X, Q1Z)
            Q2 = self.to_points(self.E2, Q2X, Q2Z)

            return CouplePoint(Q1, Q2)
        else:
            return [(Q1X, Q1Z), (Q2X, Q2Z)]
