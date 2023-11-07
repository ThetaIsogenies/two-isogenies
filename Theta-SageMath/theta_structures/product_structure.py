# Local imports
from theta_structures.dimension_one import (
    montgomery_curve_to_theta_null_point,
    montgomery_point_to_theta_point,
    montgomery_torsion_to_theta_null_point,
)
from theta_structures.dimension_two import ThetaStructure
from theta_structures.couple_point import CouplePoint


class ProductThetaStructure(ThetaStructure):
    """
    Special case of a ThetaStructure, which is computed from two elliptic curves
    E1, E2 interpreted as a product: A = E1 x E2.

    The bulk of the work here is the mapping from the dimension one structure to
    the product dimension two structure. Bulk of the conversions come from

        Models of Kummer lines and Galois representation,
        Razvan Barbulescu, Damien Robert and Nicolas Sarkis
    """

    def __init__(self, E1, E2, B1=None, B2=None):
        if B1 == None or B2 == None:
            # null point is dim 2 product
            # O01, O02 are the dim 1 null points
            null_point, O01, O02 = self.montgomery_curves_to_theta_null_point(E1, E2)

            self.E1 = E1
            self.E2 = E2

            self.O01 = O01
            self.O02 = O02

            ThetaStructure.__init__(self, null_point)
        else:
            self.E1 = E1
            self.E2 = E2

            self.O01 = montgomery_torsion_to_theta_null_point(B1[0])
            self.O02 = montgomery_torsion_to_theta_null_point(B2[0])

            null_point = self.pairwise_product(self.O01, self.O02)

            ThetaStructure.__init__(self, null_point)

    def __call__(self, *args):
        """
        Create a theta point on the product from either a pair of
        points, a CouplePoint or coefficients of a theta point.
        """
        if len(args) == 2:
            P1, P2 = args
            coords = self.montgomery_points_to_theta_point(P1, P2)
            return self._point(self, coords)

        elif len(args) == 1:
            input = args[0]
            if isinstance(input, CouplePoint):
                P1, P2 = input.points()
                coords = self.montgomery_points_to_theta_point(P1, P2)
                return self._point(self, coords)

        return self._point(self, *args)

    def __repr__(self):
        return f"Product theta structure over {self.base_ring()} with null point: {self.null_point()}"

    def pairwise_product(self, O1, O2):
        """
        Given two dim-1 theta points, compute the level two theta point on the
        product.

        NOTE: The order used is implicitly 00, 10, 01, 11
        """
        a1, b1 = O1
        a2, b2 = O2

        return (a1 * a2, b1 * a2, a1 * b2, b1 * b2)

    def montgomery_curves_to_theta_null_point(self, E1, E2):
        """
        Given a pair of elliptic curves, compute the corresponding
        dim-1 theta points as well as the coordinates of the
        level-2 null point.
        """
        O1 = montgomery_curve_to_theta_null_point(E1)
        O2 = montgomery_curve_to_theta_null_point(E2)

        null_coords = self.pairwise_product(O1, O2)

        return null_coords, O1, O2

    def montgomery_points_to_theta_point(self, P1, P2):
        """
        Given two points P1, P2 on curves E1, E2 where it is
        assumed we have computed the corresponding dim-1
        theta null points, compute the product theta point
        from P1 x P2.
        """
        O1 = montgomery_point_to_theta_point(self.O01, P1)
        O2 = montgomery_point_to_theta_point(self.O02, P2)

        return self.pairwise_product(O1, O2)
