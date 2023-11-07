import unittest

from theta_structures.dimension_one import *
from theta_structures.dimension_two import *
from theta_structures.product_structure import ProductThetaStructure
from theta_structures.couple_point import *
from isogeny_diamond import *

from tests.test_utils import random_supersingular_curve, random_supersingular_curves


class DimensionOne(unittest.TestCase):
    def test_conversion(self):
        for _ in range(10):
            # Pick a random supersingular curve
            E = random_supersingular_curve()

            # Compute corresponding theta null point
            O0 = montgomery_curve_to_theta_null_point(E)

            # Ensure the conversion back works
            self.assertEqual(theta_null_point_to_montgomery_curve(O0), E)

            for _ in range(10):
                # Test conversion to and from Montgomery points
                P = E.random_point()
                if P.is_zero():
                    continue

                O = montgomery_point_to_theta_point(O0, P)
                X, Z = theta_point_to_montgomery_point(O0, O)

                self.assertEqual(P[0], X / Z)


class DimensionTwo(unittest.TestCase):
    def test_double_iter(self):
        for _ in range(10):
            E1, E2 = random_supersingular_curves()
            T = ProductThetaStructure(E1, E2)

            for _ in range(20):
                # Pick some random points on each curve
                P1 = E1.random_point()
                P2 = E2.random_point()

                # Lift these pairs of points onto T
                k = randint(1, 100)
                O1 = T(P1, P2)
                O2 = T(2**k * P1, 2**k * P2)

                # Doubling points then lifting is the same as lifting then doubling
                self.assertEqual(O1.double_iter(k), O2)

    def test_doubling(self):
        for _ in range(10):
            E1, E2 = random_supersingular_curves()
            T = ProductThetaStructure(E1, E2)

            for _ in range(20):
                # Pick some random points on each curve
                P1 = E1.random_point()
                P2 = E2.random_point()

                # Lift these pairs of points onto T
                O1 = T(P1, P2)
                O2 = T(2 * P1, 2 * P2)

                # Doubling points then lifting is the same as lifting then doubling
                self.assertEqual(2 * O1, O2)

    def test_multiplication(self):
        for _ in range(10):
            E1, E2 = random_supersingular_curves()
            T = ProductThetaStructure(E1, E2)
            p = E1.base_ring().characteristic()

            for _ in range(20):
                # Pick some random points on each curve
                P1 = E1.random_point()
                P2 = E2.random_point()

                k = 0
                while k == 0:
                    k = randint(-p, p)

                # Lift these pairs of points onto T
                O1 = T(P1, P2)
                O2 = T(k * P1, k * P2)

                # Multiplying points then lifting is the same as lifting then multiplying
                self.assertEqual(k * O1, O2)


if __name__ == "__main__" and "__file__" in globals():
    unittest.main()
