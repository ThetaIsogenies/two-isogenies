from sage.all import ZZ

from theta_structures.dimension_two import ThetaStructure, ThetaPoint
from theta_isogenies.morphism import Morphism
from utilities.batched_inversion import batched_inversion


class ThetaIsogeny(Morphism):
    def __init__(self, domain, T1_8, T2_8, hadamard=(False, True)):
        """
        Compute a (2,2)-isogeny in the theta model. Expects as input:

        - domain: the ThetaStructure from which we compute the isogeny
        - (T1_8, T2_8): points of 8-torsion above the kernel generating the isogeny

        When the 8-torsion is not available (for example at the end of a long
        (2,2)-isogeny chain), the the helper functions in isogeny_sqrt.py
        must be used.

        NOTE: on the hadamard bools:

        The optional parameter 'hadamard' controls if we are in standard or dual
        coordinates, and if the codomain is in standard or dual coordinates. By
        default this is (False, True), meaning we use standard coordinates on
        the domain A and the codomain B.

        The kernel is then the kernel K_2 where the action is by sign. Other
        possibilities: - (False, False): standard coordinates on A, dual
        coordinates on B - (True, True): start in dual coordinates on A
        (alternatively: standard coordinates on A but quotient by K_1 whose
        action is by permutation), and standard coordinates on B. - (True,
        False): dual coordinates on A and B

        These can be composed as follows for A -> B -> C:

        - (False, True) -> (False, True) (False, False) -> (True, True):
          - standard coordinates on A and C,
          - standard/resp dual coordinates on B
        - (False, True) -> (False, False) (False, False) -> (True, False):
          - standard coordinates on A,
          - dual coordinates on C,
          - standard/resp dual coordinates on B
        - (True, True) -> (False, True) (True, False) -> (True, True):
          - dual coordinates on A,
          - standard coordinates on C,
          - standard/resp dual coordiantes on B
        - (True, True) -> (False, False) (True, False) -> (True, False):
          - dual coordinates on A and C
          - standard/resp dual coordinates on B

        On the other hand, these gives the multiplication by [2] on A:

        - (False, False) -> (False, True) (False, True) -> (True, True):
          - doubling in standard coordinates on A
          - going through dual/standard coordinates on B=A/K_2
        - (True, False) -> (False, False) (True, True) -> (True, False):
          - doubling in dual coordinates on A
          - going through dual/standard coordinates on B=A/K_2
            (alternatively: doubling in standard coordinates on A going
            through B'=A/K_1)
        - (False, False) -> (False, False) (False, True) -> (True, False):
          - doubling from standard to dual coordinates on A
        - (True, False) -> (False, True) (True, True) -> (True, True):
          - doubling from dual to standard coordinates on A
        """
        if not isinstance(domain, ThetaStructure):
            raise ValueError
        self._domain = domain

        self._hadamard = hadamard
        self._precomputation = None
        self._codomain = self._compute_codomain(T1_8, T2_8)

    def _compute_codomain(self, T1, T2):
        """
        Given two isotropic points of 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2
        """
        if self._hadamard[0]:
            xA, xB, _, _ = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*T1.coords())
            )
            zA, tB, zC, tD = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*T2.coords())
            )
        else:
            xA, xB, _, _ = T1.squared_theta()
            zA, tB, zC, tD = T2.squared_theta()

        if not self._hadamard[0] and self._domain._precomputation:
            # Batch invert denominators
            xA_inv, zA_inv, tB_inv = batched_inversion(xA, zA, tB)

            # Compute A, B, C, D
            A = ZZ(1)
            B = xB * xA_inv
            C = zC * zA_inv
            D = tD * tB_inv * B

            _, _, _, BBinv, CCinv, DDinv = self._domain._arithmetic_precomputation()
            B_inv = BBinv * B
            C_inv = CCinv * C
            D_inv = DDinv * D
        else:
            # Batch invert denominators
            xA_inv, zA_inv, tB_inv, xB_inv, zC_inv, tD_inv = batched_inversion(
                xA, zA, tB, xB, zC, tD
            )

            # Compute A, B, C, D
            A = ZZ(1)
            B = xB * xA_inv
            C = zC * zA_inv
            D = tD * tB_inv * B
            B_inv = xB_inv * xA
            C_inv = zC_inv * zA
            D_inv = tD_inv * tB * B_inv

            # NOTE: some of the computations we did here could be reused for the
            # arithmetic precomputations of the codomain However, we are always
            # in the mode (False, True) except the very last isogeny, so we do
            # not lose much by not doing this optimisation Just in case we need
            # it later:
            # - for hadamard=(False, True): we can reuse the arithmetic
            #   precomputation; we do this already above
            # - for hadamard=(False, False): we can reuse the arithmetic
            #   precomputation as above, and furthermore we could reuse B_inv,
            #   C_inv, D_inv for the precomputation of the codomain
            # - for hadamard=(True, False): we could reuse B_inv, C_inv, D_inv
            #   for the precomputation of the codomain
            # - for hadamard=(True, True): nothing to reuse!

        self._precomputation = (B_inv, C_inv, D_inv)
        if self._hadamard[1]:
            a, b, c, d = ThetaPoint.to_hadamard(A, B, C, D)
            return ThetaStructure([a, b, c, d])
        else:
            return ThetaStructure([A, B, C, D])

    def __call__(self, P):
        """
        Take into input the theta null point of A/K_2, and return the image
        of the point by the isogeny
        """
        if not isinstance(P, ThetaPoint):
            raise TypeError("Isogeny evaluation expects a CouplePoint as input")

        if self._hadamard[0]:
            xx, yy, zz, tt = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*P.coords())
            )
        else:
            xx, yy, zz, tt = P.squared_theta()

        Bi, Ci, Di = self._precomputation

        yy = yy * Bi
        zz = zz * Ci
        tt = tt * Di

        image_coords = (xx, yy, zz, tt)
        if self._hadamard[1]:
            image_coords = ThetaPoint.to_hadamard(*image_coords)
        return self._codomain(image_coords)
