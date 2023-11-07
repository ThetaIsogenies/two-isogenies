from sage.all import ZZ

from theta_structures.dimension_two import ThetaStructure, ThetaPoint
from theta_isogenies.isogeny import ThetaIsogeny
from utilities.batched_inversion import batched_inversion
from utilities.fast_sqrt import sqrt_Fp2


class ThetaIsogeny4(ThetaIsogeny):
    """
    Compute the penultimate (2,2)-isogeny when the eight-torsion above the
    kernel is not known but the four torsion above the sqrt is known. The image
    computation is the same, but the codomain is altered and B, C are computed
    from square-roots.
    """

    def __init__(self, domain, T1_4, T2_4, hadamard=(None, True)):
        super().__init__(domain, T1_4, T2_4, hadamard=hadamard)

    def _compute_codomain(self, T1, T2):
        """
        Given two isotropic points of 4-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2

        As the 8-torsion above the kernel is not known, B and C are computed
        by taking square-roots of B^2 and C^2
        """
        if self._hadamard[0]:
            AA, BB, CC, DD = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*self._domain.coords())
            )
            xAB, _, xCD, _ = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*T1.coords())
            )
        else:
            AA, BB, CC, DD = self._domain.squared_theta()
            xAB, _, xCD, _ = T1.squared_theta()

        AA_inv, BB_inv, CC_inv, DD_inv = batched_inversion(AA, BB, CC, DD)

        A = ZZ(1)
        B = sqrt_Fp2(BB * AA_inv)
        C = sqrt_Fp2(CC * AA_inv)
        D = xCD * B / (xAB * C)

        B_inv = AA * BB_inv * B
        C_inv = AA * CC_inv * C
        D_inv = AA * DD_inv * D
        self._precomputation = (B_inv, C_inv, D_inv)

        if self._hadamard[1]:
            a, b, c, d = ThetaPoint.to_hadamard(A, B, C, D)
            return ThetaStructure([a, b, c, d])
        else:
            return ThetaStructure([A, B, C, D])


class ThetaIsogeny2(ThetaIsogeny4):
    """
    Compute the final (2,2)-isogeny when the eight-torsion above the kernel is
    not known. The image computation is the same, but the codomain is altered
    and B, C and D are computed from square-roots.
    """

    def __init__(self, domain, T1=None, T2=None, hadamard=(None, True)):
        super().__init__(domain, T1, T2, hadamard=hadamard)

    # NOTE: T1 and T2 are not used here. Everything is computed from the domain
    # itself. We could add an assert check checking they are compatible with our
    # theta null point, but for now these are simply dummy variables.
    def _compute_codomain(self, T1, T2):
        """ """
        if self._hadamard[0]:
            AA, BB, CC, DD = ThetaPoint.to_squared_theta(
                *ThetaPoint.to_hadamard(*self._domain.coords())
            )
        else:
            AA, BB, CC, DD = self._domain.squared_theta()

        AA_inv, BB_inv, CC_inv, DD_inv = batched_inversion(AA, BB, CC, DD)

        A = ZZ(1)
        B = sqrt_Fp2(BB * AA_inv)
        C = sqrt_Fp2(CC * AA_inv)
        D = sqrt_Fp2(DD * AA_inv)
        B_inv = AA * BB_inv * B
        C_inv = AA * CC_inv * C
        D_inv = AA * DD_inv * D
        self._precomputation = (B_inv, C_inv, D_inv)

        if self._hadamard[1]:
            a, b, c, d = ThetaPoint.to_hadamard(A, B, C, D)
            return ThetaStructure([a, b, c, d])
        else:
            return ThetaStructure([A, B, C, D])
