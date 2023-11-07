from theta_structures.split_structure import SplitThetaStructure
from theta_structures.couple_point import CouplePoint
from theta_isogenies.morphism import Morphism
from theta_isogenies.gluing_isogeny import GluingThetaIsogeny
from theta_isogenies.isomorphism import SplittingIsomorphism
from theta_isogenies.isogeny import ThetaIsogeny
from utilities.strategy import optimised_strategy


class EllipticProductIsogeny(Morphism):
    r"""
    Given (P1, P2), (Q1, Q2) in (E1 x E2)[2^(n+2)] as the generators of a kernel
    of a (2^n, 2^n)-isogeny between elliptic products, computes this isogeny

    ker(Phi) = <(P1, P2), (Q1, Q2)> : E1 x E2 -> E3 x E4

    Input:

    - kernel = CouplePoint(P1, P2), CouplePoint(Q1, Q2):
      where points are on the elliptic curves E1, E2 of order 2^(n+2)
    - n: the length of the chain
    - strategy: the optimises strategy to compute a walk through the graph of
      images and doublings with a quasli-linear number of steps
    - zeta (optional): a second root of unity

    NOTE: if only the 2^n torsion is known, the isogeny should be computed with
    `EllipticProductIsogenySqrt()` which computes the last two steps without the
    torsion above the kernel, but instead with square-root computations (which
    is slower)
    """

    def __init__(self, kernel, n, strategy=None, zeta=None):
        self.n = n
        self.E1, self.E2 = kernel[0].curves()
        self._zeta = zeta
        assert kernel[1].curves() == (self.E1, self.E2)

        self._domain = (self.E1, self.E2)

        if strategy is None:
            strategy = self.get_strategy()
        self.strategy = strategy

        self._phis = self.isogeny_chain(kernel)
        T_last = self._phis[-1].codomain()
        self._splitting = SplitThetaStructure(T_last)

        self._codomain = self._splitting.curves()

    def get_strategy(self):
        return optimised_strategy(self.n)

    def isogeny_chain(self, kernel):
        """
        Compute the codomain of the isogeny chain and store intermediate
        isogenies for evaluation
        """
        # Extract the CouplePoints from the Kernel
        Tp1, Tp2 = kernel

        # Store chain of (2,2)-isogenies
        isogeny_chain = []

        # Bookkeeping for optimal strategy
        strat_idx = 0
        level = [0]
        ker = (Tp1, Tp2)
        kernel_elements = [ker]

        for k in range(self.n):
            prev = sum(level)
            ker = kernel_elements[-1]

            while prev != (self.n - 1 - k):
                level.append(self.strategy[strat_idx])

                # Perform the doublings
                Tp1 = ker[0].double_iter(self.strategy[strat_idx])
                Tp2 = ker[1].double_iter(self.strategy[strat_idx])

                ker = (Tp1, Tp2)

                # Update kernel elements and bookkeeping variables
                kernel_elements.append(ker)
                prev += self.strategy[strat_idx]
                strat_idx += 1

            # Compute the codomain from the 8-torsion
            Tp1, Tp2 = ker
            if k == 0:
                phi = GluingThetaIsogeny(Tp1, Tp2)
            elif k == self.n - 2:
                # The next isogeny will be a splitting isogeny, so we know we
                # will have one of a,b,c,d = 0. So at this point switch to
                # dual theta coordinate
                phi = ThetaIsogeny(Th, Tp1, Tp2, hadamard=(False, False))
            elif k == self.n - 1:
                # Compute the dual isogeny, remembering that we switched to
                # dual theta coordinates at the previous step.
                # We output dual theta coordinates on the product, change
                # to hadamard=(True, True) to output standard coordinates;
                # this does not change the conversion back to Montgomery
                # coordinates so we might as well save an Hadamard
                # transform anyway
                phi = ThetaIsogeny(Th, Tp1, Tp2, hadamard=(True, False))
            else:
                phi = ThetaIsogeny(Th, Tp1, Tp2)

            # Update the chain of isogenies
            Th = phi.codomain()
            isogeny_chain.append(phi)

            # Remove elements from list
            kernel_elements.pop()
            level.pop()

            # Push through points for the next step
            kernel_elements = [(phi(T1), phi(T2)) for T1, T2 in kernel_elements]

        splitting_iso = SplittingIsomorphism(Th, zeta=self._zeta)
        isogeny_chain.append(splitting_iso)

        return isogeny_chain

    def evaluate_isogeny(self, P):
        """
        Given a point P, of type CouplePoint on the domain E1 x E2, computes the
        CouplePoint on the codomain ThetaStructure isomorphic to E3 x E4
        """
        if not isinstance(P, CouplePoint):
            raise TypeError(
                "EllipticProductIsogeny isogeny expects as input a CouplePoint on the domain product E1 x E2"
            )
        for f in self._phis:
            P = f(P)
        return P

    def __call__(self, P, lift=True):
        """
        Evaluate a CouplePoint under the action of this isogeny. If lift=True,
        then the affine coordinates of the points are returned, otherwise points
        on the Kummer line are returned.
        """
        image_P = self.evaluate_isogeny(P)
        return self._splitting(image_P, lift=lift)
