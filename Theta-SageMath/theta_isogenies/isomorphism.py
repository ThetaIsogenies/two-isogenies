from sage.all import Matrix
from theta_structures.dimension_two import ThetaStructure
from theta_structures.dimension_two import ThetaPoint
from theta_structures.couple_point import CouplePoint
from theta_isogenies.morphism import Morphism


class Isomorphism(Morphism):
    """
    Given a 4x4 matrix, compute an isomorphism described by this change of basis
    matrix
    """

    def __init__(self, N=None):
        self.N = N

    def dual(self):
        """ """
        if self.N is None:
            raise ValueError("Cannot compute the dual without the matrix N")

        if self.domain() is None or self.codomain() is None:
            raise ValueError(
                "Dual can only be computed once domain and codomain are known"
            )

        return DualIsomorphism(self.domain(), self.codomain(), N)

    def apply_isomorphism(self, P):
        """
        Applies change of theta coordinate to a point.

        Input:
        - N: matrix of the level 2 coordinates base change (in row convention).
        - P: point in level 2 theta-coordinates.

        Output: Point in level 2 theta'-coordinates after base change by N.
        """
        if self.N is None:
            raise ValueError(
                "Cannot compute an isomorphism with the corresponding matrix set."
            )

        x, y, z, t = P.coords()
        X = self.N[0, 0] * x + self.N[0, 1] * y + self.N[0, 2] * z + self.N[0, 3] * t
        Y = self.N[1, 0] * x + self.N[1, 1] * y + self.N[1, 2] * z + self.N[1, 3] * t
        Z = self.N[2, 0] * x + self.N[2, 1] * y + self.N[2, 2] * z + self.N[2, 3] * t
        T = self.N[3, 0] * x + self.N[3, 1] * y + self.N[3, 2] * z + self.N[3, 3] * t

        return (X, Y, Z, T)

    def __call__(self, Q):
        if not isinstance(Q, (ThetaPoint, CouplePoint)):
            raise TypeError(
                f"Cannot compute isomorphism on input: {Q} of type {type(Q)}"
            )

        if isinstance(Q, CouplePoint):
            # Convert to a theta point first
            Q = self.domain()(Q)

        new_coords = self.apply_isomorphism(Q)
        return self._codomain(new_coords)


class DualIsomorphism(Isomorphism):
    """
    Given a change of basis matrix N mapping from a domain to a codomain,
    compute the dual isomorphism by inverting N and switching the (co)domains.

    NOTE: assumes N is invertible.
    """

    def __init__(self, domain, codomain, N):
        self._domain = codomain
        self._codomain = domain
        self.N = N.inverse()


class SplittingIsomorphism(Isomorphism):
    """
    Given a ThetaStructure which admits a splitting, compute the isomorphism
    which sets the zero U^2 value to appear at the index (11, 11)"""

    def __init__(self, domain, zeta=None):
        # Ensure domain is of the correct type
        if not isinstance(domain, ThetaStructure):
            raise ValueError("Domain must be a Theta Structure")
        self._domain = domain

        # We need a second root of unity
        if zeta is None:
            self.zeta = domain.base_ring().gen()

        # Select the precomputed change of basis to map the zero to (11, 11) by
        # identifying the current zero index.
        self.N = self.compute_splitting_matrix()

        # Compute the codomain
        mapped_null_coords = self.apply_isomorphism(domain.null_point())

        # Set the codomain from the new null values
        self._codomain = ThetaStructure(mapped_null_coords)

    def level_22_constants_sqr(self, null_coords, chi, i):
        """
        For a given index (chi, i) compute the U^2 value. Used to find the zero value.
        """
        chi_eval = {
            (0, 0): 1,
            (0, 1): 1,
            (0, 2): 1,
            (0, 3): 1,
            (1, 0): 1,
            (1, 1): -1,
            (1, 2): 1,
            (1, 3): -1,
            (2, 0): 1,
            (2, 1): 1,
            (2, 2): -1,
            (2, 3): -1,
            (3, 0): 1,
            (3, 1): -1,
            (3, 2): -1,
            (3, 3): 1,
        }
        U_constant = 0
        for t in range(4):
            U_constant += chi_eval[(chi, t)] * null_coords[t] * null_coords[i ^ t]
        return U_constant

    def identify_even_index(self, null_coords):
        """
        Loop through all even index values to find (chi, i) such that U^2 is
        zero.
        """
        even_indices = [
            [0, 0],
            [0, 1],
            [0, 2],
            [0, 3],
            [1, 0],
            [1, 2],
            [2, 0],
            [2, 1],
            [3, 0],
            [3, 3],
        ]
        for chi, i in even_indices:
            if self.level_22_constants_sqr(null_coords, chi, i) == 0:
                return chi, i
        raise ValueError("Not a product of elliptic curves")

    def compute_splitting_matrix(self):
        """
        Find the zero value and from this select the right change of basis
        """

        null_coords = self.domain().coords()
        i = self.zeta  # self.domain().base_ring().gen()
        assert i**2 == -1

        # Computed from description in the paper.
        splitting_map = {
            (0, 2): [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
            (3, 3): [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1],
            (0, 3): [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1],
            (2, 1): [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1],
            (0, 1): [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, -1, 0, 0],
            (1, 2): [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
            (2, 0): [1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1],
            (3, 0): [1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1],
            (1, 0): [1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1],
            (0, 0): [1, i, 1, i, 1, -i, -1, i, 1, i, -1, -i, -1, i, -1, i],
        }

        zero_location = self.identify_even_index(null_coords)
        return Matrix(4, 4, splitting_map[zero_location])
