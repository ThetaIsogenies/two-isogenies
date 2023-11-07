"""
For testing, we need an easy way to generate a kernel which generates and
isogeny between elliptic products. This file, for a range of characteristics,
does exactly this by picking primes where we can compute (2^n, 2^n)-isogenies
between elliptic products by using that there is an endomorphism of degree 
    A - B = X^2 + Y^2
"""

from sage.all import *

from richelot_isogenies.richelot_isogenies import (
    compute_richelot_chain,
    evaluate_richelot_chain,
)
from montgomery_isogenies.isogenies_x_only import (
    isogeny_from_scalar_x_only,
    evaluate_isogeny_x_only,
)
from utilities.supersingular import torsion_basis, fix_torsion_basis_renes
from utilities.order import has_order_D
from utilities.discrete_log import BiDLP
from utilities.strategy import optimised_strategy

# Some precomputed
# ea, eb, X, Y
# Such that
# A = 2**ea
# B = 3**eb
# C = A - B = X^2 + Y^2
DIAMONDS = [
    (1, 9, 5, 10, 13),
    (1, 72, 41, 28930079307, 62039858138),
    (105, 126, 75, 4153566748924546849, 8198183225858618766),
    (1, 176, 87, 62060540753699060752145450, 303198714682554770479637193),
    (15, 208, 105, 3624033553412861679508963751522, 19956014635540685051671510447227),
    (
        1,
        602,
        363,
        1150828796773097449350776400976394707767314979521270112577233719378561676605426313207017359,
        3908152402269014765841876930143810585616099277004914795408350071193595839964749987840751014,
    ),
]


def generate_splitting_kernel(param_index=0):
    """ """
    f, ea, eb, X, Y = DIAMONDS[param_index]

    # Configure a prime which allows an easy diamond configuration:
    A = ZZ(2**ea)
    B = ZZ(3**eb)
    p = f * 4 * A * B - 1
    F = GF(p**2, name="i", modulus=[1, 0, 1])

    # Cofactor isogeny is the sum of squares
    C = A - B
    assert C == X**2 + Y**2

    # Starting curve
    E0 = EllipticCurve(F, [1, 0])
    iota = E0.automorphisms()[2]

    # Do Bob's SIDH Key-Exchange
    bob_secret = randint(0, B)

    # Compute torsion basis
    P2, Q2 = torsion_basis(E0, 4 * A)
    P2, Q2 = fix_torsion_basis_renes(P2, Q2, ea + 2)
    P3, Q3 = torsion_basis(E0, B)

    # Check automorphism
    assert iota(iota(P2)) == -P2
    assert iota(iota(Q2)) == -Q2

    phiB, _ = isogeny_from_scalar_x_only(E0, B, bob_secret, basis=(P3, Q3))
    phi_P0, phi_Q0 = evaluate_isogeny_x_only(phiB, P2, Q2, 4 * A, B)

    # We pick values such that the aux. is easy to compute
    def aux_endomorphism(P):
        return X * P + Y * iota(P)

    # Kernel which generates the split isogeny
    P1 = aux_endomorphism(P2)
    Q1 = aux_endomorphism(Q2)
    ker_Phi = (P1, Q1, phi_P0, phi_Q0)

    return ker_Phi, (E0, bob_secret)


if __name__ == "__main__":
    # Note: Tests use the slower Mumford isogenies as we know they work from
    # previous code.

    def test_generate_splitting_kernel(param_index=1):
        """ """
        ea, eb, X, Y = DIAMONDS[param_index]

        # Configure a prime which allows an easy diamond configuration:
        A = ZZ(2**ea)
        B = ZZ(3**eb)
        p = 4 * A * B - 1
        F = GF(p**2, name="i", modulus=[1, 0, 1])

        # Cofactor isogeny is the sum of squares
        C = A - B
        assert C == X**2 + Y**2

        # Starting curve
        E0 = EllipticCurve(F, [1, 0])
        iota = E0.automorphisms()[2]

        # Do an SIDH Key-Exchange
        alice_secret = randint(0, A)
        bob_secret = randint(0, B)

        # Compute torsion basis
        P2, Q2 = torsion_basis(E0, 4 * A)
        P2, Q2 = fix_torsion_basis_renes(
            P2, Q2, ea + 2
        )  # Ensure (0,0) will not be in the kernel
        P3, Q3 = torsion_basis(E0, B)

        # Check automorphism
        assert iota(iota(P2)) == -P2
        assert iota(iota(Q2)) == -Q2

        # Compute first isogenies
        phiA, EA = isogeny_from_scalar_x_only(
            E0, A, alice_secret, basis=(4 * P2, 4 * Q2)
        )
        PB, QB = evaluate_isogeny_x_only(phiA, P3, Q3, B, A)

        phiB, EB = isogeny_from_scalar_x_only(E0, B, bob_secret, basis=(P3, Q3))
        PA, QA = evaluate_isogeny_x_only(phiB, P2, Q2, 4 * A, B)

        # Compute the second isogenies
        _, ESA = isogeny_from_scalar_x_only(EB, A, alice_secret, basis=(4 * PA, 4 * QA))
        _, ESB = isogeny_from_scalar_x_only(EA, B, bob_secret, basis=(PB, QB))

        # Make sure the secrets align
        assert ESA.j_invariant() == ESB.j_invariant()

        # We pick values such that the aux. is easy to compute
        def aux_endomorphism(P):
            return X * P + Y * iota(P)

        # Kernel which generates the split isogeny
        P2c = aux_endomorphism(P2)
        Q2c = aux_endomorphism(Q2)

        # Constants needed for the faster 2,2 chain
        N1 = B
        N2 = C
        N_constant = F(N1 + N2) / F(N1 - N2)
        strategy = optimised_strategy(ea - 1, mul_c=0.175)

        # Compute the chain and make sure it splits
        ker_Phi_scaled = (4 * P2c, 4 * Q2c, 4 * PA, 4 * QA)

        Phi, (F1, F2) = compute_richelot_chain(ker_Phi_scaled, ea, N_constant, strategy)

        # Find which of the two codomain curves is our starting curve
        if F1.is_isomorphic(E0):
            E_start = F1
            index = 0
        else:
            E_start = F2
            index = 1

        # Push the torsion basis of EB through the (2,2) isogeny
        PB3, QB3 = torsion_basis(EB, B)

        # One of these points will have full order, which we can
        # recover the secret from
        K_img = evaluate_richelot_chain(Phi, (None, PB3))[index]
        if not has_order_D(K_img, B):
            K_img = evaluate_richelot_chain(Phi, (None, QB3))[index]

        # Isomorphisms back to original curve E0
        isomorphisms = E_start.isomorphisms(E0)
        for iso in isomorphisms:
            K = iso(K_img)

            # Recover secret by solving dlogs
            a, b = BiDLP(K, P3, Q3, B)

            try:
                secret = (Mod(ZZ(b), B) / a).lift()

                # Ensure the collected secret is correct
                _, EB_test = isogeny_from_scalar_x_only(E0, B, secret, basis=(P3, Q3))

                if EB == EB_test:
                    break
            except:
                pass

        print(f"Recovered secret: {secret}")
        print(f"{bob_secret = }")

    for _ in range(10):
        test_generate_splitting_kernel()
