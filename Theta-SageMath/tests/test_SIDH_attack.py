from sage.all import ZZ, Mod, gcd

# Python imports
import time

# Local imports
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only
from utilities.order import has_order_D
from utilities.discrete_log import BiDLP
from utilities.supersingular import torsion_basis
from utilities.strategy import optimised_strategy
from utilities.utils import speed_up_sagemath, verbose_print
from isogeny_diamond import generate_splitting_kernel, DIAMONDS


def check_result(E0, EA, EB, B, Phi):
    # Push the torsion basis of EB through the (2,2) isogeny
    P3, Q3 = torsion_basis(E0, B)
    PB3, QB3 = torsion_basis(EB, B)

    E3, E4 = Phi.codomain()

    # Find which of the two codomain curves is our starting curve
    if E3.is_isomorphic(E0):
        E_start = E3
        index = 0
    else:
        assert E4.is_isomorphic(E0)
        E_start = E4
        index = 1

    # One of these points will have full order, which we can
    # recover the secret from
    L1 = CouplePoint(EA(0), PB3)
    L2 = CouplePoint(EA(0), QB3)

    t0 = time.process_time()
    K_img = Phi(L1)[index]
    verbose_print(
        f"Computing an image took: {time.process_time() - t0:.5f}", verbose=verbose
    )

    # Ensure we have a point in the full order
    if not has_order_D(K_img, B):
        K_img = Phi(L2)[index]
    assert has_order_D(K_img, B)

    # Isomorphisms back to original curve E0
    isomorphisms = E_start.isomorphisms(E0)
    for iso in isomorphisms:
        K = iso(K_img)

        # Recover secret by solving dlogs
        a, b = BiDLP(K, P3, Q3, B)

        # fix due to the fact we use E1728 as a starting curve
        if gcd(a, B) != 1:
            iota = E0.automorphisms()[2]
            K = iota(K)
            a, b = BiDLP(K, P3, Q3, B)

        # Recover secret from the BiDLP
        secret = (Mod(ZZ(b), B) / a).lift()

        # Ensure the collected secret is correct
        _, EB_test = isogeny_from_scalar_x_only(E0, B, secret, basis=(P3, Q3))

        if EB == EB_test:
            break

    return secret


def test_SIDH_attack(test_index=-1, verbose=False):
    (P1, Q1, P2, Q2), (E0, bob_secret) = generate_splitting_kernel(test_index)

    # Create kernel from CouplePoint data
    ker_Phi = (CouplePoint(P1, P2), CouplePoint(Q1, Q2))
    EA, EB = ker_Phi[0].curves()

    _, ea, eb, _, _ = DIAMONDS[test_index]
    B = ZZ(3**eb)

    t0 = time.process_time()
    strategy = optimised_strategy(ea)
    Phi = EllipticProductIsogeny(ker_Phi, ea, strategy=strategy)
    verbose_print(f"(2,2)-chain took: {time.process_time() - t0:.5f}s", verbose=verbose)

    ker_Phi_scaled = (4 * ker_Phi[0], 4 * ker_Phi[1])
    strategy_sqrt = optimised_strategy(ea - 2)
    t0 = time.process_time()
    Phi2 = EllipticProductIsogenySqrt(ker_Phi_scaled, ea, strategy=strategy_sqrt)
    verbose_print(
        f"(2,2)-sqrt-chain took: {time.process_time() - t0:.5f}s", verbose=verbose
    )

    verbose_print(f"Original secret: {bob_secret}", verbose=verbose)
    verbose_print("Recovering secret via isogeny chain:", verbose=verbose)
    secret = check_result(E0, EA, EB, B, Phi)
    verbose_print(f"Recovered secret: {secret}", verbose=verbose)
    verbose_print("Recovering secret via sqrt isogeny chain:", verbose=verbose)
    secret2 = check_result(E0, EA, EB, B, Phi2)
    verbose_print(f"Recovered secret via sqrt chain: {secret2}", verbose=verbose)

    assert secret == bob_secret, "Secrets do not match!"
    assert secret2 == bob_secret, "Secrets do not match!"


if __name__ == "__main__":
    speed_up_sagemath()
    test_number = 1
    verbose = True

    for test_index in range(len(DIAMONDS)):
        for _ in range(test_number):
            print(f"Testing index: {test_index}")
            test_SIDH_attack(test_index=test_index, verbose=verbose)
            print("")
