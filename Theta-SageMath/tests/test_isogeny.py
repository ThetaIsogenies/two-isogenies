import time
from sage.all import ZZ, set_random_seed

from utilities.supersingular import torsion_basis
from utilities.strategy import optimised_strategy

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from isogeny_diamond import generate_splitting_kernel, DIAMONDS


def test_isogeny_chain():
    strategy = optimised_strategy(ea)

    # Compute the (2^ea,2^ea)-isogeny
    t0 = time.process_time()
    Phi = EllipticProductIsogeny(ker_Phi, ea, strategy=strategy)
    print(f"Theta Model isogeny took: {time.process_time() - t0:.5f} seconds")

    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(EA(0), PB3)
    t0 = time.process_time()
    _ = Phi(L1)
    print(f"Theta Model image took: {time.process_time() - t0:.5f} seconds")


def test_isogeny_chain_sqrt(ker_Phi):
    strategy = optimised_strategy(ea - 2)
    (P, Q) = ker_Phi
    ker_Phi_scaled = (4 * P, 4 * Q)

    # Compute the (2^ea,2^ea)-isogeny
    t0 = time.process_time()
    Phi = EllipticProductIsogenySqrt(ker_Phi_scaled, ea, strategy=strategy)
    print(f"Theta Model Sqrt isogeny took: {time.process_time() - t0:.5f} seconds")

    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(EA(0), PB3)
    t0 = time.process_time()
    _ = Phi(L1)
    print(f"Theta Model Sqrt image took: {time.process_time() - t0:.5f} seconds")


if __name__ == "__main__":
    test_number = 1
    for test_index in range(0, len(DIAMONDS) - 3):
        set_random_seed(0)
        for _ in range(test_number):
            _, ea, eb, _, _ = DIAMONDS[test_index]
            B = ZZ(3**eb)

            t0 = time.process_time()
            (P1, Q1, P2, Q2), (E0, _) = generate_splitting_kernel(test_index)
            ker_Phi = (CouplePoint(P1, P2), CouplePoint(Q1, Q2))
            print(
                f"Generating splitting kernel with parameters 'ea={ea} eb={eb}' took: {time.process_time() - t0:.5f} seconds"
            )
            EA, EB = ker_Phi[0].curves()

            t0 = time.process_time()
            P3, Q3 = torsion_basis(E0, B)
            PB3, QB3 = torsion_basis(EB, B)
            print(
                f"Generating 3^{eb}-torsion basis on the two elliptic curves took: {time.process_time() - t0:.5f} seconds"
            )
            test_isogeny_chain()
            print("")
            test_isogeny_chain_sqrt(ker_Phi)
            print()
