import cProfile
import pstats

from sage.all import ZZ
from isogeny_diamond import generate_splitting_kernel, DIAMONDS
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from utilities.supersingular import torsion_basis
from utilities.strategy import optimised_strategy_old, optimised_strategy
from utilities.utils import speed_up_sagemath

speed_up_sagemath()

test_index = -2
(P1, Q1, P2, Q2), (E0, _) = generate_splitting_kernel(test_index)

# Create kernel from CouplePoint data
ker_Phi = (CouplePoint(P1, P2), CouplePoint(Q1, Q2))
EA, EB = ker_Phi[0].curves()

_, ea, eb, _, _ = DIAMONDS[test_index]

P3, Q3 = torsion_basis(E0, 3**eb)
PB3, QB3 = torsion_basis(EB, 3**eb)

strategy = optimised_strategy_old(ea, mul_c=2)

# start profiler
pr_isogeny = cProfile.Profile()
pr_isogeny.enable()

for _ in range(1000):
    # Compute the (2^ea,2^ea)-isogeny
    Phi = EllipticProductIsogeny(ker_Phi, ea, strategy=strategy)

# End profiler
pr_isogeny.disable()
pr_isogeny.dump_stats("./benchmarks/isogeny.cProfile")
p = pstats.Stats("./benchmarks/isogeny.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(200)


# start profiler
pr_image = cProfile.Profile()
pr_image.enable()

# Compute the time to push a point through the isogeny
L1 = CouplePoint(EA(0), PB3)
for _ in range(1000):
    K_img = Phi(L1)

# End profiler
pr_image.disable()
pr_image.dump_stats("./benchmarks/image.cProfile")
p = pstats.Stats("./benchmarks/image.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(30)
