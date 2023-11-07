# Sage imports
from sage.all import ZZ

# import pari for fast dlog
import cypari2

# ===================================== #
#  Fast DLP solving using Weil pairing  #
# ===================================== #

# Make instance of Pari
pari = cypari2.Pari()


def discrete_log_pari(a, base, order):
    """
    Wrapper around pari discrete log. Works like a.log(b),
    but allows us to use the optional argument order. This
    is important as we skip Pari attempting to factor the
    full order of Fp^k, which is slow.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)


def weil_pairing_pari(P, Q, D, check=False):
    """
    Wrapper around Pari's implementation of the Weil pairing
    Allows the check of whether P,Q are in E[D] to be optional
    """
    if check:
        nP, nQ = D * P, D * Q
        if nP.is_zero() or nQ.is_zero():
            raise ValueError("points must both be n-torsion")

    return pari.ellweilpairing(P.curve(), P, Q, D)


def tate_pairing_pari(P, Q, D):
    """
    Wrapper around Pari's implementation of the Tate pairing
    NOTE: this is not the reduced Tate pairing, so to make this
    match with SageMath you need

    P.tate_pairing(Q, D, k) == pari.elltatepairing(E, P, Q, D)**((p^k - 1) / D)
    """
    E = P.curve()
    return pari.elltatepairing(E, P, Q, D)


def _precompute_baby_steps(base, step, e):
    """
    Helper function to compute the baby steps for
    pohlig_hellman_base and windowed_pohlig_hellman.
    """
    baby_steps = [base]
    for _ in range(e):
        base = base**step
        baby_steps.append(base)
    return baby_steps


def pohlig_hellman_base(a, base, e):
    """
    Solve the discrete log for a = base^x for
    elements base,a of order 2^e using the
    Pohlig-Hellman algorithm.
    """
    baby_steps = _precompute_baby_steps(base, 2, e)

    dlog = 0
    exp = 2 ** (e - 1)

    # Solve the discrete log mod 2, only 2 choices
    # for each digit!
    for i in range(e):
        if a**exp != 1:
            a /= baby_steps[i]
            dlog += 2**i

        if a == 1:
            break

        exp //= 2

    return dlog


def windowed_pohlig_hellman(a, base, e, window):
    """
    Solve the discrete log for a = base^x for
    elements base,a of order 2^e using the
    windowed Pohlig-Hellman algorithm following
    https://ia.cr/2016/963.

    Algorithm runs recursively, computing in windows
    l^wi for window=[w1, w2, w3, ...].

    Runs the base case when window = []
    """
    # Base case when all windows have been used
    if not window:
        return pohlig_hellman_base(a, base, e)

    # Collect the next window
    w, window = window[0], window[1:]
    step = 2**w

    # When the window is not a divisor of e, we compute
    # e mod w and solve for both the lower e - e mod w
    # bits and then the e mod w bits at the end.
    e_div_w, e_rem = divmod(e, w)
    e_prime = e - e_rem

    # First force elements to have order e - e_rem
    a_prime = a ** (2**e_rem)
    base_prime = base ** (2**e_rem)

    # Compute base^(2^w*i) for i in (0, ..., e/w-1)
    baby_steps = _precompute_baby_steps(base_prime, step, e_div_w)

    # Initialise some pieces
    dlog = 0
    if e_prime:
        exp = 2 ** (e_prime - w)

        # Work in subgroup of size 2^w
        s = base_prime ** (exp)

        # Windowed Pohlig-Hellman to recover dlog as
        # alpha = sum l^(i*w) * alpha_i
        for i in range(e_div_w):
            # Solve the dlog in 2^w
            ri = a_prime**exp
            alpha_i = windowed_pohlig_hellman(ri, s, w, window)

            # Update a value and dlog computation
            a_prime /= baby_steps[i] ** (alpha_i)
            dlog += alpha_i * step**i

            if a_prime == 1:
                break

            exp //= step

    # TODO:
    # I don't know if this is a nice way to do
    # this last step... Works well enough but could
    # be improved I imagine...
    exp = 2**e_prime
    if e_rem:
        base_last = base**exp
        a_last = a / base**dlog
        dlog_last = pohlig_hellman_base(a_last, base_last, e_rem)
        dlog += exp * dlog_last

    return dlog


def BiDLP(R, P, Q, D, ePQ=None):
    """
    Given a basis P,Q of E[D] finds
    a,b such that R = [a]P + [b]Q.

    Uses the fact that
        e([a]P + [b]Q, [c]P + [d]Q) = e(P,Q)^(ad-bc)

    Optional: include the pairing e(P,Q) which can be precomputed
    which is helpful when running multiple BiDLP problems with P,Q
    as input. This happens, for example, during compression.
    """
    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = weil_pairing_pari(P, Q, D)

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = weil_pairing_pari(R, Q, D)

    # e(R,-P) = e(P, Q)^b
    pair_b = weil_pairing_pari(R, -P, D)

    # Now solve the dlog in Fq
    a = discrete_log_pari(pair_a, pair_PQ, D)
    b = discrete_log_pari(pair_b, pair_PQ, D)

    return a, b


def BiDLP_power_two(R, P, Q, e, window, ePQ=None):
    r"""
    Same as the above, but uses optimisations using that
    D = 2^e.

    First, rather than use the Weil pairing, we can use
    the Tate pairing which is approx 2x faster.

    Secondly, rather than solve the discrete log naively,
    we use an optimised windowed Pohlig-Hellman.

    NOTE: this could be optimised further, following the
    SIKE key-compression algorithms and for the Tate pairing
    we could reuse Miller loops to save even more time, but
    that seemed a little overkill for a SageMath PoC

    Finally, as the Tate pairing produces elements in \mu_n
    we also have fast inversion from conjugation, but SageMath
    has slow conjugation, so this doesn't help for now.
    """
    p = R.curve().base_ring().characteristic()
    D = 2**e
    exp = (p**2 - 1) // D

    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = tate_pairing_pari(P, Q, D) ** exp

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = tate_pairing_pari(Q, -R, D) ** exp

    # e(R,-P) = e(P, Q)^b
    pair_b = tate_pairing_pari(P, R, D) ** exp

    # Now solve the dlog in Fq
    a = windowed_pohlig_hellman(pair_a, pair_PQ, e, window)
    b = windowed_pohlig_hellman(pair_b, pair_PQ, e, window)

    return a, b


def DLP_power_two(R, P, Q, e, window, ePQ=None, first=True):
    r"""
    This is the same as BiDLP but it only returns either a or b
    depending on whether first is true or false.
    This is used in compression, where we only send 3 of the 4
    scalars from BiDLP
    """
    p = R.curve().base_ring().characteristic()
    D = 2**e
    exp = (p**2 - 1) // D

    # e(P,Q)
    if ePQ:
        pair_PQ = ePQ
    else:
        pair_PQ = tate_pairing_pari(P, Q, D) ** exp

    if first:
        pair_a = tate_pairing_pari(Q, -R, D) ** exp
        x = windowed_pohlig_hellman(pair_a, pair_PQ, e, window)

    else:
        pair_b = tate_pairing_pari(P, R, D) ** exp
        x = windowed_pohlig_hellman(pair_b, pair_PQ, e, window)

    return x
