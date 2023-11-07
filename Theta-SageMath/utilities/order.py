from sage.all import(
    ZZ,
    prod,
    cached_function
)

# ================================================== #
#  Code to check whether a group element has order D #
# ================================================== #

def batch_cofactor_mul_generic(G_list, pis, group_action, lower, upper):
    """
    Input:  A list of elements `G_list`, such that
                G is the first entry and the rest is empty
                in the sublist G_list[lower:upper]
            A list `pis` of primes p such that
                their product is D
            The `group_action` of the group
            Indices lower and upper
    Output: None

    NOTE: G_list is created in place
    """

    # check that indices are valid
    if lower > upper:
        raise ValueError(f"Wrong input to cofactor_multiples()")

    # last recursion step does not need any further splitting
    if upper - lower == 1:
        return

    # Split list in two parts,
    # multiply to get new start points for the two sublists,
    # and call the function recursively for both sublists.
    mid = lower + (upper - lower + 1) // 2
    cl, cu = 1, 1
    for i in range(lower, mid):
        cu = cu * pis[i]
    for i in range(mid, upper):
        cl = cl * pis[i]
    # cl = prod(pis[lower:mid])
    # cu = prod(pis[mid:upper])

    G_list[mid] = group_action(G_list[lower], cu)
    G_list[lower] = group_action(G_list[lower], cl)

    batch_cofactor_mul_generic(G_list, pis, group_action, lower, mid)
    batch_cofactor_mul_generic(G_list, pis, group_action, mid, upper)


@cached_function
def has_order_constants(D):
    """
    Helper function, finds constants to
    help with has_order_D
    """
    D = ZZ(D)
    pis = [p for p, _ in D.factor()]
    D_radical = prod(pis)
    Dtop = D // D_radical
    return Dtop, pis


def has_order_D(G, D, multiplicative=False):
    """
    Given an element G in a group, checks if the
    element has order exactly D. This is much faster
    than determining its order, and is enough for 
    many checks we need when computing the torsion
    basis.

    We allow both additive and multiplicative groups
    which means we can use this when computing the order
    of points and elements in Fp^k when checking the 
    multiplicative order of the Weil pairing output
    """
    # For the case when we work with elements of Fp^k
    if multiplicative:
        group_action = lambda a, k: a**k
        is_identity = lambda a: a == 1
        identity = 1
    # For the case when we work with elements of E / Fp^k
    else:
        group_action = lambda a, k: k * a
        is_identity = lambda a: a.is_zero()
        identity = G.curve()(0)

    if is_identity(G):
        return False

    D_top, pis = has_order_constants(D)

    # If G is the identity after clearing the top
    # factors, we can abort early
    Gtop = group_action(G, D_top)
    if is_identity(Gtop):
        return False

    G_list = [identity for _ in range(len(pis))]
    G_list[0] = Gtop

    # Lastly we have to determine whether removing any prime 
    # factors of the order gives the identity of the group
    if len(pis) > 1:
        batch_cofactor_mul_generic(G_list, pis, group_action, 0, len(pis))
        if not all([not is_identity(G) for G in G_list]):
            return False

    return True