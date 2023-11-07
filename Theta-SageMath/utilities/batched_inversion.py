def batched_inversion(*input):
    """
    Given k values, compute the inverse using
    3(k-1) multiplications and one inversion
    """
    # for values a0, a1, ..., an compute
    # a0, a0a1, a0a1a2, ... a0...an using
    # (k-1) multiplications
    multiples = [input[0]]
    for ai in input[1:]:
        multiples.append(ai * multiples[-1])

    # Compute 1 / (a0a1a2...an)
    last_multiple = multiples[-1]
    inverses_multiples = [1 / last_multiple]

    # Compute (a0a1a2...an)^-1, (a0a1a2...a(n-1)) ... a0^-1
    # using k-1 multiplications
    for ai in reversed(input[1:]):
        inverses_multiples.append(ai * inverses_multiples[-1])

    # Reverse for easy ordering below
    # inverses_multiples = a0^-1, (a0a1)^-1 ...(a0a1...an)^-1
    inverses_multiples = inverses_multiples[::-1]
    inverses = [inverses_multiples[0]]

    # Compute the inverse of each element from multiples and
    # their inverses using k-1 multiplications
    k = len(input)
    for i in range(1, k):
        inverses.append(inverses_multiples[i] * multiples[i - 1])

    return inverses
