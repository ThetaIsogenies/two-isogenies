from sage.all import proof, cached_method, GF

# ========================== #
#     Speed up SageMath!     #
# ========================== #


def speed_up_sagemath():
    """
    First we set proof.all(False) for general speed ups which
    keeping everything correct (enough)

    Then, we apply a monkey patch to cache the vector_space which
    helps with the performance of polynomial computations
    """
    # Skips strong primality checks and other slow things
    proof.all(False)

    # Cache vector spaces to improve performance
    p = 2**127 - 1  # Arbitrary large prime
    to_patch = [GF(3), GF(3**2), GF(p), GF(p**2)]
    for x in to_patch:
        type(x).vector_space = cached_method(type(x).vector_space)


# ====================== #
#     Print Debugging    #
# ====================== #


def verbose_print(msg, verbose=False):
    if verbose:
        print(msg)


def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)
