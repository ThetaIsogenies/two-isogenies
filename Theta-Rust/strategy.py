# ==================================================== #
#       Compute optimised strategy for (2,2)-chain     #
# ---------------------------------------------------- #
# Here we deviate from the usual formula by allowing   #
# doubling and images from the elliptic product to     #
# have a distinct cost from the rest of the tree,      #
# helping limit the number of costly images from the   #
# product in favour for slightly more doubling         #
# ==================================================== #

import functools
import sys

# For the long chain, we need a lot of recursion!
sys.setrecursionlimit(1500)

# fmt: off
def optimised_strategy(n, M, S, I):
    """
    A modification of

    Algorithm 60: https://sike.org/files/SIDH-spec.pdf Shown to be appropriate
    for (l,l)-chains in https://ia.cr/2023/508

    Which allows the leftmost branch to have a different cost for the rest of
    the tree. This is partiularly useful for (2,2) isogenies, where the gluing
    doubling and images have a much higher cost than the rest of the tree.

    Thanks to Robin Jadoul for helping with the implementation of this function 
    via personal communication
    """

    
    # Define the costs and initalise the nodes which we store during doubling

    # Normalised Costs
    # left_cost = (8*M + 6*S, 12*M + 12*S)       # (regular_cost, left_branch_cost) Double
    # right_cost = (4*M + 3*S, 82*M + 18*S + I)  # (regular_cost, first_right_cost) Images
    
    # Projective costs
    left_cost = (8*M + 8*S, 12*M + 12*S)       # (regular_cost, left_branch_cost) Double
    right_cost = (4*M + 4*S, 82*M + 18*S + I)  # (regular_cost, first_right_cost) Images
    checkpoints = ({}, {})  # (inner, left edge)

    @functools.cache
    def cost(n, leftmost):
        """
        The minimal cost to get to all children of a height `n` tree.
        If `leftmost` is true, we're still on the leftmost edge of the "outermost" tree

        Updates a global "Check points" which are the points along a branch which we 
        keep for later
        """
        if n <= 1:
            return 0  # no cost here

        c = float("inf")
        for i in range(1, n):  # where to branch off
            # We need `i` moves on the left branch and `n - i` on the right branch
            # to make sure the corresponding subtrees don't overlap and everything
            # is covered exactly once
            thiscost = sum([
                cost(n - i, leftmost),    # We still need to finish off our walk to the left
                i * left_cost[leftmost],  # The cost for the moves on the left branch
                cost(i, False),           # The tree on the right side, now definitely not leftmost
                right_cost[leftmost] + (n - i - 1) * right_cost[False],  # The cost of moving right, maybe one at the first right cost
            ])
            # If a new lower cost has been found, update values
            if thiscost < c:
                c = thiscost
                checkpoints[leftmost][n] = i
        return c

    def convert(n, checkpoints):
        """
        Given a list of checkpoints, convert this to a list of
        the number of doublings to compute and keep before 
        pushing everything through an isogeny. This forces the
        output to match the more usual implementation, e.g.
        https://crypto.stackexchange.com/a/58377

        Warning! Everything about this function is very hacky, but does the job!
        """
        kernels = [n]
        doubles = []
        leftmost = 1

        # We always select the last point in our kernel
        while kernels != []:
            point = kernels[-1]
            if point == 1:
                # Remove this point and push everything through the isogeny
                kernels.pop()
                kernels = [k - 1 for k in kernels]
                leftmost = 0
            else:
                # checkpoints tells us to double this d times
                d = checkpoints[leftmost][point]
                # Remember that we did this
                doubles.append(d)
                kernels.append(point - d)
        return doubles

    # Compute the cost and populate the checkpoints
    c = cost(n, True)

    # Use the checkpoints to compute the list
    l = convert(n, checkpoints)

    return l

# fmt: on
# n = length of chain
# M = Cost of Multiplication
# S = Cost of Squaring
# I = Cost of Inversion
# data = [n, M, S, I]

# NOTE:
# Costings are computed from the msi.rs benchmark
# where the cost is the time in ns for one operation
data_sml = [126, 74, 52, 3314]
data_med = [208, 188, 153, 5939]
data_big = [632, 2717, 2265, 54_823]

strat_sml = optimised_strategy(*data_sml)
print(strat_sml)
strat_med = optimised_strategy(*data_med)
print(strat_med)
strat_big = optimised_strategy(*data_big)
print(strat_big)
