# ====================================== #
#  Compute f^(-1) mod g for f,g in R[X]  #
# ====================================== #

def invert_mod_polynomial_quadratic(f, g):
    """
    Given polynomials f, g with deg(g) = 2
    and deg(f) < deg(g) compute f^(-1) mod g

    Inverse is found using linear algebra.
    Cost: 9M 1I
    """
    R = f.parent()

    f0, f1 = f[0], f[1]
    g0, g1, g2 = g[0], g[1], g[2]

    f0_g2 = f0*g2

    A =  f0_g2
    B = -f1*g0
    C =  f1*g2
    D =  f0_g2 - f1*g1

    inv_det = ((A*D - C*B)).inverse()
    inv_det *= g2

    h0 =  D * inv_det
    h1 = -C * inv_det

    return R([h0, h1])

def invert_mod_polynomial_quartic(f, g):
    """
    Given polynomials f, g with deg(g) = 4
    and deg(f) < deg(g) compute f^(-1) mod g

    Inverse is found using linear algebra.
    Cost: 48M 2I

    - 4M 1I (ensure g is monic)
    - 12M   (compute matrix coefficients)
    - 28M   (compute matrix determinant)
    - 4M 1I (solve for inverse)

    - 48M 2I (total)
    """
    R = f.parent()

    # Extract out the coefficients
    f0, f1, f2, f3 = f.list()
    g0, g1, g2, g3, g4 = g.monic().list()

    # First, make g monic
    # 1I 4M
    if not g4.is_one():
        inv_g4 = 1/g4
        g0 *= inv_g4
        g1 *= inv_g4
        g2 *= inv_g4
        g3 *= inv_g4

    # Compute the matrix coefficients
    # M = [[a1, a2, a3, a4]
    #      [b1, b2, b3, b4]
    #      [c1, c2, c3, c4]
    #      [d1, d2, d3, d4]

    # Linear shift and feedback
    # 12M
    a1, b1, c1, d1 = f0, f1, f2, f3

    a2 =    - d1*g0
    b2 = a1 - d1*g1
    c2 = b1 - d1*g2
    d2 = c1 - d1*g3

    a3 =    - d2*g0
    b3 = a2 - d2*g1
    c3 = b2 - d2*g2
    d3 = c2 - d2*g3

    a4 =    - d3*g0
    b4 = a3 - d3*g1
    c4 = b3 - d3*g2
    d4 = c3 - d3*g3


    # Compute the determinant of the matrix
    # 28M
    b1_c2 = b1*c2
    b1_c3 = b1*c3
    b1_c4 = b1*c4

    b2_c1 = b2*c1
    b2_c3 = b2*c3
    b2_c4 = b2*c4

    b3_c1 = b3*c1    
    b3_c2 = b3*c2
    b3_c4 = b3*c4
    
    b4_c1 = b4*c1  
    b4_c2 = b4*c2
    b4_c3 = b4*c3

    D1 = (b3_c4 - b4_c3)*d2 + (b4_c2 - b2_c4)*d3 + (b2_c3 - b3_c2)*d4
    D2 = (b4_c3 - b3_c4)*d1 - (b4_c1 - b1_c4)*d3 + (b3_c1 - b1_c3)*d4
    D3 = (b2_c4 - b4_c2)*d1 + (b4_c1 - b1_c4)*d2 + (b1_c2 - b2_c1)*d4
    D4 = (b3_c2 - b2_c3)*d1 - (b3_c1 - b1_c3)*d2 + (b2_c1 - b1_c2)*d3

    det  = a1*D1 + a2*D2 + a3*D3 + a4*D4

    # Invert the determinant
    # 1I
    det_inv = 1/det

    # Compute solution
    # 4M
    h0 = D1*det_inv
    h1 = D2*det_inv
    h2 = D3*det_inv
    h3 = D4*det_inv

    return R([h0,h1,h2,h3])