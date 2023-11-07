"""
The file implements the (2^b,2^b)-isogeny between elliptic products and some helper
functions which FESTA call to compute and evaluate this isogeny.

The (2^b,2^b)-isogeny chain is computed in three steps:

1. A gluing isogeny which takes a pair of supersingular elliptic curves (E, E') 
   and pairs of points (P, Q), (P', Q') and outputs a Jacobian of a hyperelliptic 
   curve J(H) and a pair of divisors (D, D').

2. An isogeny between Jacobians of hyperelliptic curves with domain J(H) and kernel
   (D, D') and outputs the codomain J(H')

3. A splitting isogeny which takes J(H) and outputs a pair of supersingular elliptic
   curves.

The formula used follow https://ia.cr/2022/1283 for the compact explicit formula 
for the gluing isogeny and follow Benjamin Smith's thesis for the Richelot 
correspondence and Splitting isogeny. http://iml.univ-mrs.fr/~kohel/phd/thesis_smith.pdf

This code has been adapted from the file:

  https://github.com/jack4818/Castryck-Decru-SageMath/blob/main/richelot_aux.py

which was originally used for the SageMath implementation of the Castryck-Decru
attack, as well as the variant of the attack described by Oudompheng. See:

  https://ia.cr/2022/975
  https://www.normalesup.org/~oudomphe/textes/202208-castryck-decru-shortcut.pdf
"""

# Sage imports
from sage.all import (
    PolynomialRing,
    EllipticCurve,
    Matrix,
)

# local imports
from richelot_isogenies.divisor_arithmetic import affine_dbl_iter, affine_add
from utilities.supersingular import weil_pairing_pari
from utilities.polynomial_inversion import invert_mod_polynomial_quadratic, invert_mod_polynomial_quartic

def FromProdToJac(P2, Q2, R2, S2):
    """
    Compact explicit formula for the gluing isogeny is derived in
    https://ia.cr/2022/1283
    """
    Fp2 = R2.curve().base()

    # NOTE: we use a generic implementation of the PolynomialRing
    # as this is faster than NTL in this scenario due to the slowness
    # of coefficient extraction and polynomial construction!
    R = PolynomialRing(Fp2, name="x", implementation="generic")
    x = R.gens()[0]

    # Extract roots.
    ai = (P2[0], Q2[0], (P2 + Q2)[0])
    bi = (R2[0], S2[0], (R2 + S2)[0])
    a, a_inv = (a for a in ai if not a.is_zero())
    b, b_inv = (b for b in bi if not b.is_zero())

    # Compute the first two roots of the
    # hyperelliptic curve
    # 2M
    alpha_2 = a*b_inv
    alpha_3 = b*a_inv

    # values to invert
    z1 = a - a_inv
    z2 = b - b_inv
    d  = alpha_2 - alpha_3

    # Montgomery trick for inversions
    # 2I = 3M + 1I
    z1d = z1*d
    z1d_inv = 1/z1d

    d_inv  = z1*z1d_inv
    z1_inv = d*z1d_inv

    # Compute the last root of the 
    # hyperelliptic curve
    # 1M
    alpha_1 = z2 * z1_inv

    # Compute s1, s2 and s1_inv for
    # twisting and eval
    # 3M
    s1 =  z1 * d_inv
    s2 = -z2 * d_inv
    s1_inv = d * z1_inv
    
    # Compute codomain curve equation
    # 2M
    h1, h3, h5 = 0, 0, 0
    h0 = s2
    h2 = s1 - s2*(alpha_2 + alpha_3)
    h4 = -s1*(alpha_1 + alpha_2 + alpha_3)
    h6 = s1

    h_coeffs = [h0,h1,h2,h3,h4,h5,h6]
    h = R(h_coeffs)

    # We need the image of (P1, P2) and (Q1, Q2) in J
    # The image of (P1, P2) is the image of P1 as a divisor on H
    # plus the image of P2 as a divisor on H.
    def isogeny(pair):
        # Argument members may be None to indicate the zero point.
        P1, P2 = pair

        # The projection maps are:
        # H->C: (xC = s1/x²+s2, yC = s1 y)
        # so we compute Mumford coordinates of the divisor f^-1(P_c): a(x), y-b(x)
        if P1:
            xP1, yP1 = P1.xy()
            uP1 = x**2 + (s2 - xP1) * s1_inv
            vP1 = R(yP1 * s1_inv)

            DP1 = (uP1, vP1)

        # Same for E
        # H->E: (xE = s2 x² + s1, yE = s2 y/x^3)
        if P2:
            xP2, yP2 = P2.xy()
            shift_inv = 1 / (xP2 - s1)
            
            uP2 =  x**2 - s2 * shift_inv
        
            # mod uP2 we have that
            # x^2 = s2 / (xP2 - s1)
            # vP2 = yP2 * x**3 * s2_inv
            vP2 = x * yP2 * shift_inv

            DP2 = (uP2, vP2)

        # Now we perform addition of the two divisors
        if P1 and P2:
            return affine_add(h, DP1, DP2)
        if P1:
            return DP1
        if P2:
            return DP2

    return h, isogeny

class RichelotCorr:
    """
    The Richelot correspondence between hyperelliptic
    curves y²=g1*g2*g3 and y²=h1*h2*h3=hnew(x)

    It is defined by equations:
        g1(x1) h1(x2) + g2(x1) h2(x2) = 0
    and y1 y2 = g1(x1) h1(x2) (x1 - x2)

    Given a divisor D in Mumford coordinates:
        U(x) = x^2 + u1 x + u0 = 0
        y = V(x) = v1 x + v0
    Let xa and xb be the symbolic roots of U.
    Let s, p by the sum and product (s=-u1, p=u0)

    Then on x-coordinates, the image of D is defined by equation:
        (g1(xa) h1(x) + g2(xa) h2(x))
      * (g1(xb) h1(x) + g2(xb) h2(x))
    which is a symmetric function of xa and xb.
    This is a non-reduced polynomial of degree 4.

    Write gred = g-U = g1*x + g0
    then gred(xa) gred(xb) = g1^2*p + g1*g0*s + g0^2
    and  g1red(xa) g2red(xb) + g1red(xb) g2red(xa)
       = 2 g11 g21 p + (g11*g20+g10*g21) s + 2 g10*g20

    On y-coordinates, the image of D is defined by equations:
           V(xa) y = Gred1(xa) h1(x) (xa - x)
        OR V(xb) y = Gred1(xb) h1(x) (xb - x)
    If we multiply:
    * y^2 has coefficient V(xa)V(xb)
    * y has coefficient h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
      (x-degree 3)
    * 1 has coefficient Gred1(xa) Gred1(xb) h1(x)^2 (x-xa)(x-xb)
                      = Gred1(xa) Gred1(xb) h1(x)^2 U(x)
      (x-degree 4)
    """
    def __init__(self, G1, G2, H1, H2, hnew):
        self.G1 = G1
        self.G2 = G2
        self.H1 = H1
        self.H11 = H1*H1
        self.H12 = H1*H2
        self.H22 = H2*H2
        self.hnew = hnew
        self.x = hnew.parent().gen()

    def map(self, D):
        "Computes (non-monic) Mumford coordinates for the image of D"
        U, V = D
        U = U.monic()
        V = V % U
        
        # Sum and product of (xa, xb)
        s, p = -U[1], U[0]
        
        # Compute X coordinates (non reduced, degree 4)
        g1red = self.G1 - U
        g2red = self.G2 - U
        g11, g10 = g1red[1], g1red[0]
        g21, g20 = g2red[1], g2red[0]

        # Precompute and reuse some multiplications
        tt = g11*p
        t0 = g11*tt + g11*g10*s + g10*g10
        t1 = g21*tt
        t2 = g10*g20

        # see above
        Px = t0 * self.H11 \
           + (t1 + t1 + (g11*g20 + g21*g10)*s + t2 + t2) * self.H12 \
           + (g21*(g21*p + g20*s) + g20*g20) * self.H22

        # Compute Y coordinates (non reduced, degree 3)
        v1, v0 = V[1], V[0]
        # coefficient of y^2 is V(xa)V(xb)
        Py2 = v1*v1*p + v1*v0*s + v0*v0
        
        # coefficient of y is h1(x) * (V(xa) Gred1(xb) (x-xb) + V(xb) Gred1(xa) (x-xa))
        # so we need to symmetrize:
        # V(xa) Gred1(xb) (x-xb)
        # = (v1*xa+v0)*(g11*xb+g10)*(x-xb)
        # = (v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)*x
        # - xb*(v1*g11*p + v1*g10*xa + v0*g11*xb + v0*g10)
        # Symmetrizing xb^2 gives u1^2-2*u0

        # Precomputing some values we will reuse
        z0 = v1*g11*p
        z1 = v1*g10
        z2 = v0*g11
        z3 = v0*g10

        Py1 = (z0 + z0 + z3 + z3 + s*(z1 + z2))*self.x \
            - (p*(z1 + z1 - z2 - z2) + s*(z0 + s*z2 + z3))
        Py1 *= self.H1
        # coefficient of 1 is Gred1(xa) Gred1(xb) h1(x)^2 U(x)
        Py0 = self.H11 * U * t0

        # Now reduce the divisor, and compute Cantor reduction.
        # Py2 * y^2 + Py1 * y + Py0 = 0
        # y = - (Py2 * hnew + Py0) / Py1
        Py1inv = invert_mod_polynomial_quartic(Py1, Px)
        Py = (- Py1inv * (Py2 * self.hnew + Py0)) % Px

        Dx = ((self.hnew - Py*Py) // Px)
        Dy = (-Py) % Dx
        return (Dx, Dy)

def FromJacToJac(h, D1, D2):
    """
    Isogeny between J(H) -> J(H') with kernel ((D11, D12), (D21, D22))
    where D1, D2 are divisors on J(H) and Dij are the Mumford coordinates
    of the divisors Di
    """
    x = h.parent().gen()

    G1, G2 = D1[0].monic(), D2[0].monic()
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0

    # H1 = 1/det (G2[1]*G3[0] - G2[0]*G3[1])
    #        +2x (G2[2]*G3[0] - G3[2]*G2[0])
    #        +x^2(G2[1]*G3[2] - G3[1]*G2[2])
    # The coefficients correspond to the inverse matrix of delta.

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    delta = delta.inverse()
    H1 = -delta[0][0]*x**2 + 2*delta[1][0]*x - delta[2][0]
    H2 = -delta[0][1]*x**2 + 2*delta[1][1]*x - delta[2][1]
    H3 = -delta[0][2]*x**2 + 2*delta[1][2]*x - delta[2][2]

    # New hyperelliptic curve H'
    hnew = H1*H2*H3

    # Class to compute the evaluation of the isogeny
    R = RichelotCorr(G1, G2, H1, H2, hnew)
    return hnew, R.map

def FromJacToProd(G1, G2, G3, N_constant=None):
    """
    Construct the "split" isogeny from Jac(y^2 = G1*G2*G3)
    to a product of elliptic curves.

    This computation is the same as Benjamin Smith
    see 8.3 in http://iml.univ-mrs.fr/~kohel/phd/thesis_smith.pdf

    Cost: 42M 1S 1I + 1 sqrt
    """
    R = G1.parent()
    x = R.gen()

    # make monic for SL2 transform
    b1, a1, g1 = G1.list()
    b2, a2, g2 = G2.list()
    b3, a3, g3 = G3.list()

    # Montgomery trick to compute
    # inverse of gi 
    # Cost 6M 1I
    g12 = g1*g2
    g123 = g12*g3
    g123_inv = 1/g123
    g12_inv = g3 * g123_inv

    g1_inv  = g2 * g12_inv
    g2_inv  = g1 * g12_inv
    g3_inv  = g12 * g123_inv

    # Hi = x^2 + ai*x + bi
    # Make Hi monic
    # Cost: 6M
    b1, a1 = b1 * g1_inv, a1 * g1_inv
    b2, a2 = b2 * g2_inv, a2 * g2_inv
    b3, a3 = b3 * g3_inv, a3 * g3_inv

    # D^2 = Res(H1, H2)
    r = a1 - a2
    q = b1 - b2

    # Compute resultant of H1, H2
    # Cost: 1M 1S
    if N_constant is None:
        DD = r*(a1*b2 - b1*a2) + q**2
        DD = q**2
        D  = DD.sqrt()
    else:
        D  = q*N_constant
        DD = D * D

    # Mapping to remove linear terms
    u1 = q + D
    u2 = q - D
    v1 = r
    v2 = r
    u_map = R([-u1, u2])
    v_map = R([v1, -v2])

    # Compute coefficients of 
    # Fi = beta_i x^2 + gamma_i
    # Cost: 8M
    X  = DD + DD
    Y1 = D*(q + q - a1*r)
    Y2 = D*(q + q - a2*r)
    Y3 = D*(q + q - a3*r)

    beta1 = g123*(X - Y1)
    beta2 = (X - Y2)
    beta3 = (X - Y3)
    
    gamma1 = g123*(X + Y1)
    gamma2 = (X + Y2)
    gamma3 = (X + Y3)

    # Precompute products of the coefficients
    # to construct the polynomials
    # Cost: 9M
    beta12 = beta1*beta2
    beta23 = beta2*beta3
    beta13 = beta1*beta3
    beta123 = beta1*beta23

    gamma12 = gamma1*gamma2
    gamma23 = gamma2*gamma3
    gamma13 = gamma1*gamma3
    gamma123 = gamma1*gamma23

    betagamma123 = beta123*gamma123

    # Coefficients of the even terms of the transformed
    # sextic
    # Cost: 6M
    c3 = beta123
    c2 = beta23*gamma1 + beta13*gamma2 + beta12*gamma3
    c1 = beta3*gamma12 + beta2*gamma13 + beta1*gamma23
    c0 = gamma123

    # Applying the projection to the above, we get two 
    # cubics
    E1_poly = R([c0,c1,c2,c3])
    E2_poly = E1_poly.reverse()

    # For SageMath, we need this cubic to be monic
    # we work with the scaled coefficients and then
    # map to the curves with `morphE1` and `morphE2`
    # Cost: 4M
    e12 = c2
    e11 = beta123*c1
    e10 = beta123*betagamma123

    e22 = c1
    e21 = gamma123*c2
    e20 = gamma123*betagamma123

    E1 = EllipticCurve([0, e12, 0, e11, e10])
    E2 = EllipticCurve([0, e22, 0, e21, e20])

    def morphE1(x, y):
        # from y^2=p1 to y^2=p1norm
        return (beta123*x, beta123*y)
    
    def morphE2(x, y):
        # from y^2=p2 to y^2=p2norm
        return (gamma123*x, gamma123*y)

    def isogeny(D):
        # To map a divisor, perform the change of coordinates
        # on Mumford coordinates
        U_input, V_input = D
        
        # apply homography
        # y = v1 x + v0 =>
        U = U_input[0] * v_map**2 + U_input[1]*u_map*v_map + U_input[2]*u_map**2
        V = V_input[0] * v_map**3 + V_input[1]*u_map*v_map**2
        V = V % U

        # Extract coefficents from V
        v1, v0 = V[1], V[0]
        
        # Prepare symmetric functions
        s = - U[1] / U[2]
        p = U[0] / U[2]

        # Compute Mumford coordinates on E1
        # Points x1, x2 map to x1^2, x2^2
        U1 = x**2 - (s*s - 2*p)*x + p**2
        # y = v1 x + v0 becomes (y - v0)^2 = v1^2 x^2
        # so 2v0 y-v0^2 = p1 - v1^2 xH^2 = p1 - v1^2 xE1
        V1 = (E1_poly - v1**2 * x + v0**2) / (2*v0)
        # Reduce Mumford coordinates to get a E1 point
        V1 = V1 % U1
        U1red = (E1_poly - V1**2) // U1
        xP1 = -U1red[0] / U1red[1]
        yP1 = V1(xP1)

        # Same for E2
        # Points x1, x2 map to 1/x1^2, 1/x2^2
        U2 = x**2 - (s*s-2*p)/p**2*x + 1/p**2
        # yE = y1/x1^3, xE = 1/x1^2
        # means yE = y1 x1 xE^2
        # (yE - y1 x1 xE^2)(yE - y2 x2 xE^2) = 0
        # p2 - yE (x1 y1 + x2 y2) xE^2 + (x1 y1 x2 y2 xE^4) = 0
        V21 = x**2 * (v1 * (s*s-2*p) + v0*s)
        V20 = E2_poly + x**4 * (p*(v1**2*p + v1*v0*s + v0**2))
        # V21 * y = V20
        V21 = V21 % U2
        V21inv = invert_mod_polynomial_quadratic(V21, U2)
        V2 = (V21inv * V20) % U2

        # Reduce coordinates
        U2red = (E2_poly - V2**2) // U2
        xP2 = -U2red[0] / U2red[1]
        yP2 = V2(xP2)

        return E1(morphE1(xP1, yP1)), E2(morphE2(xP2, yP2))

    return (E1, E2), isogeny

def _check_maximally_isotropic(P, Q, R, S, a):
    """
    Ensures that the kernel ((P, R), (Q, S)) is maximally
    isotropic by ensuring that:

    Isotropic:
    (P, Q) and (R, S) are valid torsion bases for E[2^a] and
    E'[2^a]

    Maximally isotropic:
    e(P, Q) * e(R, S) = 1

    Returns the scaled points 2^(a-1)(P, Q, R, S) on success
    for use in the gluing isogeny
    """
    # Scale to put points in E[2] and E'[2]
    k = 2**(a-1)
    P2 = k*P
    Q2 = k*Q
    R2 = k*R
    S2 = k*S

    two_torsion = (P2, Q2, R2, S2)

    # Ensure both bases are linearly independent
    if P2 == Q2 or R2 == S2:
        return False, [None]*4

    # Ensure all points have order dividing 2^a
    if not all((2*X).is_zero() for X in two_torsion):
        return False, [None]*4
        
    # Ensure all points have order exactly 2^a
    if any(X.is_zero() for X in two_torsion):
        return False, [None]*4

    # Ensure kernel is maximally isotropic
    ePQ = weil_pairing_pari(P, Q, 2*k) # 2k = 2^a
    eRS = weil_pairing_pari(R, S, 2*k)
    if ePQ*eRS != 1:
        return False, [None]*4

    return True, two_torsion

def split_richelot_chain(P, Q, R, S, a, N_constant, strategy):
    r"""
    Given curves C, E and points (P, Q) \in E
                                 (R, S) \in E'
    
    Computes a (2^a, 2^a)-chain of isogenies with kernel:

    ker(Phi) = <(P, R), (Q, S)> : E x E' -> E'' x E''' 

    If the codomain splits, returns the chain and the codomain, 
    and None otherwise

    We expect this to fail only on malformed ciphertexts!
    """
    # We will output the final codomain as well as the isogeny
    # Phi : (E1, E2) -> J -> J -> ... -> J -> J -> (E3, E4)
    # Phi : Glue -> Richelot Chain -> Split
    richelot_chain = []

    # ========================= #
    #  Gluing step              
    #  (E, C) --> Jacobian      #
    # ========================= #

    check, (P2, Q2, R2, S2) = _check_maximally_isotropic(P, Q, R, S, a)
    if not check:
        return None, None

    # Compute the gluing map f : E x E' -> Jac(h)
    h, f = FromProdToJac(P2, Q2, R2, S2)

    # Evaluate the action of f on the points of
    # order 2^b
    D1 = f((P, R))
    D2 = f((Q, S))

    # Store the isogeny for later evaluations
    richelot_chain.append(f)

    # Bookkeeping variables
    strat_idx = 0
    level = [0]

    # We will precompute and push through elements
    # of the kernel, keep track of them with `ker`
    # and `kernel_elements`
    ker = (D1, D2)
    kernel_elements = [ker]

    # ======================================= #
    #  Middle Steps                           #
    #  Jacobian ---> Jacobian                 #
    # ======================================= #
    for i in range(1, a - 1):
        prev = sum(level)
        ker = kernel_elements[-1]
        while prev != (a - 1 - i):
            level.append(strategy[strat_idx])
            # Perform repeated doublings to compute
            # D_new = 2^strategy[strat_idx] D
            D1, D2 = ker
            D1 = affine_dbl_iter(h, D1[0], D1[1], strategy[strat_idx])
            D2 = affine_dbl_iter(h, D2[0], D2[1], strategy[strat_idx])
            ker = (D1, D2)

            # Update kernel elements and bookkeeping variables
            kernel_elements.append(ker)
            prev += strategy[strat_idx]
            strat_idx += 1

        # Compute the next step in the isogeny with the divisors D1, D2
        D1, D2 = ker
        h, f = FromJacToJac(h, D1, D2)
        
        # Update the chain of isogenies
        richelot_chain.append(f)

        # Remove elements from lists
        kernel_elements.pop()
        level.pop()

        # Push the kernel elements through the last step in the isogeny chain
        kernel_elements = [(f(D1), f(D2)) for D1, D2 in kernel_elements]

    # Now we are left with a quadratic splitting: is it singular?
    D1, D2 = kernel_elements[-1]
    G1, G2 = D1[0], D2[0]
    G3, r3 = h.quo_rem(G1 * G2)
    assert r3 == 0

    delta = Matrix(G.padded_list(3) for G in (G1,G2,G3))
    if delta.determinant():
        # Determinant is non-zero, no splitting
        return None, None

    # ======================================= #
    #  Splitting Step                         #
    #  Jacobian ---> (E3, E4)                 #
    # ======================================= #
    h, f = FromJacToProd(G1, G2, G3, N_constant)
    richelot_chain.append(f)
    return richelot_chain, h

def compute_richelot_chain(ker_Phi, b, N_constant, strategy):
    """
    Helper function which takes as input a kernel for
    a (2^b,2^b)-isogeny and returns the isogeny which
    can be used to evaluate points.
    """
    # Unpack kernel generator
    glue_P1, glue_Q1, glue_P2, glue_Q2 = ker_Phi

    chain, domain = split_richelot_chain(
        glue_P1, glue_Q1, glue_P2, glue_Q2, b, N_constant, strategy
    )
    if chain is None:
        raise ValueError("No splitting, ciphertext must be malformed")
    return chain, domain

def evaluate_richelot_chain(Phi, X):
    """
    Evaluate the isogeny Phi : E1 x E2 -> E3 x E4 on a
    pair of points ((P1, P2), (Q1, Q2)) and obtain the
    points ((P3, P4), (Q3, Q4))
    """
    for f in Phi:
        X = f(X)
    return X