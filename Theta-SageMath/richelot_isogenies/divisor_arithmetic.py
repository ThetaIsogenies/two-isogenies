"""
The file contains explicit formula for the addition and doubling 
of divisors of Jacobians of hyperelliptic curves with a generic 
sextic model

The work is a generalisation of https://ia.cr/2011/306
    "Group Law Computations on Jacobians of Hyperelliptic Curves"
    by Craig Costello and Kristin Lauter

NOTE: this code is still a work in progress to reduce the total number of 
Fq operations in both addition and doubling

Current costs:

Addition: 25M 4S 1I
Doubling: 32M 6S 1I

NOTE: these formula have been unrolled to avoid the fact that scalar 
multiplication and field element multiplication have the same cost in
SageMath, so things are particularly ugly.
"""


# ================================================ #
#      Helper functions for Jacobian Addition      #
# ================================================ #

def affine_add(f, D1, D2):
    """
    Wrapper function to perform an addition of two divisors
    D1, D2 in Jac(H) where Di = (ui,vi) is expressed in 
    terms of its Mumford coordinates and the hyperelliptic 
    curve is H : y^2 = f(x) and f(x) is a (non-monic) sextic
    polynomial
    """
    f4 = f[4]
    f5 = f[5]
    f6 = f[6]

    u, v = D1
    ud, vd = D2

    u  = u.monic()
    ud = ud.monic()

    u0, u1 = u[0], u[1]
    u0d, u1d = ud[0], ud[1]
    U1 = u1*u1
    U0 = u1*u0

    v0, v1 = v[0], v[1]
    v0d, v1d = vd[0], vd[1]
    U1d = u1d*u1d
    U0d = u1d*u0d

    u1dd, u0dd, v1dd, v0dd, _, _ = _add_generic_compact(u1, u0, v1, v0, U1, U0, u1d, u0d, v1d, v0d, U1d, U0d, f4, f5, f6)

    x = u.parent().gen()
    u_add = x*x + u1dd*x + u0dd
    v_add = x*v1dd + v0dd

    return u_add, v_add


# ================================================ #
#      Helper functions for Jacobian doubling      #
# ================================================ #


def affine_dbl(f, u, v):
    """
    Wrapper function to perform a single doubling
    of a divisor in Jac(H) where D = (u,v) is expressed in 
    terms of its Mumford coordinates and the hyperelliptic 
    curve is H : y^2 = f(x) and f(x) is a (non-monic) sextic
    polynomial
    """
    f2 = f[2]
    f3 = f[3]
    f4 = f[4]
    f5 = f[5]
    f6 = f[6]

    u = u.monic()

    u0, u1 = u[0], u[1]
    v0, v1 = v[0], v[1]

    U1 = u1*u1
    U0 = u1*u0

    u1dd, u0dd, v1dd, v0dd, _, _ = _dbl_divisor_generic(u1, u0, v1, v0, U0, U1, f2, f3, f4, f5, f6)

    x = u.parent().gen()
    u_add = x*x + u1dd*x + u0dd
    v_add = x*v1dd + v0dd

    return u_add, v_add

def affine_dbl_iter(f, u, v, n):
    """
    Wrapper function to perform a chain of n doublings of a 
    divisors D in Jac(H) where D = (u,v) is expressed in 
    terms of its Mumford coordinates and the hyperelliptic 
    curve is H : y^2 = f(x) and f(x) is a (non-monic) sextic
    polynomial
    """
    u = u.monic()
    
    u0, u1 = u[0], u[1]
    v0, v1 = v[0], v[1]

    f2 = f[2]
    f3 = f[3]
    f4 = f[4]
    f5 = f[5]
    f6 = f[6]

    U1 = u1**2
    U0 = u0*u1

    for _ in range(n):
        u1, u0, v1, v0, U1, U0 = _dbl_divisor_generic(u1, u0, v1, v0, U0, U1, f2, f3, f4, f5, f6)

    x = u.parent().gen()
    udd = x*x + u1*x + u0
    vdd = v1*x + v0
    return udd, vdd

# ================================== #
#      Generic Addition Formula      #
# ================================== #

def _add_generic_compact(u1, u0, v1, v0, U1, U0, u1d, u0d, v1d, v0d, U1d, U0d, f4, f5, f6):
    """
    Explicit addition law for D1, D2 in Jac(H) where the Mumford coordinates
    are expressed in terms of their coefficients (u1, u0) and (v1, v0)

    u = x^2 + u1*x + u0
    v = x^2 + v1*x + v0

    The output are the coefficients of the new Mumford coordinates

    The extra coordinates U1 = u1^2, U0 = u1*u0 are included as they appear
    naturally in the multiplications and hence save time for repeated addition
    or doubling
    
    Cost: 25M 4S 1I
    """
    u1S = u1 + u1d
    v0D = v0 - v0d
    v1D = v1 - v1d
    
    M1 = U1 - U1d - u0 + u0d
    M2 = U0d - U0
    M3 = u1  - u1d
    M4 = u0d - u0

    # 4M
    t1 = (M2-v0D)*(v1D-M1)
    t2 = (-v0D-M2)*(v1D+M1)
    t3 = (-v0D+M4)*(v1D-M3)
    t4 = (-v0D-M4)*(v1D+M3)

    # Technically this is 2 * li_num
    l2_num = t1-t2
    l3_num = t3-t4

    # 1M
    # Technically this is 2 * d
    # d = 2*(M4-M2)*(M1+M3) + t3 + t4 - t1 - t2
    d = (M4-M2)*(M1+M3)
    d = d + d
    d = d + t3 + t4 - t1 - t2

    # Montgomery inverse trick
    # 2S 5M
    A = d*d
    B = (l3_num*l3_num - f6*A)
    C = 1/(d*B)
    d_inv = B*C
    d_shifted_inv = d*A*C

    # Compute l2, l3
    # 2M
    l2 = l2_num * d_inv
    l3 = l3_num * d_inv

    # u1''
    # 2M
    u1dd = - u1S - (f5 - l2*l3 - l2*l3) * d_shifted_inv    

    # u0''
    # 6M 1S 
    u0dd = l3*(l3 * (u0 - U1) + l2*u1 + v1) 
    u0dd = u0dd + u0dd
    u0dd = u0dd + l2*l2 - f4
    u0dd = u0dd * d_shifted_inv
    u0dd = u0dd - u1*u1d - u0 - u0d - u1S * u1dd

    # 1S 1M
    U1dd = u1dd*u1dd
    U0dd = u1dd*u0dd

    # v1'', v0''
    # 4M
    v1dd = l3*(u0dd - U1dd + U1 - u0) + l2*(u1dd - u1) - v1
    v0dd = l3*(U0 - U0dd) + l2*(u0dd - u0) - v0

    # Total cost
    # 25M 4S 1I
    return u1dd, u0dd, v1dd, v0dd, U1dd, U0dd

# ================================== #
#      Generic Doubling Formula      #
# ================================== #

def _dbl_divisor_generic(u1, u0, v1, v0, U0, U1, f2, f3, f4, f5, f6):
    """
    Explicit doubling law for D in Jac(H) where the Mumford coordinates
    are expressed in terms of their coefficients (u1, u0) and (v1, v0)

    Expects mumford coordinates

    u = x^2 + u1*x + u0
    v = v1*x + v0

    For an element of the Jacobian with underlying curve model
    H : y^2 = f(x) = f6 x^6 +  ... + f1 x^1 + f0
    Where f(x) is a generic (non-monic) sextic

    The output are the coefficients of the new Mumford coordinates

    The extra coordinates U1 = u1^2, U0 = u1*u0 are included as they appear
    naturally in the multiplications and hence save time for repeated addition
    or doubling

    Cost: 32M 6S 1I
    """
    # Precomputation
    # 2S
    vv  = v1 * v1
    va  = v1 + u1
    va *= va
    va -= vv
    va -= U1

    # Matrix Coeffs
    # 1M
    # 2(v0-va)
    M1  = v0 + v0
    M1 -= va
    M1 -= va
    # 2*v1*(u0+2*U1)
    M2  = U1 + U1
    M2 += u0
    M2 *= v1
    M2 += M2
    # -2*v1
    M3  = -(v1 + v1)
    # va+2*v0
    M4  = v0 + v0
    M4 += va

    # Precomp multiplications
    # 4M
    f6u0 = f6*u0
    f6U1 = f6*U1
    f5u0 = u0*f5
    f5u1 = f5*u1

    # Precomp additions
    f6U1_3  = f6U1 + f6U1 
    f6U1_3 += f6U1            #  3 * f6*U1
    f6U1_4  = f6U1 + f6U1_3   #  4 * f6*U1
    f5u0_2  = f5u0 + f5u0
    f5u1_2  = f5u1 + f5u1     #  2 * f5*u1
    f5u1_3  = f5u1_2 + f5u1   #  3 * f5*u1
    f6u0_3  = f6u0 + f6u0
    f6u0_3 += f6u0            #  3 * f6*u0
    f6u0_6  = f6u0_3 + f6u0_3 #  6 * f6*u0
    f4_2    = f4 + f4         #  2 * f4

    # z1 = U1*(-3*f6U1 + 2*f5u1 - f4) + u0*(2*f5u1 + 3*f6u0 - 2*f4) - vv + f2
    # 2M
    z11  = f5u1_2 - f6U1_3
    z11 -= f4
    z11 *= U1
    z12  = f5u1_2 + f6u0_3
    z12 -= f4_2
    z12 *= u0

    z1  = z11 + z12
    z1 -= vv
    z1 += f2
    
    # z2 = u1*( 6*f6u0 - 4*f6U1 + 3*f5u1 - 2*f4) - 2*f5u0 + f3
    # 1M
    z2  = f6u0_6 - f6U1_4
    z2 += f5u1_3
    z2 -= f4_2
    z2 *= u1
    z2 -= f5u0_2
    z2 += f3

    # Additions
    t11 = M2 - z1
    t12 = z2 - M1
    t21 = -(z1 + M2)
    t22 = z2 + M1
    t31 = -z1 + M4
    t32 =  z2 - M3
    t41 = -(z1 + M4)
    t42 = z2 + M3

    # 4M
    t1 = t11 * t12
    t2 = t21 * t22
    t3 = t31 * t32
    t4 = t41 * t42

    # ell_2 and ell_3 numerators
    l2_num = t1-t2
    l3_num = t3-t4

    # 1M
    # d = t3+t4-t1-t2-2*(M2-M4)*(M1+M3)
    d11 = M4 - M2
    d12 = M1 + M3
    d   = d11 * d12
    d  += d
    d  += t3
    d  += t4
    d  -= t1 
    d  -= t2

    # Montgomery inverse trick
    # 2S 5M
    A = d * d
    # B = l3_num*l3_num - f6*A
    B1 = l3_num * l3_num
    B2 = f6 * A
    B  = B1 - B2
    C  = d * B
    C  = 1 / C

    d_inv = B * C
    d_shifted_inv  = A * C
    d_shifted_inv *= d
    
    # Compute l2, l3
    # 2M
    l2 = l2_num * d_inv
    l3 = l3_num * d_inv

    # u1'' 
    # 2M
    # u1dd = - 2*u1 + (2*l2*l3 - f5) * d_shifted_inv
    u1dd  = l2 * l3
    u1dd += u1dd
    u1dd -= f5
    u1dd *= d_shifted_inv
    u1dd -= u1
    u1dd -= u1

    # 1M precomp
    u1dd_u1 = u1dd * u1

    # u0''
    # 4M 1S 
    # u0dd  = (2*l3*(l3 * (u0 - U1) + l2*u1 + v1) + l2*l2 - f4)
    l2_u1 = l2*u1
    l2_l2 = l2*l2

    u0dd  = u0 - U1
    u0dd *= l3
    u0dd += l2_u1
    u0dd += v1
    u0dd *= l3
    u0dd += u0dd
    u0dd += l2_l2
    u0dd -= f4
    u0dd *= d_shifted_inv

    # u0dd = u0dd - 2*(u1dd*u1 + u0) - U1
    u0dd  -= u0
    u0dd  -= u0
    u0dd  -= u1dd_u1
    u0dd  -= u1dd_u1
    u0dd  -= U1

    # Next precomputations:
    # 1M 1S
    U1dd = u1dd * u1dd
    U0dd = u1dd * u0dd

    # v1'', v0''
    # 4M
    # v1dd = l3*(u0dd - U1dd + U1 - u0) + l2*(u1dd - u1) - v1
    v1dd_1  = u0dd - U1dd
    v1dd_1 += U1
    v1dd_1 -= u0
    v1dd_1 *= l3
    v1dd_2  = u1dd - u1
    v1dd_2 *= l2
    v1dd    = v1dd_1 + v1dd_2
    v1dd   -= v1

    # v0dd = l3*(U0dd - U0) + l2*(u0dd - u0) - v0
    v0dd_1  = U0 - U0dd
    v0dd_1 *= l3
    v0dd_2  = u0dd - u0
    v0dd_2 *= l2
    v0dd  = v0dd_1 + v0dd_2
    v0dd -= v0

    # 32M 6S 1I
    return u1dd, u0dd, v1dd, v0dd, U1dd, U0dd