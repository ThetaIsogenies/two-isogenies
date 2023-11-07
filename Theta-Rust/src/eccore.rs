#![allow(non_snake_case)]

// NOTE:
// The majority of this code was written by Thomas Pornin, which is part of a
// larger work in progress implementing another project which is still ongoing.
// It is similar in form to the macros in the cryptographic library for rust:
// crrl. https://github.com/pornin/crrl

// A macro to create the following Elliptic Curve types:
// - Point:            Type of a point on an elliptic curve (X : Y : Z) with coordinates as elements of Fq
// - PointX:           Type of a point on the Kummer line of an elliptic curve (X : Z) with coordinates as elements of Fq
// - Curve:            Type of an elliptic curve in Montgomery form: E : y^2 = x(x^2 + Ax + 1)
// - CouplePoint       Type of a pair of `Point` P = (P1, P2) in E1 x E2
// - EllipticProduct   Type of a pair of `Curve` representing the product of two curves E1 x E2

// Macro expections:
//    Fq    type of field element
macro_rules! define_ec_core {
    () => {
        use core::ops::Neg;
        use rand_core::{CryptoRng, RngCore};
        use std::fmt;

        /// Curve point.
        /// Points do not know which curve they are on! The caller must ensure
        /// that only proper curve points are used on a given curve; the Rust
        /// type system does not enforce it.
        #[derive(Clone, Copy, Debug)]
        pub struct Point {
            X: Fq,
            Y: Fq,
            Z: Fq,
        }

        impl Point {
            /// The point-at-infinity (neutral element of the group law).
            pub const INFINITY: Self = Self {
                X: Fq::ZERO,
                Y: Fq::ONE,
                Z: Fq::ZERO,
            };

            /// Create a new point: WARNING no check is made
            pub fn new_xy(X: &Fq, Y: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: Fq::ONE,
                }
            }

            /// 0xFFFFFFFF for the point-at-infinity; 0x00000000 otherwise.
            pub fn isinfinity(self) -> u32 {
                self.Z.iszero()
            }

            /// Returns None for the point-at-infinity; otherwise, the (x,y)
            /// affine coordinates.
            pub fn to_xy_vartime(self) -> Option<(Fq, Fq)> {
                if self.Z.iszero() != 0 {
                    return None;
                }
                let t = self.Z.invert();
                Some((self.X * &t, self.Y * &t))
            }

            /// Get the (x,y) affine coordinates. For the point-at-infinity,
            /// this returns (0,0).
            pub fn to_xy(self) -> (Fq, Fq) {
                let t = self.Z.invert();
                (self.X * &t, self.Y * &t)
            }

            /// Returns the X and Z coordinates of the projective point
            pub fn to_xz(self) -> (Fq, Fq) {
                (self.X, self.Z)
            }

            /// Copy rhs into self if ctl == 0xFFFFFFFF.
            /// Do nothing is ctl == 0x00000000.
            /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
            pub fn set_cond(&mut self, rhs: &Self, ctl: u32) {
                self.X.set_cond(&rhs.X, ctl);
                self.Y.set_cond(&rhs.Y, ctl);
                self.Z.set_cond(&rhs.Z, ctl);
            }

            /// Negate this point.
            pub fn set_neg(&mut self) {
                self.Y.set_neg();
            }

            /// Negate this point if ctl == 0xFFFFFFFF.
            /// Do nothing is ctl == 0x00000000.
            /// ctl MUST be either 0xFFFFFFFF or 0x00000000.
            pub fn set_condneg(&mut self, ctl: u32) {
                self.Y.set_condneg(ctl);
            }

            /// Return 0xFFFFFFFF if self and rhs represent the same point.
            /// Otherwise, return 0x00000000.
            pub fn equals(self, rhs: &Self) -> u32 {
                // P1 == P2 if and only if:
                //    P1 == inf AND P2 == inf
                //  OR:
                //    P1 != inf AND P2 != inf AND X1*Z2 = X2*Z1 AND Y1*Z2 = Y2*Z1
                let lz = self.Z.iszero();
                let rz = rhs.Z.iszero();
                let vx = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                let vy = (&self.Y * &rhs.Z).equals(&(&rhs.Y * &self.Z));
                (lz & rz) | (!lz & !rz & vx & vy)
            }
        }

        impl fmt::Display for Point {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(
                    f,
                    "Elliptic Curve Point: ({} : {} : {})",
                    self.X, self.Y, self.Z
                )
            }
        }

        impl Neg for Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = self;
                r.set_neg();
                r
            }
        }

        impl Neg for &Point {
            type Output = Point;

            #[inline(always)]
            fn neg(self) -> Point {
                let mut r = *self;
                r.set_neg();
                r
            }
        }

        /// Special X-only representation of a point (or a pair of points,
        /// since two Y coordinates may match a given X).
        #[derive(Clone, Copy, Debug)]
        pub struct PointX {
            X: Fq,
            Z: Fq,
        }

        impl PointX {
            /// The point-at-infinity (neutral element of the group law).
            #[allow(dead_code)]
            const INFINITY: Self = Self {
                X: Fq::ZERO,
                Z: Fq::ZERO,
            };

            #[allow(dead_code)]
            fn from_point(P: &Point) -> Self {
                Self { X: P.X, Z: P.Z }
            }

            #[allow(dead_code)]
            fn isinfinity(self) -> u32 {
                self.Z.iszero()
            }

            pub fn new_xz(X: &Fq, Z: &Fq) -> Self {
                Self { X: *X, Z: *Z }
            }

            #[allow(dead_code)]
            fn equals(self, rhs: &PointX) -> u32 {
                let inf1 = self.isinfinity();
                let inf2 = rhs.isinfinity();
                let e = (&self.X * &rhs.Z).equals(&(&rhs.X * &self.Z));
                (inf1 & inf2) | (!inf1 & !inf2 & e)
            }
        }

        /// Curve y^2 = x^3 + A*x^2 + x, for a given constant A
        /// (special case of a Montgomery curve).
        #[derive(Clone, Copy, Debug)]
        pub struct Curve {
            A: Fq,   // curve parameter
            A24: Fq, // (A+2)/4
        }

        impl Curve {
            /// Create a new curve instance, with the provided constant.
            pub fn new(A: &Fq) -> Self {
                // We check that the curve is not singular, i.e. A^2 != 4.
                // FIXME: do we want to keep that check?
                assert!(A.equals(&Fq::TWO) == 0);
                assert!((A + &Fq::TWO).iszero() == 0);

                Self {
                    A: *A,
                    A24: (A + &Fq::TWO).half().half(),
                }
            }

            /// Set the point to the provided affine coordinate.
            /// IMPORTANT: this function does NOT check that the point is
            /// really part of the curve.
            pub fn set_xy_nocheck(self, P: &mut Point, x: &Fq, y: &Fq) {
                P.X = *x;
                P.Y = *y;
                P.Z = Fq::ONE;
            }

            /// Set the point P to the provided projective coordinates (use
            /// Z = 0 for the point-at-infinity).
            /// IMPORTANT: this function does NOT check that the point is
            /// really part of the curve.
            pub fn set_xyz_nocheck(self, P: &mut Point, x: &Fq, y: &Fq, z: &Fq) {
                P.X = *x;
                P.Y = *y;
                P.Z = *z;
            }

            /// Set the point P to the provided affine coordinates. This
            /// function checks that the point is on the curve. If the
            /// point is on the curve, then this function returns 0xFFFFFFFF;
            /// otherwise, P is set to the infinity point and this function
            /// returns 0x00000000.
            pub fn set_xy(self, P: &mut Point, x: &Fq, y: &Fq) -> u32 {
                let mut t = x + &self.A;
                t *= &x.square();
                t += x;
                let r = y.square().equals(&t);
                self.set_xy_nocheck(P, x, y);
                P.Z.set_cond(&Fq::ZERO, !r);
                r
            }

            /// Set the point P to the provided projective coordinates (use
            /// Z = 0 for the point-at-infinity). This function checks that
            /// the point is on the curve. If the point is on the curve, then
            /// this function returns 0xFFFFFFFF; otherwise, P is set to the
            /// infinity point and this function returns 0x00000000.
            pub fn set_xyz(self, P: &mut Point, x: &Fq, y: &Fq, z: &Fq) -> u32 {
                let mut t1 = (x + z).square();
                t1 += (&self.A - &Fq::TWO) * &(x * z);
                t1 *= x;
                let t2 = &y.square() * z;
                let r = t1.equals(&t2) | z.iszero();
                self.set_xyz_nocheck(P, x, y, z);
                P.Z.set_cond(&Fq::ZERO, !r);
                r
            }

            /// Get a Point instance set to the provided affine coordinates.
            /// IMPORTANT: this function does NOT check whether the point is on
            /// curve or not.
            pub fn point_xy_nocheck(self, x: &Fq, y: &Fq) -> Point {
                Point {
                    X: *x,
                    Y: *y,
                    Z: Fq::ONE,
                }
            }

            /// Get a Point instance set to the provided projective coordinates
            /// (use Z = 0 for the point-at-infinity).
            /// IMPORTANT: this function does NOT check whether the point is on
            /// curve or not.
            pub fn point_xyz_nocheck(self, x: &Fq, y: &Fq, z: &Fq) -> Point {
                Point {
                    X: *x,
                    Y: *y,
                    Z: *z,
                }
            }

            /// Get a Point instance from the provided affine coordinates.
            /// This returns None if the coordinates do not designate a point
            /// on the curve.
            /// CT: whether the point is on the curve or not may leak
            pub fn point_xy_vartime(self, x: &Fq, y: &Fq) -> Option<Point> {
                let mut P = Point::INFINITY;
                if self.set_xy(&mut P, x, y) == 0 {
                    None
                } else {
                    Some(P)
                }
            }

            /// Get a Point instance from the provided projective coordinates
            /// (use Z = 0 for the point-at-infinity). This returns None if the
            /// coordinates do not designate a point on the curve.
            /// CT: whether the point is on the curve or not may leak
            pub fn point_xyz_vartime(self, x: &Fq, y: &Fq, z: &Fq) -> Option<Point> {
                let mut P = Point::INFINITY;
                if self.set_xyz(&mut P, x, y, z) == 0 {
                    None
                } else {
                    Some(P)
                }
            }

            /// P3 <- P1 + P2
            pub fn add_into(self, P3: &mut Point, P1: &Point, P2: &Point) {
                // Complete routine, to handle all edge cases:
                //   if Z1 == 0:            # P1 == inf
                //       return P2
                //   if Z2 == 0:            # P2 == inf
                //       return P1
                //   L <- Y2*Z1 - Y1*Z2
                //   T <- X2*Z1 - X1*Z2
                //   if T == 0:             # x1 == x2
                //       if L == 0:         # ... and y1 == y2: doubling case
                //           L <- 3*X1^2 + 2*A*X1*Z1 + Z1^2
                //           T <- 2*Y1*Z1
                //       else:              # ... but y1 != y2, thus P2 = -P1
                //           return inf
                //   U <- Z1*Z2*L^2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                //   X3 <- U*T
                //   Y3 <- L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                //   Z3 <- Z1*Z2*T^3
                //
                // Constant-time processing:
                //   Cases P1 == inf and P2 == inf are handled at the end.
                //   (L,T) are always computed for both normal and doubling cases.
                //   If P1 == -P2 then we can let T == 0 and L != 0, this will
                //   properly lead to Z3 == 0.
                //
                // Formulas from https://eprint.iacr.org/2015/1060 are faster
                // but do not cover the case when P1 - P2 is a point of order 2,
                // which can happen in all generality.
                //
                // : current formulas have cost 16M+5S; this can probably
                // be improved. Main issues to tackle:
                //   - Multiplications by A are expensive (since A can be any value)
                //   - There are three points of order 2; this makes finding
                //     complete formulas challenging.

                // T = X2*Z1 - X1*Z2
                // L = Y2*Z1 - Y1*Z2
                let x1z2 = &P1.X * &P2.Z;
                let x2z1 = &P2.X * &P1.Z;
                let mut T = &x2z1 - &x1z2;
                let y1z2 = &P1.Y * &P2.Z;
                let y2z1 = &P2.Y * &P1.Z;
                let mut L = &y2z1 - &y1z2;

                // Alternate (T,L) for doubling:
                //   Td = 2*Y1*Z1
                //   Ld = 3*X1^2 + 2*A*X1*Z1 + Z1^2
                let dbl = T.iszero() & L.iszero();
                let Td = (&P1.Y * &P1.Z).mul2();
                let x1x1 = P1.X.square();
                let z1z1 = P1.Z.square();
                let dx1z1 = &(&P1.X + &P1.Z).square() - &x1x1 - &z1z1;
                let Ld = &x1x1.mul3() + &z1z1 + &self.A * &dx1z1;
                T.set_cond(&Td, dbl);
                L.set_cond(&Ld, dbl);

                // U = L^2*Z1*Z2 - (X1*Z2 + X2*Z1 + A*Z1*Z2)*T^2
                let T2 = T.square();
                let T3 = &T * &T2;
                let z1z2 = &P1.Z * &P2.Z;
                let U = &(&L.square() * &z1z2) - &(&(&x1z2 + &x2z1 + &(&self.A * &z1z2)) * &T2);

                // X3 = U*T
                // Y3 = L*(X1*Z2*T^2 - U) - Y1*Z2*T^3
                // Z3 = Z1*Z2*T^3
                P3.X = &U * &T;
                P3.Y = &(&L * &(&(&x1z2 * &T2) - &U)) - &(&y1z2 * &T3);
                P3.Z = &z1z2 * &T3;

                // Corrective action in case one of the inputs was the
                // point-at-infinity.
                let inf1 = P1.Z.iszero();
                let inf2 = P2.Z.iszero();
                P3.set_cond(&P2, inf1);
                P3.set_cond(&P1, inf2);
            }

            /// P1 <- P1 + P2
            pub fn addto(self, P1: &mut Point, P2: &Point) {
                let mut P3 = Point::INFINITY;
                self.add_into(&mut P3, P1, P2);
                *P1 = P3;
            }

            /// Return P1 + P2 as a new point
            pub fn add(self, P1: &Point, P2: &Point) -> Point {
                let mut P3 = Point::INFINITY;
                self.add_into(&mut P3, P1, P2);
                P3
            }

            /// P3 <- P1 - P2
            pub fn sub_into(self, P3: &mut Point, P1: &Point, P2: &Point) {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.add_into(P3, P1, &nP2);
            }

            /// P1 <- P1 - P2
            pub fn subfrom(self, P1: &mut Point, P2: &Point) {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.addto(P1, &nP2);
            }

            /// Return P1 - P2 as a new point
            pub fn sub(self, P1: &Point, P2: &Point) -> Point {
                let mut nP2 = *P2;
                nP2.set_neg();
                self.add(P1, &nP2)
            }

            #[inline(always)]
            pub fn double_from_coords(self, X: &Fq, Y: &Fq, Z: &Fq) -> (Fq, Fq, Fq) {
                // Doubling formulas in cost 6M+6S
                // These formulas are complete.
                // Formulas from https://eprint.iacr.org/2015/1060 would be
                // more expensive, because multiplications by A are not cheap
                // in the general case.
                //
                // V <- X^2 - Z^2
                // M <- X^2 + Z^2
                // X' <- 2*Y*Z*V^2
                // Y' <- V*(M*(M + 2*A*X*Z) + 4*X^2*Z^2)
                // Z' <- 8*(Y*Z)^3
                let xx = X.square();
                let zz = Z.square();
                let dxz = &(X + Z).square() - &xx - &zz;
                let dyz = (Y * Z).mul2();
                let v = &xx - &zz;
                let m = &xx + &zz;
                let X2 = &dyz * &v.square();
                let Y2 = &v * (&(&m * &(&m + &(&self.A * &dxz))) + &dxz.square());
                let Z2 = &dyz * &dyz.square();

                (X2, Y2, Z2)
            }

            /// P3 <- 2*P1
            pub fn double_into(self, P3: &mut Point, P1: &Point) {
                let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
                P3.X = X2;
                P3.Y = Y2;
                P3.Z = Z2;
            }

            // This is essentially reuse of the above, but allowing the point
            // itself to be mutated... Maybe it would be better to redesign the
            // below functions to stop the code duplication...
            pub fn double_self(self, P1: &mut Point) {
                let (X2, Y2, Z2) = self.double_from_coords(&P1.X, &P1.Y, &P1.Z);
                P1.X = X2;
                P1.Y = Y2;
                P1.Z = Z2;
            }

            /// Return 2*P as a new point
            pub fn double(self, P: &Point) -> Point {
                let mut P3 = Point::INFINITY;
                self.double_into(&mut P3, P);
                P3
            }

            /// Return [2^n]*P as a new point
            pub fn double_iter(self, P: &Point, n: usize) -> Point {
                let mut P3 = *P;
                for _ in 0..n {
                    self.double_self(&mut P3);
                }
                P3
            }

            /// Compute the x-only double of a given point and return
            /// the X-coords
            #[inline(always)]
            pub fn x_dbl_coords(self, X: &Fq, Z: &Fq) -> (Fq, Fq) {
                let mut V1 = (&*X + &*Z).square();
                let V2 = (&*X - &*Z).square();
                let X_new = &V1 * &V2;
                V1 -= &V2;
                let mut Z_new = V1;
                Z_new *= &self.A24;
                Z_new += &V2;
                Z_new *= &V1;

                (X_new, Z_new)
            }

            /// P3 <- n*P
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            pub fn mul_into(self, P3: &mut Point, P: &Point, n: &[u8], nbitlen: usize) {
                // Montgomery ladder: see https://eprint.iacr.org/2017/212

                #[inline(always)]
                fn xdbl(curve: &Curve, X: &mut Fq, Z: &mut Fq) {
                    let mut V1 = (&*X + &*Z).square();
                    let V2 = (&*X - &*Z).square();
                    *X = &V1 * &V2;
                    V1 -= &V2;
                    *Z = V1;
                    *Z *= &curve.A24;
                    *Z += &V2;
                    *Z *= &V1;
                }

                #[inline(always)]
                fn xadd(
                    _curve: &Curve,
                    Xp: &Fq,
                    Zp: &Fq,
                    X0: &Fq,
                    Z0: &Fq,
                    X1: &mut Fq,
                    Z1: &mut Fq,
                ) {
                    let V1 = &(X0 - Z0) * &(&*X1 + &*Z1);
                    let V2 = &(X0 + Z0) * &(&*X1 - &*Z1);
                    *X1 = Zp * &(&V1 + &V2).square();
                    *Z1 = Xp * &(&V1 - &V2).square();
                }

                #[inline(always)]
                fn xadd_aff(_curve: &Curve, Xp: &Fq, X0: &Fq, Z0: &Fq, X1: &mut Fq, Z1: &mut Fq) {
                    let V1 = &(X0 - Z0) * &(&*X1 + &*Z1);
                    let V2 = &(X0 + Z0) * &(&*X1 - &*Z1);
                    *X1 = (&V1 + &V2).square();
                    *Z1 = Xp * &(&V1 - &V2).square();
                }

                // We will need the complete 2*P at the end, to handle some
                // special cases of the formulas.
                let dP = self.double(P);
                let mut X0 = Fq::ONE;
                let mut Z0 = Fq::ZERO;
                let mut X1 = P.X;
                let mut Z1 = P.Z;
                let mut cc = 0u32;
                if nbitlen > 21 {
                    // If n is large enough then it is worthwhile to
                    // normalize the source point to affine.
                    // We do not care if P = inf, since that is handled at
                    // the end in the corrective steps.
                    let Xp = &P.X / &P.Z;
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        xadd_aff(&self, &Xp, &X0, &Z0, &mut X1, &mut Z1);
                        xdbl(&self, &mut X0, &mut Z0);
                        cc = ctl;
                    }
                } else {
                    for i in (0..nbitlen).rev() {
                        let ctl = (((n[i >> 3] >> (i & 7)) as u32) & 1).wrapping_neg();
                        Fq::condswap(&mut X0, &mut X1, ctl ^ cc);
                        Fq::condswap(&mut Z0, &mut Z1, ctl ^ cc);
                        xadd(&self, &P.X, &P.Z, &X0, &Z0, &mut X1, &mut Z1);
                        xdbl(&self, &mut X0, &mut Z0);
                        cc = ctl;
                    }
                }
                Fq::condswap(&mut X0, &mut X1, cc);
                Fq::condswap(&mut Z0, &mut Z1, cc);

                // Special cases:
                //  - ladder fails if P = (0,0) (a point of order 2)
                //  - y is not reconstructed correctly if P has order 2,
                //    or if (n+1)*P = P, -P or infinity.
                let z0z = Z0.iszero();
                let z1z = Z1.iszero();
                let x1ex = (&X1 * &P.Z).equals(&(&P.X * &Z1));

                // (X0/Z0) is the X coordinate of P0 = n*P
                // (X1/Z1) is the X coordinate of P1 = (n + 1)*P
                // We recompute the Y coordinate of n*P (formulas from
                // Okeya and Sakurai).
                let xxzz = &(&P.X * &X0) + &(&P.Z * &Z0);
                let xpz0 = &P.X * &Z0;
                let x0zp = &X0 * &P.Z;
                let zz = &P.Z * &Z0;
                let zzdA = &self.A.mul2() * &zz;
                let u = &(&xxzz * &(&xpz0 + &x0zp + &zzdA)) - &(&zzdA * &zz);
                let v = &P.Y.mul2() * &zz * &Z1;
                P3.X = &X0 * &v;
                P3.Y = &(&u * &Z1) - &(&(&xpz0 - &x0zp).square() * &X1);
                P3.Z = &Z0 * &v;

                // Fix result for the special cases.
                //  P = inf                          -> inf
                //  P != inf, 2*P = inf              -> inf or P (depending on n_0)
                //  2*P != inf, P0 = inf             -> inf
                //  2*P != inf, P0 != inf, P1 = inf  -> -P
                //  2*P != inf, P0 != inf, P1 = -P   -> -2*P
                let order1 = P.Z.iszero();
                let order2 = !order1 & P.Y.iszero();
                let z0inf = !order1 & !order2 & z0z;
                let z1inf = !order1 & !order2 & !z0z & z1z;
                let p1mp = !order1 & !order2 & !z0z & !z1z & x1ex;

                let n_odd = ((n[0] as u32) & 1).wrapping_neg();
                P3.Z.set_cond(&Fq::ZERO, order1 | (order2 & !n_odd) | z0inf);
                P3.set_cond(&P, z1inf | (order2 & n_odd));
                P3.set_cond(&dP, p1mp);
                P3.set_condneg(z1inf | p1mp);
            }

            /// Return n*P as a new point.
            /// Integer n is encoded as unsigned little-endian, with length
            /// nbitlen bits. Bits beyond that length are ignored.
            pub fn mul(self, P: &Point, n: &[u8], nbitlen: usize) -> Point {
                let mut P3 = Point::INFINITY;
                self.mul_into(&mut P3, P, n, nbitlen);
                P3
            }

            /// P3 <- n*P
            /// CT: constant-time for the points, but the integer n may leak.
            pub fn mul_small_into(self, P3: &mut Point, P: &Point, n: u64) {
                match n {
                    0 => {
                        *P3 = Point::INFINITY;
                    }
                    1 => {
                        *P3 = *P;
                    }
                    2 => {
                        self.double_into(P3, P);
                    }
                    3 => {
                        self.add_into(P3, &self.double(P), P);
                    }
                    _ => {
                        self.mul_into(P3, P, &n.to_le_bytes(), (64 - n.leading_zeros()) as usize);
                    }
                }
            }

            /// Return n*P as a new point.
            /// CT: constant-time for the points, but the integer n may leak.
            pub fn mul_small(self, P: &Point, n: u64) -> Point {
                let mut P3 = Point::INFINITY;
                self.mul_small_into(&mut P3, P, n);
                P3
            }

            /// Set P to a random curve point.
            pub fn set_rand_point<T: CryptoRng + RngCore>(self, rng: &mut T, P: &mut Point) {
                // This function cannot actually return the point-at-infinity;
                // this is not a problem as long as the curve order is larger
                // than 2^128.
                P.Z = Fq::ONE;
                loop {
                    P.X.set_rand(rng);
                    P.Y = &(&(&(&P.X + &self.A) * &P.X) * &P.X) + &P.X;
                    if P.Y.legendre() >= 0 {
                        P.Y.set_sqrt();

                        // Randomly chooses the square root to use.
                        let mut tmp = [0u8; 1];
                        rng.fill_bytes(&mut tmp);
                        let ctl = 0u32.wrapping_sub((tmp[0] as u32) & 1);
                        P.Y.set_condneg(ctl);
                        return;
                    }
                }
            }

            /// Return a new random curve point.
            pub fn rand_point<T: CryptoRng + RngCore>(self, rng: &mut T) -> Point {
                let mut P = Point::INFINITY;
                self.set_rand_point(rng, &mut P);
                P
            }

            /// Complete an X-only point into a full point;
            /// (an error is returned if there is no matching Y coordinate).
            /// On error, P3 is set to the point-at-infinity.
            fn complete_pointX_into(self, P3: &mut Point, P: &PointX) -> u32 {
                let XZ = &P.X * &P.Z;
                let V = &(&P.X + &P.Z).square() + &(&(&self.A - &Fq::TWO) * &XZ);
                P3.X = XZ;
                P3.Y = &V * &XZ;
                let ok = P3.Y.set_sqrt();
                P3.Z = P.Z.square();

                // Set to inf on error.
                P3.Z.set_cond(&Fq::ZERO, !ok);

                ok
            }

            /// Complete an X-only point into a full point;
            /// On error, the output point is set to the point-at-infinity.
            pub fn complete_pointX(self, P: &PointX) -> (Point, u32) {
                let mut P3 = Point::INFINITY;
                let ok = self.complete_pointX_into(&mut P3, P);
                (P3, ok)
            }
        }

        impl fmt::Display for Curve {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Montgomery Curve with coefficient: {}", self.A)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct CouplePoint {
            P1: Point,
            P2: Point,
        }

        impl CouplePoint {
            /// Return the pair of points at infinity: O1, O2 on E1 x E2
            pub const INFINITY: Self = Self {
                P1: Point::INFINITY,
                P2: Point::INFINITY,
            };

            /// Create a CouplePoint given a pair of points P1, P2 on E1 x E2
            pub fn new(P1: &Point, P2: &Point) -> Self {
                Self { P1: *P1, P2: *P2 }
            }

            /// Return the points P1, P2
            pub fn points(self) -> (Point, Point) {
                (self.P1, self.P2)
            }
        }

        /// Print debugging, not used within computations.
        impl fmt::Display for CouplePoint {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "Couple point with points:\n{}\n{}", self.P1, self.P2)
            }
        }

        #[derive(Clone, Copy, Debug)]
        pub struct EllipticProduct {
            E1: Curve,
            E2: Curve,
        }

        impl EllipticProduct {
            /// Create an EllipticProduct given a pair of elliptic curves of
            /// type Curve
            pub fn new(E1: &Curve, E2: &Curve) -> Self {
                Self { E1: *E1, E2: *E2 }
            }

            /// Return the pair of curves as a tuple
            pub fn curves(self) -> (Curve, Curve) {
                (self.E1, self.E2)
            }

            /// Addition of elements (P1, P2) and (Q1, Q2) on E1 x E2 is defined
            /// as (P1 + Q1, P2 + Q2). This function calls the add function for
            /// the pair of curves on the EllipticProduct
            pub fn add(self, C1: &CouplePoint, C2: &CouplePoint) -> CouplePoint {
                let mut C3 = CouplePoint::INFINITY;
                C3.P1 = self.E1.add(&C1.P1, &C2.P1);
                C3.P2 = self.E2.add(&C1.P2, &C2.P2);
                C3
            }

            /// Doubles the pair of points (P1, P2) on E1 x E2 as ([2]P1, [2]P2)
            pub fn double(self, C: &CouplePoint) -> CouplePoint {
                let mut C3 = *C;
                C3.P1 = self.E1.double(&C3.P1);
                C3.P2 = self.E2.double(&C3.P2);
                C3
            }

            /// Repeatedly doubles the pair of points (P1, P2) on E1 x E2 to get
            /// ([2^n]P1, [2^n]P2)
            pub fn double_iter(self, C: &CouplePoint, n: usize) -> CouplePoint {
                let mut C3 = *C;
                C3.P1 = self.E1.double_iter(&C3.P1, n);
                C3.P2 = self.E2.double_iter(&C3.P2, n);
                C3
            }
        }
    };
} // End of macro: define_ec_core

pub(crate) use define_ec_core;
