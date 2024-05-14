#![allow(non_snake_case)]

// Macro for defining the following types:
// - ThetaPoint: An element of the level-2 theta structure, which is encoded by
// four projective coordinates: (Fq : Fq : Fq : Fq)
//
// - ThetaStructure: The parent of ThetaPoint, the identity point is the theta
// null point of type ThetaPoint. For arithmetic, this type also has an
// arithmetic precom. which can be reused for both doublings and isogenies.
//
// - product_isogeny: an implementation of an isogeny between elliptic products
// (E1 x E2) of type EllipticProduct = (Curve, Curve) with a kernel of type
// (CouplePoint, CouplePoint) where each CouplePoint represents a pair of points
// P1, P2 on E1 and E2 respectively

// Macro expectations:
// Fq      type of field element Fp^2
// Curve   type of elliptic curve in Montgomery model
// Point   type of point on a montgomery curve
// EllipticProduct    type of E1 x E2
// CouplePoint        type of point on E1 x E2
macro_rules! define_theta_structure {
    () => {
        use std::fmt;

        /// Given four elements of Fq, compute the hadamard transform using recursive
        /// addition.
        /// Cost: 8a
        #[inline(always)]
        fn to_hadamard(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
            let t1 = X + Y;
            let t2 = X - Y;
            let t3 = Z + T;
            let t4 = Z - T;

            let A = &t1 + &t3;
            let B = &t2 + &t4;
            let C = &t1 - &t3;
            let D = &t2 - &t4;
            (A, B, C, D)
        }

        /// Given four elements of Fq, first square each coordinate and
        /// then compute the hadamard transform
        /// Cost: 4S, 8a
        #[inline(always)]
        fn to_squared_theta(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> (Fq, Fq, Fq, Fq) {
            let XX = X.square();
            let YY = Y.square();
            let ZZ = Z.square();
            let TT = T.square();

            to_hadamard(&XX, &YY, &ZZ, &TT)
        }

        // ========================================================
        // Functions for working with ThetaPoints
        // ========================================================

        /// Theta Point Struct
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaPoint {
            X: Fq,
            Y: Fq,
            Z: Fq,
            T: Fq,
        }

        impl ThetaPoint {
            /// Use for initalisation, probably stupid, or at least should have
            /// a different name!
            pub const ZERO: Self = Self {
                X: Fq::ZERO,
                Y: Fq::ZERO,
                Z: Fq::ZERO,
                T: Fq::ZERO,
            };

            /// Compile time, create a new theta point from Fq elements
            pub const fn new(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                    T: *T,
                }
            }

            /// Create a new theta point from Fq elements
            pub fn from_coords(X: &Fq, Y: &Fq, Z: &Fq, T: &Fq) -> Self {
                Self {
                    X: *X,
                    Y: *Y,
                    Z: *Z,
                    T: *T,
                }
            }

            /// Recover the coordinates of the element
            pub fn coords(self) -> (Fq, Fq, Fq, Fq) {
                (self.X, self.Y, self.Z, self.T)
            }

            /// Recover the coordinates of the element
            pub fn list(self) -> [Fq; 4] {
                [self.X, self.Y, self.Z, self.T]
            }

            /// Compute the Hadamard transform of the point's coordinates
            pub fn hadamard(self) -> (Fq, Fq, Fq, Fq) {
                to_hadamard(&self.X, &self.Y, &self.Z, &self.T)
            }

            /// Square each of the point's coordinates and then
            /// compute the hadamard transform
            pub fn squared_theta(self) -> (Fq, Fq, Fq, Fq) {
                to_squared_theta(&self.X, &self.Y, &self.Z, &self.T)
            }
        }

        /// For debugging, pretty print out the coordinates of a point
        impl fmt::Display for ThetaPoint {
            fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "{}\n{}\n{}\n{}\n", self.X, self.Y, self.Z, self.T)
            }
        }

        // ========================================================
        // Functions for working with ThetaStructures
        // ========================================================

        /// Theta Structure
        #[derive(Clone, Copy, Debug)]
        pub struct ThetaStructure {
            null_point: ThetaPoint,
            arithmetic_precom: [Fq; 6],
        }

        impl ThetaStructure {
            /// Given the coordinates of a null point, create a null point and
            /// precompute 6 Fp2 elements which are used for doubling and isogeny
            /// computations.
            pub fn new_from_coords(X: &Fq, Z: &Fq, U: &Fq, V: &Fq) -> Self {
                let null_point = ThetaPoint::new(X, Z, U, V);
                Self {
                    null_point,
                    arithmetic_precom: ThetaStructure::precomputation(&null_point),
                }
            }

            /// Given a null point, store the null point and precompute 6 Fp2
            /// elements which are used for doubling and isogeny computations.
            pub fn new_from_point(null_point: &ThetaPoint) -> Self {
                Self {
                    null_point: *null_point,
                    arithmetic_precom: ThetaStructure::precomputation(null_point),
                }
            }

            /// Return the null point of the ThetaStructure
            pub fn null_point(self) -> ThetaPoint {
                self.null_point
            }

            /// For doubling and also computing isogenies, we need the following
            /// constants, which we can precompute once for each ThetaStructure.
            /// We require 6 Fq elements in total and the cost is
            /// 4S (sqr theta) + 1I + 15M (batch inversion) + 6M (calc)
            /// Cost: 1I + 21M + 4S
            #[inline]
            pub fn precomputation(O0: &ThetaPoint) -> [Fq; 6] {
                let (a, b, c, d) = O0.coords();
                let (AA, BB, CC, DD) = O0.squared_theta();

                // Use Montgomery's trick to match invert k values for a cost
                // of a single inversion and 3(k - 1) multiplications. Inversion
                // is done in place
                let mut inverses = [b, c, d, BB, CC, DD];
                Fq::batch_invert(&mut inverses);

                let y0 = &a * &inverses[0];
                let z0 = &a * &inverses[1];
                let t0 = &a * &inverses[2];

                let Y0 = &AA * &inverses[3];
                let Z0 = &AA * &inverses[4];
                let T0 = &AA * &inverses[5];

                [y0, z0, t0, Y0, Z0, T0]
            }

            /// Given a point P, compute it's double [2]P in place.
            /// Cost 8S + 6M
            #[inline(always)]
            pub fn set_double_self(self, P: &mut ThetaPoint) {
                let (mut xp, mut yp, mut zp, mut tp) = P.squared_theta();

                // Compute temp. coordinates, 8S and 3M
                xp = xp.square();
                yp = &self.arithmetic_precom[3] * &yp.square();
                zp = &self.arithmetic_precom[4] * &zp.square();
                tp = &self.arithmetic_precom[5] * &tp.square();

                // Compute the final coordinates, 3M
                let (X, mut Y, mut Z, mut T) = to_hadamard(&xp, &yp, &zp, &tp);
                Y *= &self.arithmetic_precom[0];
                Z *= &self.arithmetic_precom[1];
                T *= &self.arithmetic_precom[2];

                P.X = X;
                P.Y = Y;
                P.Z = Z;
                P.T = T;
            }

            /// Compute [2] * self
            #[inline]
            pub fn double_point(self, P: &ThetaPoint) -> ThetaPoint {
                let mut P2 = *P;
                self.set_double_self(&mut P2);
                P2
            }

            /// Compute [2^n] * self
            #[inline]
            pub fn double_iter(self, P: &ThetaPoint, n: usize) -> ThetaPoint {
                let mut R = *P;
                for _ in 0..n {
                    self.set_double_self(&mut R)
                }
                R
            }
        }

        // ========================================================
        // Compting the gluing (2,2)-isogeny from a product of
        // elliptic curves to a level 2 theta structure
        //
        // A lot of the code below is to compute a 4x4 matrix,
        // represented as an array [Fq; 16] to compute a symplectic
        // basis transformation to for the points into a form
        // compatible with the isogeny formula
        // ========================================================

        /// Given a point in the four torsion, compute the 2x2 matrix needed
        /// for the basis change
        /// M = [[a, b], [c, d]] represented as an array [a, b, c, d]
        /// Cost: 14M + 2S + 1I
        fn get_base_submatrix(E: &Curve, T: &Point) -> (Fq, Fq, Fq, Fq) {
            let (x, z) = T.to_xz();
            let (u, w) = E.x_dbl_coords(&x, &z); // Cost 3M 2S

            // Precompute some pieces
            let wx = &w * &x;
            let wz = &w * &z;
            let ux = &u * &x;
            let uz = &u * &z;
            let det = &wx - &uz;

            // Batch inversion
            let mut inverse = [det, z];
            Fq::batch_invert(&mut inverse);

            // Compute the matrix coefficients
            let d = &uz * &inverse[0]; // Computing d then a saves one negation
            let a = -&d;
            let b = -&(&wz * &inverse[0]);
            let c = &ux * &inverse[0] - &x * &inverse[1];

            (a, b, c, d)
        }

        /// Given the four torsion below the isogeny kernel, compute the
        /// compatible symplectic transform to allow the theta points to have
        /// a good representation for the gluing isogeny
        ///
        /// Input is expected to be K1 = (P1, P2), K2 = (Q1, Q2) in E1 x E2
        /// Inside (E1 x E2)[4].
        /// Cost 100M + 8S + 4I
        fn get_base_matrix(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
        ) -> [Fq; 16] {
            // First compute the submatrices from each point
            let (E1, E2) = E1E2.curves();
            let (P1, P2) = P1P2.points();
            let (Q1, Q2) = Q1Q2.points();

            // TODO: if these were submatrix computations were done together, we
            // could save 3 inversions... It would make the code harder to read
            // but would be an optimisation for the gluing.
            // Something to think about for when cost REALLY matters.
            // Cost: 4 x 14M + 2S + 1I = 56M + 8S + 4I
            let (g00_1, g01_1, g10_1, g11_1) = get_base_submatrix(&E1, &P1);
            let (g00_2, g01_2, g10_2, g11_2) = get_base_submatrix(&E2, &P2);
            let (h00_1, _, h10_1, _) = get_base_submatrix(&E1, &Q1);
            let (h00_2, h01_2, h10_2, h11_2) = get_base_submatrix(&E2, &Q2);

            // Compute the product of g1 * h1 and g2 * h2 as 2x2 matricies
            // and extract out the first column

            // first col of g1 * h1 = [[gh00_1, *], [gh10_1, *]]
            let gh00_1 = &g00_1 * &h00_1 + &g01_1 * &h10_1;
            let gh10_1 = &g10_1 * &h00_1 + &g11_1 * &h10_1;

            // first col of g2 * h2 = [[gh00_2, *], [gh10_2, *]]
            let gh00_2 = &g00_2 * &h00_2 + &g01_2 * &h10_2;
            let gh10_2 = &g10_2 * &h00_2 + &g11_2 * &h10_2;

            // start the trace with the identity
            let mut a = Fq::ONE;
            let mut b = Fq::ZERO;
            let mut c = Fq::ZERO;
            let mut d = Fq::ZERO;

            // T1
            a += &g00_1 * &g00_2;
            b += &g00_1 * &g10_2;
            c += &g10_1 * &g00_2;
            d += &g10_1 * &g10_2;

            // T2
            a += &h00_1 * &h00_2;
            b += &h00_1 * &h10_2;
            c += &h10_1 * &h00_2;
            d += &h10_1 * &h10_2;

            // T1+T2
            a += &gh00_1 * &gh00_2;
            b += &gh00_1 * &gh10_2;
            c += &gh10_1 * &gh00_2;
            d += &gh10_1 * &gh10_2;

            // Now we act by (0, Q2)
            let a1 = &h00_2 * &a + &h01_2 * &b;
            let b1 = &h10_2 * &a + &h11_2 * &b;
            let c1 = &h00_2 * &c + &h01_2 * &d;
            let d1 = &h10_2 * &c + &h11_2 * &d;

            // Now we act by (P1, 0)
            let a2 = &g00_1 * &a + &g01_1 * &c;
            let b2 = &g00_1 * &b + &g01_1 * &d;
            let c2 = &g10_1 * &a + &g11_1 * &c;
            let d2 = &g10_1 * &b + &g11_1 * &d;

            // Now we act by (P1, Q2)
            let a3 = &g00_1 * &a1 + &g01_1 * &c1;
            let b3 = &g00_1 * &b1 + &g01_1 * &d1;
            let c3 = &g10_1 * &a1 + &g11_1 * &c1;
            let d3 = &g10_1 * &b1 + &g11_1 * &d1;
            // 44M

            [a, b, c, d, a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3]
        }

        /// Apply the base change described by M on a ThetaPoint in-place
        /// Cost: 16M
        #[inline]
        fn apply_base_change(P: &mut ThetaPoint, M: [Fq; 16]) {
            let (x, y, z, t) = P.coords();
            P.X = &M[0] * &x + &M[1] * &y + &M[2] * &z + &M[3] * &t;
            P.Y = &M[4] * &x + &M[5] * &y + &M[6] * &z + &M[7] * &t;
            P.Z = &M[8] * &x + &M[9] * &y + &M[10] * &z + &M[11] * &t;
            P.T = &M[12] * &x + &M[13] * &y + &M[14] * &z + &M[15] * &t;
        }

        /// Given a couple point as input, compute the corresponding ThetaPoint on
        /// the level two structure and then apply the basis change on this point
        /// Cost: 20M
        fn base_change_couple_point(P1P2: &CouplePoint, M: [Fq; 16]) -> ThetaPoint {
            let (P1, P2) = P1P2.points();
            let (mut X1, mut Z1) = P1.to_xz();
            let (mut X2, mut Z2) = P2.to_xz();

            // If we have the point (0, 0) swap to (1, 0)
            let P1_check = X1.iszero() & Z1.iszero();
            X1.set_cond(&Fq::ONE, P1_check);
            Z1.set_cond(&Fq::ZERO, P1_check);

            // If we have the point (0, 0) swap to (1, 0)
            let P2_check = X2.iszero() & Z2.iszero();
            X2.set_cond(&Fq::ONE, P2_check);
            Z2.set_cond(&Fq::ZERO, P2_check);

            // Take all products to get level-2 theta point
            let X = &X1 * &X2;
            let Y = &X1 * &Z2;
            let Z = &Z1 * &X2;
            let T = &Z1 * &Z2;
            let mut P = ThetaPoint::from_coords(&X, &Y, &Z, &T);

            // Finally apply the base change on the point
            apply_base_change(&mut P, M);
            P
        }

        /// For a theta point which is on an elliptic product,
        /// one of the dual coordinates will be zero. We need
        /// to identify the index of this zero element for the
        /// gluing isogeny codomain and evaluation functions
        fn zero_index(dual_coords: &[Fq; 4]) -> usize {
            let mut z_idx = 0;
            for (i, el) in dual_coords.iter().enumerate() {
                // When el is zero, the result is 0xFF...FF
                // and zero otherwise, so we can use this as
                // a mask for each step.
                let el_is_zero = el.iszero();
                z_idx |= (i as u32 & el_is_zero);
            }
            z_idx as usize
        }

        /// Given the 8-torsion above the kernel of order 2, computes the
        /// codomain ThetaStructure (2,2)-isogenous from a product of elliptic
        /// curves
        ///
        /// NOTE: this function is a little fussy as we need to avoid the
        /// zero-dual coordinate. There's a chance refactoring this could make
        /// it appear more friendly
        ///
        /// Cost: 8S 13M 1I
        fn gluing_codomain(
            T1_8: &ThetaPoint,
            T2_8: &ThetaPoint,
        ) -> (ThetaStructure, (Fq, Fq), usize) {
            // First construct the dual coordinates of the kernel and look
            // for the element which is zero
            // For convenience we pack this as an array instead of a tuple:
            let xAxByCyD: [Fq; 4] = T1_8.squared_theta().into();
            let zAtBzYtD: [Fq; 4] = T2_8.squared_theta().into();

            // One element for each array above will be zero. Identify this
            // element to get the right permutation below for filling arrays
            let z_idx = zero_index(&xAxByCyD);

            // Compute intermediate values for codomain
            let t1 = zAtBzYtD[1 ^ z_idx];
            let t2 = xAxByCyD[2 ^ z_idx];
            let t3 = zAtBzYtD[3 ^ z_idx];
            let t4 = xAxByCyD[3 ^ z_idx];

            // Invert all four values for codomain and images
            let mut inverse = [t1, t2, t3, t4];
            Fq::batch_invert(&mut inverse);

            // Codomain coefficients
            let mut ABCD = [Fq::ZERO; 4];
            ABCD[0 ^ z_idx] = Fq::ZERO;
            ABCD[1 ^ z_idx] = &t1 * &inverse[2];
            ABCD[2 ^ z_idx] = &t2 * &inverse[3];
            ABCD[3 ^ z_idx] = Fq::ONE;

            // Used for the image computation
            let a_inverse = &t3 * &inverse[0];
            let b_inverse = &t4 * &inverse[1];

            let (A, B, C, D) = to_hadamard(&ABCD[0], &ABCD[1], &ABCD[2], &ABCD[3]);
            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            (codomain, (a_inverse, b_inverse), z_idx)
        }

        /// Given a point and it's shifted value, compute the image of this
        /// point using the Fq elements precomputed during the codomain
        /// computation
        /// Note: The shift value is needed as one element of the dual coordinates
        /// is zero, so we take values from both and use linear algebra to recover
        /// the correct image.
        ///
        /// Cost: 8S + 10M + 1I
        fn gluing_image(
            T: &ThetaPoint,
            T_shift: &ThetaPoint,
            a_inv: &Fq,
            b_inv: &Fq,
            z_idx: usize,
        ) -> ThetaPoint {
            // Find dual coordinates of point to push through
            let AxByCzDt: [Fq; 4] = T.squared_theta().into();

            // We are in the case where at most one of A, B, C, D is
            // zero, so we need to account for this
            // To recover values, we use the translated point to get
            let AyBxCtDz: [Fq; 4] = T_shift.squared_theta().into();

            // We can always directly compute three elements
            let y = AxByCzDt[1 ^ z_idx] * a_inv;
            let z = AxByCzDt[2 ^ z_idx] * b_inv;
            let t = AxByCzDt[3 ^ z_idx];

            // To compute the `x` value, we need to compute a scalar, lambda,
            // which we use to normalise the given `xb`. We can always do this,
            // but we are required to compute an inverse which means checking
            // whether z is zero. If z is zero, we can compute lambda by simply
            // extracting the scaled zb and dividing by b from above. However,
            // when z is zero, we instead have to use t. To ensure that this is
            // constant time, we compute both using that inverting zero just
            // gives zero and conditionally swapping lanbda with lambda_t
            let zb = AyBxCtDz[3 ^ z_idx];
            let tb = &AyBxCtDz[2 ^ z_idx] * b_inv;

            let mut inverse = [zb, tb];
            Fq::batch_invert(&mut inverse);

            // Potentially one of these inverses are zero, but we do both
            // to avoid branching.
            let mut lam = &z * &inverse[0];
            let lam_t = &t * &inverse[1];
            lam.set_cond(&lam_t, z.iszero());

            // Finally we recover x
            let xb = AyBxCtDz[1 ^ z_idx] * a_inv;
            let x = xb * lam;

            // We now have values for `x,y,z,t` but to order them we need to use
            // the xor trick as above, so we pack them into an array with the
            // right ordering and then extract them back out
            let mut xyzt = [Fq::ZERO; 4];
            xyzt[0 ^ z_idx] = x;
            xyzt[1 ^ z_idx] = y;
            xyzt[2 ^ z_idx] = z;
            xyzt[3 ^ z_idx] = t;

            let (x, y, z, t) = to_hadamard(&xyzt[0], &xyzt[1], &xyzt[2], &xyzt[3]);

            ThetaPoint::from_coords(&x, &y, &z, &t)
        }

        /// Compute the gluing (2,2)-isogeny from a ThetaStructure computed
        /// from an elliptic product.
        fn gluing_isogeny(
            E1E2: &EllipticProduct,
            P1P2_8: &CouplePoint,
            Q1Q2_8: &CouplePoint,
            image_points: &[CouplePoint],
        ) -> (ThetaStructure, Vec<ThetaPoint>) {
            // First recover the four torsion below the 8 torsion
            let P1P2_4 = E1E2.double(&P1P2_8);
            let Q1Q2_4 = E1E2.double(&Q1Q2_8);

            // Use the four torsion to deterministically find basis change
            let M = get_base_matrix(&E1E2, &P1P2_4, &Q1Q2_4);

            // Take the points P1, P2 in E1 x E2 and represent them
            // as a theta point on a level 2 structure and map them
            // through the above basis change
            let T1_8 = base_change_couple_point(&P1P2_8, M);
            let T2_8 = base_change_couple_point(&Q1Q2_8, M);

            // Now it's time to compute the codomain and image of the isogeny
            // with kernel below T1, and T2.
            // We save the zero index, as we use it for the images, and we also
            // can precompute a few inverses to save time for evaluation.
            let (codomain, (a_inv, b_inv), z_idx) = gluing_codomain(&T1_8, &T2_8);

            // We now want to push through a set of points by evaluating each of them
            // under the action of this isogeny. As the domain is an elliptic product,
            // with elements of type CouplePoint, and the codomain is a ThetaStructure
            // with elements of type ThetaPoint, we need a new vector here and as we
            // iteratate through each CouplePoint, we can compute its image and push it
            // to the new vector.
            let mut theta_images: Vec<ThetaPoint> = Vec::new();

            // Per image cost =
            // 2 * (16M + 5S) for the CouplePoint addition
            // 2 * 20M for the base change
            // 8S + 4M + 1I for the gluing image
            // Total:
            // 76M + 18S + 1I per point
            for P in image_points.iter() {
                // Need affine coordinates here to do an add, if we didn't we
                // could use faster x-only... Something to think about but no
                // obvious solution.
                let P_sum_T = E1E2.add(P, &P1P2_4);

                // After we have added the points, we can use the gluing formula
                // to recover theta points on the level 2 theta structure. First we
                // must compute the basis change as we did for the kernel:
                let T = base_change_couple_point(&P, M);
                let T_shift = base_change_couple_point(&P_sum_T, M);

                // With a point and the shift value from the kernel, we can find
                // the image
                let T_image = gluing_image(&T, &T_shift, &a_inv, &b_inv, z_idx);
                theta_images.push(T_image);
            }

            (codomain, theta_images)
        }

        // ===================================================================
        // Compting general (2,2)-isogenies between theta structures
        //
        // NOTE: For the two steps before a product structure is reached, we
        // need additional symplectic transforms which is controlled by the
        // `hadamard` array of `bool`s. The purpose of these is to avoid null
        // points (or dual null points) which have zero elements, which are
        // incompatible with the doubling formula.
        // ===================================================================

        /// Given the 8-torsion above the kernel, compute the codomain of the
        /// (2,2)-isogeny and the image of all points in `image_points`
        /// Cost:
        /// Codomain: 8S + 13M + 1I
        /// Image: 4S + 3M
        fn two_isogeny(
            domain: &ThetaStructure,
            T1: &ThetaPoint,
            T2: &ThetaPoint,
            image_points: &mut [ThetaPoint],
            hadamard: [bool; 2],
        ) -> ThetaStructure {
            // Compute the squared theta transform of both elements
            // of the kernel
            let (xA, xB, _, _) = T1.squared_theta();
            let (zA, tB, zC, tD) = T2.squared_theta();

            // Batch invert denominators
            let mut inverse = [xA, zA, tB];
            Fq::batch_invert(&mut inverse);

            // Compute the codomain coordinates
            let mut A = Fq::ONE;
            let mut B = &xB * &inverse[0];
            let mut C = &zC * &inverse[1];
            let mut D = &tD * &inverse[2] * &B;

            // Inverses will be used for evaluation below
            let B_inv = &domain.arithmetic_precom[3] * &B;
            let C_inv = &domain.arithmetic_precom[4] * &C;
            let D_inv = &domain.arithmetic_precom[5] * &D;

            // Finish computing the codomain coordinates
            // For the penultimate case, we skip the hadamard transformation
            if hadamard[1] {
                (A, B, C, D) = to_hadamard(&A, &B, &C, &D);
            }
            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            // Now push through each point through the isogeny
            for P in image_points.iter_mut() {
                let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();
                if hadamard[0] {
                    (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                    (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
                } else {
                    (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);
                }
                YY *= &B_inv;
                ZZ *= &C_inv;
                TT *= &D_inv;

                if hadamard[1] {
                    (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                }

                P.X = XX;
                P.Y = YY;
                P.Z = ZZ;
                P.T = TT;
            }

            codomain
        }

        /// Special function for the case when we are (2,2)-isogenous to a
        /// product of elliptic curves. Essentially the same as above, but with
        /// some small changes to deal with that a dual coordinate is now zero.
        /// Computes the codomain of the (2,2)-isogeny and the image of all
        /// points in `image_points`
        /// Cost:
        /// Codomain: 8S + 23M + 1I
        /// Image: 4S + 3M
        fn two_isogeny_to_product(
            T1: &ThetaPoint,
            T2: &ThetaPoint,
            image_points: &mut [ThetaPoint],
        ) -> ThetaStructure {
            // Compute the squared theta transform of both elements
            // of the kernel
            let (mut xA, mut xB, yC, yD) = T1.hadamard();
            (xA, xB, _, _) = to_squared_theta(&xA, &xB, &yC, &yD);

            let (mut zA, mut tB, mut zC, mut tD) = T2.hadamard();
            (zA, tB, zC, tD) = to_squared_theta(&zA, &tB, &zC, &tD);

            // Batch invert denominators
            let mut inverse = [xA, zA, tB, xB, zC, tD];
            Fq::batch_invert(&mut inverse);

            // Compute the codomain coordinates as well as precomputations for
            // the images
            let A = Fq::ONE;
            let B = &xB * &inverse[0];
            let C = &zC * &inverse[1];
            let D = &tD * &inverse[2] * &B;
            let B_inv = &xA * &inverse[3];
            let C_inv = &zA * &inverse[4];
            let D_inv = &tB * &inverse[5] * &B_inv;

            let codomain = ThetaStructure::new_from_coords(&A, &B, &C, &D);

            for P in image_points.iter_mut() {
                let (mut XX, mut YY, mut ZZ, mut TT) = P.coords();

                (XX, YY, ZZ, TT) = to_hadamard(&XX, &YY, &ZZ, &TT);
                (XX, YY, ZZ, TT) = to_squared_theta(&XX, &YY, &ZZ, &TT);

                YY *= B_inv;
                ZZ *= C_inv;
                TT *= D_inv;

                P.X = XX;
                P.Y = YY;
                P.Z = ZZ;
                P.T = TT;
            }

            codomain
        }

        // ========================================================
        // Compting the symplectic transform to expose the
        // product structure and then compute the correct
        // splitting to Montgomery curves.
        // ========================================================

        // This function is a bit of a mess. Ultimately, we want to know whether
        // given some pair of indices whether we should multiply by minus one.
        // We do this by returning either 0: do nothing, or 0xFF...FF: negate
        // the value, which concretely is performed with set_negcond() on the
        // field element.
        //
        // Mathematically we have a few things to juggle. Firstly, although the
        // index should really be tuples (x, y) for x,y in {0,1} we simply index
        // from {0, ..., 3}. So there is first the identification of:
        //
        // 0 : (0, 0)
        // 1 : (1, 0)
        // 2 : (0, 1)
        // 3 : (1, 1)
        //
        // The next thing we need is the dot product of these indices
        // For example:
        // Take i . j is the dot product, so input (x, y) = (1, 3)
        // corresponds to computing:
        // (1, 0) . (1, 1) = 1*1 + 0*1 = 1
        //
        // This means evaluation of chi means the sign is dictated by
        // => (-1)^(i.j) = (-1)^1 = -1
        //
        // A similar thing is done for all pairs of indices below.
        //
        // TODO: there may be a nicer way to organise this function, but
        // I couldn't find a nice closed form for ±1 from a pair (i, j)
        // which i could compute on the fly without first matching from
        // x,y in {0,..,3} to i,j in {(0,0)...(1,1)} (which would mean
        // using a match anyway!!).
        fn chi_eval(x: &usize, y: &usize) -> u32 {
            match (x, y) {
                (0, 0) => 0,
                (0, 1) => 0,
                (0, 2) => 0,
                (0, 3) => 0,
                (1, 0) => 0,
                (1, 1) => u32::MAX,
                (1, 2) => 0,
                (1, 3) => u32::MAX,
                (2, 0) => 0,
                (2, 1) => 0,
                (2, 2) => u32::MAX,
                (2, 3) => u32::MAX,
                (3, 0) => 0,
                (3, 1) => u32::MAX,
                (3, 2) => u32::MAX,
                (3, 3) => 0,
                _ => 1,
            }
        }

        /// For a given index (chi, i) compute the level 2,2 constants (square).
        /// The purpose of this is to identify for which (chi, i) this constant
        /// is zero.
        fn level_22_constants_sqr(null_point: &ThetaPoint, chi: &usize, i: &usize) -> Fq {
            let mut U_constant = Fq::ZERO;
            let null_coords = null_point.list();

            for t in 0..4 {
                let mut U_it = &null_coords[t] * &null_coords[i ^ t];
                U_it.set_condneg(chi_eval(chi, &t));
                U_constant += &U_it;
            }
            U_constant
        }

        /// For each possible even index compute the level 2,2 constant. Return
        /// the even index for which this constant is zero. This only fails for
        /// bad input in which case the whole chain would fail. Evaluates all
        /// positions, and so should run in constant time.
        fn identify_even_index(null_point: &ThetaPoint) -> (usize, usize) {
            const EVEN_INDICIES: [(usize, usize); 10] = [
                (0, 0),
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
                (3, 0),
                (3, 3),
            ];
            // Initialise the return tuple
            let mut chi_zero = 0;
            let mut i_zero = 0;

            for (chi, i) in EVEN_INDICIES.iter() {
                let U_sqr = level_22_constants_sqr(null_point, chi, i);

                // When U_sqr is zero, U_sqr_is_zero = 0xFF...FF
                // and 0 otherwise, so we can use this as a mask
                // to select the non-zero index through the loop
                let U_sqr_is_zero = U_sqr.iszero();
                chi_zero |= (*chi as u32 & U_sqr_is_zero);
                i_zero |= (*i as u32 & U_sqr_is_zero);
            }
            (chi_zero as usize, i_zero as usize)
        }

        /// We can precompute 10 different symplectic transforms which
        /// correspond to each of the possible 10 even indicies which could be
        /// zero. We can select the right change of basis by using the above
        /// functions and then selecting the correct map accordingly.
        fn compute_splitting_matrix(null_point: &ThetaPoint) -> [Fq; 16] {
    #[rustfmt::skip]
            const MAPS: [[Fq; 16]; 10] = [
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::MINUS_ONE, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                    Fq::ZERO, Fq::MINUS_ONE, Fq::ZERO, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ZERO, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ONE, Fq::ZERO, Fq::ZERO,
                    Fq::ZERO, Fq::ZERO, Fq::ZERO, Fq::ONE,
                    Fq::ZERO, Fq::ZERO, Fq::ONE, Fq::ZERO,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::ONE, Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::ONE, Fq::MINUS_ONE,
                ],
                [
                    Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE, Fq::ONE,
                    Fq::MINUS_ONE, Fq::MINUS_ONE, Fq::ONE,
                    Fq::ONE, Fq::ONE, Fq::MINUS_ONE, Fq::MINUS_ONE,
                    Fq::MINUS_ONE, Fq::ONE, Fq::MINUS_ONE, Fq::ONE,
                ],
                [
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::ZETA,
                    Fq::ONE, Fq::MINUS_ZETA, Fq::ONE, Fq::ZETA,
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::MINUS_ZETA,
                    Fq::ONE, Fq::ZETA, Fq::ONE, Fq::ZETA,
                ],
            ];

            // Identity the current location of the zero
            let zero_location = identify_even_index(null_point);

            // Compute the corresponding matrix to map the zero to
            // the desired place
            // TODO: is a match like this the best thing to do in Rust??
            let M: [Fq; 16];
            match zero_location {
                (0, 2) => M = MAPS[0],
                (3, 3) => M = MAPS[1],
                (0, 3) => M = MAPS[2],
                (2, 1) => M = MAPS[3],
                (0, 1) => M = MAPS[4],
                (1, 2) => M = MAPS[5],
                (2, 0) => M = MAPS[6],
                (3, 0) => M = MAPS[7],
                (1, 0) => M = MAPS[8],
                (0, 0) => M = MAPS[9],
                // The above locations are an exhaustive list of possible inputs, not sure how to tell rust this...
                _ => panic!("Unreachable"),
            }

            M
        }

        /// Map from a theta point to one which admits a splitting to elliptic
        /// products. Essentially requires computing the correct splitting
        /// matrix and then applying the isomorphism
        fn splitting_isomorphism(
            Th: ThetaStructure,
            image_points: &mut [ThetaPoint],
        ) -> ThetaStructure {
            // Compute the correct splitting matrix
            let mut O0 = Th.null_point();
            let M = compute_splitting_matrix(&O0);

            // Map the Theta Structure through the symplectic transform
            apply_base_change(&mut O0, M);

            // Map the points through the symplectic transform
            for P in image_points.iter_mut() {
                apply_base_change(P, M);
            }

            ThetaStructure::new_from_point(&mut O0)
        }

        /// Given a Theta point in the correct representation, compute two
        /// dimension 1 theta points.
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn split_theta_point(P: &ThetaPoint) -> ((Fq, Fq), (Fq, Fq)) {
            let (a, b, _, d) = P.coords();

            let P1 = (a, b);
            let P2 = (b, d);

            (P1, P2)
        }

        /// Given a dimension one null theta point, compute the corresponding
        /// elliptic curve in the Montgomery model by recovering the Montgomery
        /// coefficient A
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn null_point_to_montgomery_curve(O0: &(Fq, Fq)) -> Curve {
            let (a, b) = O0;

            let aa = a.square();
            let bb = b.square();

            let T1 = &aa + &bb;
            let T2 = &aa - &bb;

            let A = -&(&T1.square() + &T2.square()) / &(&T1 * &T2);

            Curve::new(&A)
        }

        /// Given a dimension one theta point, compute the corresponding
        /// elliptic curve point on the Kummer line (X : Z)
        /// Algorithm from:
        /// Models of Kummer lines and Galois representation,
        /// Razvan Barbulescu, Damien Robert and Nicolas Sarkis
        fn theta_point_to_montgomery_point(O0: &(Fq, Fq), P: &(Fq, Fq)) -> PointX {
            let (a, b) = O0;
            let (U, V) = P;

            let X = a * V + b * U;
            let Z = a * V - b * U;

            // TODO: rather than use this Type we could directly lift here instead...
            // exposing PointX to be public was the easiest way to keep everything from
            // eccore the same, which will help in the future
            PointX::new_xz(&X, &Z)
        }

        /// Given a ThetaStructure and set of points in the compatible form,
        /// compute the product of elliptic curves and affine points on these
        /// curves.
        fn split_to_product(
            Th: &ThetaStructure,
            image_points: &[ThetaPoint],
            num_image_points: usize,
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // First we take the domain theta null point and
            // split this to two level-1 theta null points
            let null_point = Th.null_point();
            let (O1, O2) = split_theta_point(&null_point);

            // Compute Montgomery curve from dimension one
            // null points
            let E3 = null_point_to_montgomery_curve(&O1);
            let E4 = null_point_to_montgomery_curve(&O2);
            let E3E4 = EllipticProduct::new(&E3, &E4);

            // Now compute points on E3 x E4
            let mut C: CouplePoint;
            let mut couple_points: Vec<CouplePoint> = vec![];

            for P in image_points.iter().take(num_image_points) {
                // Split to level 1
                let (P1, P2) = split_theta_point(P);
                // Compute the XPoint (X : Z) from each theta point
                let Q1X = theta_point_to_montgomery_point(&O1, &P1);
                let Q2X = theta_point_to_montgomery_point(&O2, &P2);

                // Lift these points to (X : Y : Z) on the curves
                let (Q1, _) = E3.complete_pointX(&Q1X);
                let (Q2, _) = E4.complete_pointX(&Q2X);

                // Package these points into a CouplePoint on
                // E3 x E4
                C = CouplePoint::new(&Q1, &Q2);

                // Push this into the output
                couple_points.push(C);
            }

            (E3E4, couple_points)
        }

        // ========================================================
        // Main Method! Compute the isogeny between elliptic
        // products
        // ========================================================

        // A general comment about "optimal strategies" -- For the isogeny we
        // have a kernel of order 2^n, and to make a step on the (2,2)-isogeny
        // graph what we want to do is scale this point to get 2^(n-1) * P,
        // which is a point of order two, which then allows us to compute the
        // codomain and images. However, as doubling is more expensive than
        // computing an image (in the general case) it's worth storing many
        // values along the way while doubling and pushing them all through the
        // isogeny. As these pushed through points have the same order, the
        // subsequent steps will need to be doubled less (at the cost of needing
        // to push through more points.)
        //
        // For a proper reference, see Sec 4.2 of
        // https://eprint.iacr.org/2011/506.pdf
        //
        // Gluing: Doubling cost: 16M 16S (We have to double two elliptic curve
        // points) Image cost: 76M + 18S + 1I (Very expensive!)
        //
        // All other steps: Doubling cost: 8S + 6M Image cost: 4S + 3M
        //
        // So ignoring the gluing step, we see images have 1/2 the cost
        // (mathematically this is expected as our doubling formula is
        // essentially just two isogeny images) and so the optimised strategy
        // is computed with a weight that doubling is 2X images.
        //
        // For a function to see how strategies are computed, see strategy.py
        // The current implementation "knows" that the gluing is more expensive
        // and so has extra costs for the leftmost branch of the tree.

        /// Compute an isogeny between elliptic products, naive method with no
        /// optimised strategy. Only here for benchmarking
        pub fn product_isogeny_no_strategy(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
            image_points: &[CouplePoint],
            n: usize,
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // Store the number of image points we wish to evaluate to
            // ensure we return them all from the points we push through
            let num_image_points = image_points.len();

            // Convert the &[...] to a vector so we can add points to this
            // dynamically during the optimal strategy
            let mut kernel_couple_pts = image_points.to_vec();

            // Include the kernel inside the vector of points
            // to evaluate. At each step, every element of the
            // vector should be evaluated
            kernel_couple_pts.push(*P1P2);
            kernel_couple_pts.push(*Q1Q2);

            // Compute points of order 8
            let P1P2_8 = E1E2.double_iter(&P1P2, n - 1);
            let Q1Q2_8 = E1E2.double_iter(&Q1Q2, n - 1);

            // Compute Gluing isogeny
            let (mut domain, mut kernel_pts) =
                gluing_isogeny(&E1E2, &P1P2_8, &Q1Q2_8, &kernel_couple_pts);

            // Do all remaining steps
            let mut Tp1: ThetaPoint;
            let mut Tp2: ThetaPoint;
            for k in 1..n {
                // Repeatedly double to obtain points in the 8-torsion below the kernel
                Tp1 = domain.double_iter(&kernel_pts[num_image_points], n - k - 1);
                Tp2 = domain.double_iter(&kernel_pts[num_image_points + 1], n - k - 1);

                // For the last two steps, we need to be careful because of the zero-null
                // coordinates appearing from the product structure. To avoid these, we
                // use the hadamard transform to avoid them,
                if k == (n - 2) {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, false])
                } else if k == (n - 1) {
                    domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
                } else {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, true])
                }
            }

            // Use a symplectic transform to first get the domain into a compatible form
            // for splitting
            domain = splitting_isomorphism(domain, &mut kernel_pts);

            // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
            // onto this product
            let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

            (product, couple_points)
        }

        /// Compute an isogeny between elliptic products, use an optimised
        /// strategy for all steps assuming doubling is always more expensive
        /// that images, which is not true for gluing.
        pub fn product_isogeny(
            E1E2: &EllipticProduct,
            P1P2: &CouplePoint,
            Q1Q2: &CouplePoint,
            image_points: &[CouplePoint],
            n: usize,
            strategy: &[usize],
        ) -> (EllipticProduct, Vec<CouplePoint>) {
            // Store the number of image points we wish to evaluate to
            // ensure we return them all from the points we push through
            let num_image_points = image_points.len();

            // Convert the &[...] to a vector so we can add points to this
            // dynamically during the optimal strategy
            let mut kernel_couple_pts = image_points.to_vec();

            // Include the kernel inside the vector of points
            // to evaluate. At each step, every element of the
            // vector should be evaluated
            kernel_couple_pts.push(*P1P2);
            kernel_couple_pts.push(*Q1Q2);

            // Bookkeeping for optimised strategy
            let mut strat_idx = 0;
            let mut level: Vec<usize> = vec![0];
            let mut prev: usize = level.iter().sum();

            // =======================================================
            // Gluing Step
            // TODO:
            // Because of type differences there's annoying code reuse
            // for the optimal strategy here and again for every step
            // in the chain thereafter. Which is bothersome. Maybe there
            // is a better way to write this...
            // =======================================================
            let mut ker1 = *P1P2;
            let mut ker2 = *Q1Q2;

            while prev != (n - 1) {
                // Add the next strategy to the level
                level.push(strategy[strat_idx]);

                // Double the points according to the strategy
                ker1 = E1E2.double_iter(&ker1, strategy[strat_idx]);
                ker2 = E1E2.double_iter(&ker2, strategy[strat_idx]);

                // Add these points to the image points
                kernel_couple_pts.push(ker1);
                kernel_couple_pts.push(ker2);

                // Update the strategy bookkeepping
                prev += strategy[strat_idx];
                strat_idx += 1;
            }

            // Clear out the used kernel point and update level
            kernel_couple_pts.pop();
            kernel_couple_pts.pop();
            level.pop();

            // Compute Gluing isogeny
            let (mut domain, mut kernel_pts) =
                gluing_isogeny(&E1E2, &ker1, &ker2, &kernel_couple_pts);

            // ======================================================
            // All other steps
            // Compute the (2^n-1, 2^n-1)-chain in the theta model
            // =======================================================

            let mut Tp1: ThetaPoint;
            let mut Tp2: ThetaPoint;
            let mut kernel_len: usize;

            // Do all remaining steps
            for k in 1..n {
                prev = level.iter().sum();
                kernel_len = kernel_pts.len();

                Tp1 = kernel_pts[kernel_len - 2];
                Tp2 = kernel_pts[kernel_len - 1];

                while prev != (n - 1 - k) {
                    // Add the next strategy to the level
                    level.push(strategy[strat_idx]);

                    // Double the points according to the strategy
                    Tp1 = domain.double_iter(&Tp1, strategy[strat_idx]);
                    Tp2 = domain.double_iter(&Tp2, strategy[strat_idx]);

                    // Add these points to the image points
                    kernel_pts.push(Tp1);
                    kernel_pts.push(Tp2);

                    // Update the strategy bookkeepping
                    prev += strategy[strat_idx];
                    strat_idx += 1;
                }

                // Clear out the used kernel point and update level
                kernel_pts.pop();
                kernel_pts.pop();
                level.pop();

                // For the last two steps, we need to be careful because of the zero-null
                // coordinates appearing from the product structure. To avoid these, we
                // use the hadamard transform to avoid them,
                if k == (n - 2) {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, false])
                } else if k == (n - 1) {
                    domain = two_isogeny_to_product(&Tp1, &Tp2, &mut kernel_pts)
                } else {
                    domain = two_isogeny(&domain, &Tp1, &Tp2, &mut kernel_pts, [false, true])
                }
            }

            // Use a symplectic transform to first get the domain into a compatible form
            // for splitting
            domain = splitting_isomorphism(domain, &mut kernel_pts);

            // Split from the level 2 theta model to the elliptic product E3 x E4 and map points
            // onto this product
            let (product, couple_points) = split_to_product(&domain, &kernel_pts, num_image_points);

            (product, couple_points)
        }
    };
} // End of macro: define_theta_structure

pub(crate) use define_theta_structure;
