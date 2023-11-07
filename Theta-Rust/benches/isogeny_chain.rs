// Benchmark the performance of the isogeny chain of length 208 between elliptic
// products with a base field Fp^2 with characteristic of 381 bits. This is an
// isogeny where the base field has approx 3*128 bits, which is of interest to
// certain protocols and should help with approximations of benchmarks.

#![allow(non_snake_case)]

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_rs::theta381::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};

fn criterion_product_chain(c: &mut Criterion) {
    let A1 = Fq::ZERO;
    let A2_str = "346bc58bbc832c1c6c4c8b6cfeb58a483385fcb98594ac4d22e85ce21a520f12138a05fa5cb716068ef218ee1fefb012476841a87b5755af5d1c7041ce205487d5a98922d6181f68b8ec6a41058b8501baab5bc935017509142d7e1a212b2209";
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "3a1111b79a2ac3a148c08c10cfe6f67ad7a9bd16aa56458013a99f266cbe4146e183bbc865881aad06022551b0e87f0aac35b5ca7f605e29cf7fa11f130a77a63ccc7a3930ffb0624af8856c829f90f70b21a87295b14b7e3b8cec6b23fd3f06";
    let P1_Y_str = "8cd2e58f9fa3f3678914731daf18d4fc33ac1fc491dfcd2bb0f51f4029627bd4e324205cda127a52a4fbab7d3979920b215ad30c48034e000579fecc8e26ec908c4348c98b0a2fbce0a54556e6d3489db11bfd7e2804077e18448f7dc945d008";
    let P2_X_str = "c055d51a89fbd10363d07613b1a71c925750811a11207d583642b71dd296c217e576205442b87c7d0ecf00f29d67e90b8f6898353649e1d53cd75601a1a5148384fa2a9fc89084746fc0d10b68e3d74e0f1ec66c74a582f98ffa8ba2497e080b";
    let P2_Y_str = "542497699a437f71c004bb51a0033e5199257811be7b397b84a6195532f8ec90de59953a1cd9b17bad40fa5c7897e5086ca8159f41a867455f53a1ea6766319bd3286eae4c30c37c1bdedbaa873ab854962f6526f82469b23d73194e1d9f9908";
    let Q1_X_str = "aed83312515b93af8e0fda085d3f7f96c455fe3a7a60e8c45239abb008f9c483eb61225a4ff0410be385c559766df60e28bae5c902ee1c46685717687ab6b5d67e5a87402092578d47842cd2ff943104ca7adff59544f7124d011919c564d607";
    let Q1_Y_str = "f8edffdaa53ac760c85b496ed64acdda15aed38d4b1e729f171caa77f482c406db35ac53cb94b8f9631d66c7b1b75310735f6a1d09815dadf8a1f849cdbee90f8d0cd0d2ef77f1a4b7ea303f00454f7b718e90433a487bc1f7e75d1d63106911";
    let Q2_X_str = "d9e22eb75836ea5dd0e5f08cc939b63976c948f6c1bf71b3568598c7c752df2a4ea0768d978d62895c584b1588fc9c111f36bacc042a17fea60fdfeb10608120a7afad2f61064b7fb7dd4e688089348f00d51e5182c51a88bfeff67f68564909";
    let Q2_Y_str = "79c00be165d30dcc11c5a46f512cdec224a09dc3213a0834318aba1ca33e99d4e1b992e91ec91faf714a85e7a1cfe9138e36f63a53f494223dc07a38f79c9edb090806959f0ad779b432768a9274ede5db92390a9040fd1b57e78baee33fc405";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let PA_X_str = "49c59ebe43e505e387d8306fcba301dbec0499bb5a0518e06c1feb4cabf78a8d5588089f6f7f724776e484efd54f5f0270399cd2318d17876bd7a6e0b75891d0b2b4489a85af988b7518f8b05bb8bdeb2fbc2433dfb739861ca2471fccef4d0c";
    let PA_Y_str = "47d0dde303fbc0e86869647e6c28e82bb82250761f5dce33bc4bb1a895721832270e7a063d98fb1eb1bade0760cf3811269db26ba866299486ae92a8182afcccbacaa6bce40d2905949ff30ef598e1f2a3f1ba01aad7aaac98f256263da0c70b";
    let (PA_X, _) = Fq::decode(&hex::decode(PA_X_str).unwrap());
    let (PA_Y, _) = Fq::decode(&hex::decode(PA_Y_str).unwrap());

    // Curves which define elliptic product
    let E1 = Curve::new(&A1);
    let E2 = Curve::new(&A2);
    let E1E2 = EllipticProduct::new(&E1, &E2);

    // Points on E1 x E2
    let P1 = Point::new_xy(&P1_X, &P1_Y);
    let P2 = Point::new_xy(&P2_X, &P2_Y);
    let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
    let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
    let P1P2 = CouplePoint::new(&P1, &P2);
    let Q1Q2 = CouplePoint::new(&Q1, &Q2);

    // Point to push through isogeny
    let PA = Point::INFINITY;
    let PB = Point::new_xy(&PA_X, &PA_Y);
    let PAPB = CouplePoint::new(&PA, &PB);
    let image_points = [];

    // Length of isogeny chain
    let n = 208;

    // Precomputed with strategy.py
    let strategy: [usize; 207] = [
        84, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 29, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 8, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
    ];

    c.bench_function(
        "Product isogeny of 208 steps with optimised strat with 381-bit prime",
        |b| b.iter(|| product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy)),
    );
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(15));
    targets = criterion_product_chain
}
criterion_main!(benches);
