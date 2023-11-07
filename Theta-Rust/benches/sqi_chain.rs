// Benchmark the performance of the isogeny chain of length 126 between elliptic
// products with a base field Fp^2 with characteristic of 254 bits. This is a
// very similar isogeny as that proposed in SQIsignHD, however it is a dimension
// two rather than dimesnion four isogeny so a performance of SQIsignHD can only
// be guessed from these results.

#![allow(non_snake_case)]

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_rs::theta254::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};

fn criterion_product_chain(c: &mut Criterion) {
    let A1 = Fq::ZERO;
    let A2_str = "2687f041b47d1ea9c00ae938b2a761aac6bad9ea80dd9dbe1c24d9ef697d7d0475898661998dd3a7b186e2558d1cf0dd771fb49483988c2ff578547815e8f00e";
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "63385515142ea2f7a0c3747c22668aef99b31f43987354116da1915c24ff2f0a279051962feace976834986fc955b11bb8e1ffea47d8ce994ad22ee86c7f7a00";
    let P1_Y_str = "38243631804d307b72ecd037da591d1f06ca606bf1bb71d77ce10467b00f7b082472068bfea9eb9d80b68e04bb194a23e5214ba41625915d8e590024e5dcf611";
    let P2_X_str = "3e218d8b18cf09ce29cefaa35467225134910411fe33625136e50f7b9c59e51c517629d8786a98603cc06470a10dea83f7eca03b9b378297c21755bf0aee1324";
    let P2_Y_str = "06d56830cea82c91cf4059078566145d4b90d992177916ffe1380060057d75277610f27fc5558e5d028699493a300d84521f0c077c6e52c6adc1820c8f53f11a";
    let Q1_X_str = "3b569c320301eb5aa1aa078b7399d31e2e0e6c70e91223d1d3346be7145c100b16e4c042048452157133174122c04b1c9a17f38c28e959828933a95eebf6a305";
    let Q1_Y_str = "04a2a994555e67a3ddbb5f87dbef7903a9fc9724cb36e51924d28522222c3d2bd07a4d44b06b176b63d78733a1a4839606276b8c523bd6dc8f23de01e2817e17";
    let Q2_X_str = "4076f17bf841c5d6acc194e75a4cd020fe3ea03b0914cf3f3db9cc882f8b4724a18ce4bb13c2b3c46abd8e6cdb502dd7f48e58ee49d6c1d632532f6ed995e12b";
    let Q2_Y_str = "5b88c588b1f27496800acbef34817c5fa5cbcd728de00e31a46fc7aa6ef9af0e173d1ba96e1b2c2ebc5bc3dd3f980344b508b9df1863fb624855dc1a8cc17b12";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let PA_X_str = "95a2c9abfff13d34db8c9d79ecdeb8cf8cbdb2de119f5fc3c51e69dcffab602e79fc55b299bdc279f8c3c20eac06322b43cb71718e367e1f795b4c8fc59fbb0b";
    let PA_Y_str = "b5b5fdfca0cc051994cbeb4338a5e629714df6cc7496cdd98900fdf281ff342990f21fa62afcf762a1c6f3e65635ab87f5c90269722434f2479482b2ec089a05";
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
    let n = 126;

    // Precomputed with strategy.py
    let strategy: [usize; 125] = [
        55, 34, 21, 15, 5, 3, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
        5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1,
        1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3,
        2, 1, 1, 1, 1, 1,
    ];

    c.bench_function(
        "Product isogeny of 126 steps with optimised strat with 254-bit prime",
        |b| b.iter(|| product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy)),
    );
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(15));
    targets = criterion_product_chain
}
criterion_main!(benches);
