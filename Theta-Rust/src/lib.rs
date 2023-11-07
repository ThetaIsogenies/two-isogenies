#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

#[allow(unused_macros)]
macro_rules! static_assert {
    ($condition:expr) => {
        let _ = &[()][1 - ($condition) as usize];
    };
}

pub mod eccore;
pub mod fields;
pub mod finitefield;
pub mod theta;

pub mod ec254 {
    pub type Fq = crate::fields::Fp254Ext::Fp2;
    crate::eccore::define_ec_core! {}
}

pub mod theta254 {
    pub use crate::ec254::{CouplePoint, Curve, EllipticProduct, Fq, Point, PointX};
    crate::theta::define_theta_structure! {}
}

pub mod ec381 {
    pub type Fq = crate::fields::Fp381Ext::Fp2;
    crate::eccore::define_ec_core! {}
}

pub mod theta381 {
    pub use crate::ec381::{CouplePoint, Curve, EllipticProduct, Fq, Point, PointX};
    crate::theta::define_theta_structure! {}
}

pub mod ecFESTA {
    pub type Fq = crate::fields::FpFESTAExt::Fp2;
    crate::eccore::define_ec_core! {}
}

pub mod thetaFESTA {
    pub use crate::ecFESTA::{CouplePoint, Curve, EllipticProduct, Fq, Point, PointX};
    crate::theta::define_theta_structure! {}
}
