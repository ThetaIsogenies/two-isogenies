// Fairly ugly-hack benchmark to compute the costs of multiplications, squares
// and inversions for the three base fields. These values are needed to compute
// the correct optimisation costs for the strategies.

use theta_rs::fields::Fp254Ext::Fp2 as Fp2Sml;
use theta_rs::fields::Fp381Ext::Fp2 as Fp2Med;
use theta_rs::fields::FpFESTAExt::Fp2 as Fp2Big;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use std::time::Duration;

mod util;
use util::DRNG;

fn m_sml(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Sml = Fp2Sml::rand(&mut rng);
    let y: Fp2Sml = Fp2Sml::rand(&mut rng);

    c.bench_function("New Multiplication (254 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn old_m_sml(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Sml = Fp2Sml::rand(&mut rng);
    let y: Fp2Sml = Fp2Sml::rand(&mut rng);

    c.bench_function("Old Multiplication (254 bit)", |b| {
        b.iter(|| black_box(x).mul_old(black_box(y)))
    });
}

fn s_sml(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Sml = Fp2Sml::rand(&mut rng);

    c.bench_function("Square (254 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_sml(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Sml = Fp2Sml::rand(&mut rng);

    c.bench_function("Invert (254 bit)", |b| b.iter(|| black_box(x).invert()));
}

fn m_med(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Med = Fp2Med::rand(&mut rng);
    let y: Fp2Med = Fp2Med::rand(&mut rng);

    c.bench_function("New Multiplication (381 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn old_m_med(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Med = Fp2Med::rand(&mut rng);
    let y: Fp2Med = Fp2Med::rand(&mut rng);

    c.bench_function("Old Multiplication (381 bit)", |b| {
        b.iter(|| black_box(x).mul_old(black_box(y)))
    });
}

fn s_med(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Med = Fp2Med::rand(&mut rng);

    c.bench_function("Square (381 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_med(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Med = Fp2Med::rand(&mut rng);

    c.bench_function("Invert (381 bit)", |b| b.iter(|| black_box(x).invert()));
}

fn m_festa(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Big = Fp2Big::rand(&mut rng);
    let y: Fp2Big = Fp2Big::rand(&mut rng);

    c.bench_function("New Multiplication (1293 bit)", |b| {
        b.iter(|| black_box(x).mul_new(black_box(y)))
    });
}

fn old_m_festa(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Big = Fp2Big::rand(&mut rng);
    let y: Fp2Big = Fp2Big::rand(&mut rng);

    c.bench_function("Old Multiplication (1293 bit)", |b| {
        b.iter(|| black_box(x).mul_old(black_box(y)))
    });
}

fn s_festa(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Big = Fp2Big::rand(&mut rng);

    c.bench_function("Square (1293 bit)", |b| b.iter(|| black_box(x).square()));
}

fn i_festa(c: &mut Criterion) {
    let mut rng = DRNG::new();

    let x: Fp2Big = Fp2Big::rand(&mut rng);

    c.bench_function("Invert (1293 bit)", |b| b.iter(|| black_box(x).invert()));
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(3));
    targets = m_sml, old_m_sml, s_sml, i_sml, m_med, old_m_med, s_med, i_med, m_festa, old_m_festa, s_festa, i_festa
}
criterion_main!(benches);
