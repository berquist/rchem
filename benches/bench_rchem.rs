#[macro_use]
extern crate criterion;

use criterion::Criterion;

use rchem;

fn criterion_benchmark(c: &mut Criterion) {
    let za = 1.8;
    let zb = 2.8;
    let ra = [0.0, 0.0, 0.0];
    let rb = [0.5, 0.8, -0.2];

    let powers = [0, 0, 0, 0, 0, 0];
    c.bench_function("get_overlap_0_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_overlap(za, zb, &ra, &rb, &powers))
    });
    let powers = [1, 0, 0, 0, 0, 0];
    c.bench_function("get_overlap_1_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_overlap(za, zb, &ra, &rb, &powers))
    });
    let powers = [2, 2, 2, 0, 0, 0];
    c.bench_function("get_overlap_2_2_2_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_overlap(za, zb, &ra, &rb, &powers))
    });
    let powers = [1, 1, 1, 1, 1, 1];
    c.bench_function("get_overlap_1_1_1_1_1_1", move |b| {
        b.iter(|| rchem::integrals::get_overlap(za, zb, &ra, &rb, &powers))
    });
    let powers = [2, 2, 2, 2, 2, 2];
    c.bench_function("get_overlap_2_2_2_2_2_2", move |b| {
        b.iter(|| rchem::integrals::get_overlap(za, zb, &ra, &rb, &powers))
    });

    let powers = [0, 0, 0, 0, 0, 0];
    c.bench_function("get_kinetic_0_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_kinetic(za, zb, &ra, &rb, &powers))
    });
    let powers = [1, 0, 0, 0, 0, 0];
    c.bench_function("get_kinetic_1_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_kinetic(za, zb, &ra, &rb, &powers))
    });
    let powers = [2, 2, 2, 0, 0, 0];
    c.bench_function("get_kinetic_2_2_2_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_kinetic(za, zb, &ra, &rb, &powers))
    });
    let powers = [1, 1, 1, 1, 1, 1];
    c.bench_function("get_kinetic_1_1_1_1_1_1", move |b| {
        b.iter(|| rchem::integrals::get_kinetic(za, zb, &ra, &rb, &powers))
    });
    let powers = [2, 2, 2, 2, 2, 2];
    c.bench_function("get_kinetic_2_2_2_2_2_2", move |b| {
        b.iter(|| rchem::integrals::get_kinetic(za, zb, &ra, &rb, &powers))
    });

    let rc = [0.5, 0.8, 0.2];

    let powers = [0, 0, 0, 0, 0, 0];
    c.bench_function("get_nuclear_0_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_nuclear(za, zb, &ra, &rb, &rc, &powers))
    });
    let powers = [1, 0, 0, 0, 0, 0];
    c.bench_function("get_nuclear_1_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_nuclear(za, zb, &ra, &rb, &rc, &powers))
    });
    let powers = [2, 2, 2, 0, 0, 0];
    c.bench_function("get_nuclear_2_2_2_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_nuclear(za, zb, &ra, &rb, &rc, &powers))
    });
    let powers = [1, 1, 1, 1, 1, 1];
    c.bench_function("get_nuclear_1_1_1_1_1_1", move |b| {
        b.iter(|| rchem::integrals::get_nuclear(za, zb, &ra, &rb, &rc, &powers))
    });
    // let powers = [2, 2, 2, 2, 2, 2];
    // c.bench_function("get_nuclear_2_2_2_2_2_2", move |b| {
    //     b.iter(|| rchem::integrals::get_nuclear(za, zb, &ra, &rb, &rc, &powers))
    // });

    let za = 1.1;
    let zb = 1.2;
    let zc = 1.3;
    let zd = 1.4;

    let ra = [1.0, 0.0, 1.0];
    let rb = [0.0, 1.0, 2.0];
    let rc = [0.0, 0.0, 3.0];
    let rd = [0.0, 0.0, 4.0];

    let powers = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    c.bench_function("get_coulomb_0_0_0_0_0_0_0_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_coulomb(za, zb, zc, zd, &ra, &rb, &rc, &rd, &powers))
    });
    let powers = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    c.bench_function("get_coulomb_1_0_0_0_0_0_0_0_0_0_0_0", move |b| {
        b.iter(|| rchem::integrals::get_coulomb(za, zb, zc, zd, &ra, &rb, &rc, &rd, &powers))
    });
    let powers = [2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0];
    c.bench_function("get_coulomb_2_1_0_1_0_0_1_0_0_0_1_0", move |b| {
        b.iter(|| rchem::integrals::get_coulomb(za, zb, zc, zd, &ra, &rb, &rc, &rd, &powers))
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
