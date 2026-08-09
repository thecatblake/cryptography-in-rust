use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mersenne31::{
    Field, FpBackend, Mersenne31, Mersenne31Backend, Mersenne31Fp2, Mersenne31Fp6, Mersenne31Fp12,
};

fn fm(v: u32) -> Mersenne31 {
    Mersenne31::new(v % Mersenne31Backend::MODULUS)
}

fn bench_fp2(c: &mut Criterion) {
    let a = Mersenne31Fp2::new(fm(0x1234_5678), fm(0x0abc_def0));
    let b = Mersenne31Fp2::new(fm(0x7edc_ba98), fm(0x0fed_cba9));

    let mut group = c.benchmark_group("Mersenne31Fp2");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp6(c: &mut Criterion) {
    let a = Mersenne31Fp6::new(
        Mersenne31Fp2::new(fm(1), fm(2)),
        Mersenne31Fp2::new(fm(3), fm(4)),
        Mersenne31Fp2::new(fm(5), fm(6)),
    );
    let b = Mersenne31Fp6::new(
        Mersenne31Fp2::new(fm(7), fm(8)),
        Mersenne31Fp2::new(fm(9), fm(10)),
        Mersenne31Fp2::new(fm(11), fm(12)),
    );

    let mut group = c.benchmark_group("Mersenne31Fp6");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp12(c: &mut Criterion) {
    let a = Mersenne31Fp12::new(
        Mersenne31Fp6::new(
            Mersenne31Fp2::new(fm(1), fm(2)),
            Mersenne31Fp2::new(fm(3), fm(4)),
            Mersenne31Fp2::new(fm(5), fm(6)),
        ),
        Mersenne31Fp6::new(
            Mersenne31Fp2::new(fm(7), fm(8)),
            Mersenne31Fp2::new(fm(9), fm(10)),
            Mersenne31Fp2::new(fm(11), fm(12)),
        ),
    );
    let b = Mersenne31Fp12::new(
        Mersenne31Fp6::new(
            Mersenne31Fp2::new(fm(13), fm(14)),
            Mersenne31Fp2::new(fm(15), fm(16)),
            Mersenne31Fp2::new(fm(17), fm(18)),
        ),
        Mersenne31Fp6::new(
            Mersenne31Fp2::new(fm(19), fm(20)),
            Mersenne31Fp2::new(fm(21), fm(22)),
            Mersenne31Fp2::new(fm(23), fm(24)),
        ),
    );

    let mut group = c.benchmark_group("Mersenne31Fp12");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

criterion_group!(benches, bench_fp2, bench_fp6, bench_fp12);
criterion_main!(benches);
