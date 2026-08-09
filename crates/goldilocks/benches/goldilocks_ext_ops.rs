use criterion::{black_box, criterion_group, criterion_main, Criterion};
use goldilocks::{
    Field, FpBackend, Goldilocks, GoldilocksBackend, GoldilocksFp2, GoldilocksFp6, GoldilocksFp12,
};

fn fg(v: u64) -> Goldilocks {
    Goldilocks::new(v % GoldilocksBackend::MODULUS)
}

fn bench_fp2(c: &mut Criterion) {
    let a = GoldilocksFp2::new(fg(0x1234_5678_9abc_def0), fg(0x0fed_cba9_8765_4321));
    let b = GoldilocksFp2::new(fg(0xfedc_ba98_7654_3210), fg(0x1234_5678_9abc_def0));

    let mut group = c.benchmark_group("GoldilocksFp2");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp6(c: &mut Criterion) {
    let a = GoldilocksFp6::new(
        GoldilocksFp2::new(fg(1), fg(2)),
        GoldilocksFp2::new(fg(3), fg(4)),
        GoldilocksFp2::new(fg(5), fg(6)),
    );
    let b = GoldilocksFp6::new(
        GoldilocksFp2::new(fg(7), fg(8)),
        GoldilocksFp2::new(fg(9), fg(10)),
        GoldilocksFp2::new(fg(11), fg(12)),
    );

    let mut group = c.benchmark_group("GoldilocksFp6");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp12(c: &mut Criterion) {
    let a = GoldilocksFp12::new(
        GoldilocksFp6::new(
            GoldilocksFp2::new(fg(1), fg(2)),
            GoldilocksFp2::new(fg(3), fg(4)),
            GoldilocksFp2::new(fg(5), fg(6)),
        ),
        GoldilocksFp6::new(
            GoldilocksFp2::new(fg(7), fg(8)),
            GoldilocksFp2::new(fg(9), fg(10)),
            GoldilocksFp2::new(fg(11), fg(12)),
        ),
    );
    let b = GoldilocksFp12::new(
        GoldilocksFp6::new(
            GoldilocksFp2::new(fg(13), fg(14)),
            GoldilocksFp2::new(fg(15), fg(16)),
            GoldilocksFp2::new(fg(17), fg(18)),
        ),
        GoldilocksFp6::new(
            GoldilocksFp2::new(fg(19), fg(20)),
            GoldilocksFp2::new(fg(21), fg(22)),
            GoldilocksFp2::new(fg(23), fg(24)),
        ),
    );

    let mut group = c.benchmark_group("GoldilocksFp12");
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
