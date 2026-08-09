use babybear::{BabyBear, BabyBearBackend, BabyBearFp2, BabyBearFp6, BabyBearFp12, Field, FpBackend};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn fb(v: u32) -> BabyBear {
    BabyBear::new(BabyBearBackend::to_mont(v % BabyBearBackend::MODULUS))
}

fn bench_fp2(c: &mut Criterion) {
    let a = BabyBearFp2::new(fb(0x1234_5678), fb(0x0abc_def0));
    let b = BabyBearFp2::new(fb(0x7edc_ba98), fb(0x0fed_cba9));

    let mut group = c.benchmark_group("BabyBearFp2");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp6(c: &mut Criterion) {
    let a = BabyBearFp6::new(
        BabyBearFp2::new(fb(1), fb(2)),
        BabyBearFp2::new(fb(3), fb(4)),
        BabyBearFp2::new(fb(5), fb(6)),
    );
    let b = BabyBearFp6::new(
        BabyBearFp2::new(fb(7), fb(8)),
        BabyBearFp2::new(fb(9), fb(10)),
        BabyBearFp2::new(fb(11), fb(12)),
    );

    let mut group = c.benchmark_group("BabyBearFp6");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp12(c: &mut Criterion) {
    let a = BabyBearFp12::new(
        BabyBearFp6::new(
            BabyBearFp2::new(fb(1), fb(2)),
            BabyBearFp2::new(fb(3), fb(4)),
            BabyBearFp2::new(fb(5), fb(6)),
        ),
        BabyBearFp6::new(
            BabyBearFp2::new(fb(7), fb(8)),
            BabyBearFp2::new(fb(9), fb(10)),
            BabyBearFp2::new(fb(11), fb(12)),
        ),
    );
    let b = BabyBearFp12::new(
        BabyBearFp6::new(
            BabyBearFp2::new(fb(13), fb(14)),
            BabyBearFp2::new(fb(15), fb(16)),
            BabyBearFp2::new(fb(17), fb(18)),
        ),
        BabyBearFp6::new(
            BabyBearFp2::new(fb(19), fb(20)),
            BabyBearFp2::new(fb(21), fb(22)),
            BabyBearFp2::new(fb(23), fb(24)),
        ),
    );

    let mut group = c.benchmark_group("BabyBearFp12");
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
