use criterion::{black_box, criterion_group, criterion_main, Criterion};
use goldilocks::{FpBackend, Goldilocks, GoldilocksBackend};

fn bench_goldilocks(c: &mut Criterion) {
    let a = Goldilocks::new(0x1234_5678_9abc_def0 % GoldilocksBackend::MODULUS);
    let b = Goldilocks::new(0xfedc_ba98_7654_3210 % GoldilocksBackend::MODULUS);
    let exp = 0x1234_5678u64;

    let mut group = c.benchmark_group("GoldilocksBackend");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.bench_function("pow", |bch| bch.iter(|| black_box(a).pow(black_box(exp))));
    group.finish();
}

criterion_group!(benches, bench_goldilocks);
criterion_main!(benches);
