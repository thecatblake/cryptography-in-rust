use babybear::{BabyBear, BabyBearBackend, FpBackend};
use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn bench_babybear(c: &mut Criterion) {
    let a = BabyBear::new(0x1234_5678 % BabyBearBackend::MODULUS);
    let b = BabyBear::new(0x0abc_def0 % BabyBearBackend::MODULUS);
    let exp = 0x1234_5678u32;

    let mut group = c.benchmark_group("BabyBearBackend");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.bench_function("pow", |bch| bch.iter(|| black_box(a).pow(black_box(exp))));
    group.finish();
}

criterion_group!(benches, bench_babybear);
criterion_main!(benches);
