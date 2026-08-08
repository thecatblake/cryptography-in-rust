use criterion::{black_box, criterion_group, criterion_main, Criterion};
use mersenne31::{FpBackend, Mersenne31, Mersenne31Backend};

fn bench_mersenne31(c: &mut Criterion) {
    let a = Mersenne31::new(0x1234_5678 % Mersenne31Backend::MODULUS);
    let b = Mersenne31::new(0x7edc_ba98 % Mersenne31Backend::MODULUS);
    let exp = 0x1234_5678u32;

    let mut group = c.benchmark_group("Mersenne31Backend");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.bench_function("pow", |bch| bch.iter(|| black_box(a).pow(black_box(exp))));
    group.finish();
}

criterion_group!(benches, bench_mersenne31);
criterion_main!(benches);
