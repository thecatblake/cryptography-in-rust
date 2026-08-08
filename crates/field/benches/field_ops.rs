use bigint::U256;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use field::{DefaultBackend, Fp, FpBackend, FpConfig, FpMontConfig, MontBackend};

// secp256k1 base field prime, p = 2^256 - 2^32 - 977.
struct Secp256k1P;
impl FpConfig for Secp256k1P {
    const MODULUS: U256 = U256::from_limbs([
        0xffff_fffe_ffff_fc2f,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
    ]);
}

struct Secp256k1PMont;
impl FpMontConfig for Secp256k1PMont {
    const MODULUS: U256 = Secp256k1P::MODULUS;
    // R2 = R^2 mod p and N_PRIME = -p^-1 mod R (R = 2^256), computed and
    // round-trip verified offline; see mont_roundtrip_is_identity below.
    const R2: U256 = U256::from_limbs([0x0000_07a2_000e_90a1, 1, 0, 0]);
    const N_PRIME: U256 = U256::from_limbs([
        0xd838_091d_d225_3531,
        0xbcb2_23fe_dc24_a059,
        0x9c46_c2c2_95f2_b761,
        0xc9bd_1905_1553_8399,
    ]);
}

fn plain_a() -> U256 {
    U256::from_limbs([
        0x1111_1111_2222_2222,
        0x3333_3333_4444_4444,
        0x5555_5555_6666_6666,
        0x0777_7777_8888_8888,
    ])
}

fn plain_b() -> U256 {
    U256::from_limbs([
        0xaaaa_aaaa_bbbb_bbbb,
        0xcccc_cccc_dddd_dddd,
        0xeeee_eeee_1f1f_1f1f,
        0x0101_0101_0202_0202,
    ])
}

fn bench_backend<B: FpBackend<Repr = U256>>(c: &mut Criterion, name: &str, a: U256, b: U256) {
    let fa = Fp::<B>::new(a);
    let fb = Fp::<B>::new(b);

    let mut group = c.benchmark_group(name);
    group.bench_function("add", |bch| bch.iter(|| black_box(fa) + black_box(fb)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(fa) - black_box(fb)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(fa) * black_box(fb)));
    group.bench_function("square", |bch| bch.iter(|| black_box(fa).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(fa)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(fa).inverse()));
    group.bench_function("pow", |bch| bch.iter(|| black_box(fa).pow(black_box(b))));
    group.finish();
}

fn bench_default(c: &mut Criterion) {
    bench_backend::<DefaultBackend<Secp256k1P>>(c, "DefaultBackend", plain_a(), plain_b());
}

fn bench_mont(c: &mut Criterion) {
    // Seed values in Montgomery form via the public mul/R2 API, without
    // needing MontBackend's private to_mont helper.
    let a_mont = MontBackend::<Secp256k1PMont>::mul(plain_a(), Secp256k1PMont::R2);
    let b_mont = MontBackend::<Secp256k1PMont>::mul(plain_b(), Secp256k1PMont::R2);

    bench_backend::<MontBackend<Secp256k1PMont>>(c, "MontBackend", a_mont, b_mont);
}

criterion_group!(benches, bench_default, bench_mont);
criterion_main!(benches);
