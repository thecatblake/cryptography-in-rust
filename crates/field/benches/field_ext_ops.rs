use bigint::U256;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use field::{Fp, FpMontConfig, MontBackend};
use field_core::{Field, Fp2, Fp2Config, Fp6, Fp6Config, Fp12, Fp12Config, FpBackend};

// secp256k1 base field prime, p = 2^256 - 2^32 - 977. Same modulus/R2/
// N_PRIME as field_ops.rs's Secp256k1PMont fixture (see there for the
// derivation of R2/N_PRIME); this bench only adds the Fp2/Fp6/Fp12 tower on
// top, so it's redefined here rather than shared, matching field_ops.rs's
// own convention of keeping curve fixtures local to the bench that uses
// them (the field crate itself stays modulus-agnostic).
struct Secp256k1PMont;
impl FpMontConfig for Secp256k1PMont {
    const MODULUS: U256 = U256::from_limbs([
        0xffff_fffe_ffff_fc2f,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
    ]);
    const R2: U256 = U256::from_limbs([0x0000_07a2_000e_90a1, 1, 0, 0]);
    const N_PRIME: U256 = U256::from_limbs([
        0xd838_091d_d225_3531,
        0xbcb2_23fe_dc24_a059,
        0x9c46_c2c2_95f2_b761,
        0xc9bd_1905_1553_8399,
    ]);
}

type Backend = MontBackend<Secp256k1PMont>;
type Fq = Fp<Backend>;

// Fp2/Fp6/Fp12 tower over the secp256k1 base field. Every constant below is
// the Montgomery form (value * R mod p, R = 2^256) of a canonical value
// found by the accompanying offline derivation script (small-non-residue
// search + modular exponentiation over the exact same QuadExt/CubicExt
// tower arithmetic field-core builds on) -- same "verified offline"
// convention as mersenne31/goldilocks/babybear's towers, except here the
// Montgomery conversion itself also has to happen offline: field-core's
// to_mont_uNN/pow_mont_uNN compile-time helpers only exist for primitive
// uNN reprs, not U256, so there's no compile-time path available for a
// 256-bit modulus the way babybear gets one.
impl Fp2Config for Secp256k1PMont {
    type Base = Fq;

    // Canonical value: 3 (smallest quadratic non-residue mod p).
    const BETA: Self::Base = Fp::new(U256::from_limbs([
        0x0000_0003_0000_0b73,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0000,
        0x0000_0000_0000_0000,
    ]));

    // Canonical value: p - 1 (= BETA^((p-1)/2) mod p).
    const FROBENIUS_COEFF: Self::Base = Fp::new(U256::from_limbs([
        0xffff_fffd_ffff_f85e,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
    ]));
}

impl Fp6Config for Secp256k1PMont {
    // Canonical value: u (smallest cubic non-residue in Fp2 found by the
    // search).
    const XI: Fp2<Self> = Fp2 {
        c0: Fp::new(U256::ZERO),
        c1: Fp::new(U256::from_limbs([
            0x0000_0001_0000_03d1,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
            0x0000_0000_0000_0000,
        ])),
    };

    // Canonical value: XI^((p-1)/3) in Fp2.
    const FROBENIUS_COEFF_C1: Fp2<Self> = Fp2 {
        c0: Fp::new(U256::from_limbs([
            0xa75b_c9e2_717e_72e1,
            0xfc02_1e9c_e3b4_7f50,
            0x0716_7687_2fd1_c6fa,
            0x85b5_c951_4344_c2ac,
        ])),
        c1: Fp::new(U256::ZERO),
    };
    // Canonical value: XI^(2(p-1)/3) in Fp2.
    const FROBENIUS_COEFF_C2: Fp2<Self> = Fp2 {
        c0: Fp::new(U256::from_limbs([
            0xa75b_c9e1_717e_6f10,
            0xfc02_1e9c_e3b4_7f50,
            0x0716_7687_2fd1_c6fa,
            0x85b5_c951_4344_c2ac,
        ])),
        c1: Fp::new(U256::ZERO),
    };
}

impl Fp12Config for Secp256k1PMont {
    // Canonical value: v (Fp6's own generator, embedded via c1 with
    // c0 = c2 = 0), the smallest quadratic non-residue in Fp6 found by the
    // search.
    const GAMMA: Fp6<Self> = Fp6 {
        c0: Fp2 { c0: Fp::new(U256::ZERO), c1: Fp::new(U256::ZERO) },
        c1: Fp2 {
            c0: Fp::new(U256::from_limbs([
                0x0000_0001_0000_03d1,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
                0x0000_0000_0000_0000,
            ])),
            c1: Fp::new(U256::ZERO),
        },
        c2: Fp2 { c0: Fp::new(U256::ZERO), c1: Fp::new(U256::ZERO) },
    };

    // Canonical value: GAMMA^((p-1)/2) in Fp6.
    const FROBENIUS_COEFF: Fp6<Self> = Fp6 {
        c0: Fp2 {
            c0: Fp::new(U256::ZERO),
            c1: Fp::new(U256::from_limbs([
                0x37c9_434b_7b2a_264b,
                0xfeab_5f89_a13c_2a70,
                0x025c_d22d_0ff0_97a8,
                0x81e7_431b_166c_40e4,
            ])),
        },
        c1: Fp2 { c0: Fp::new(U256::ZERO), c1: Fp::new(U256::ZERO) },
        c2: Fp2 { c0: Fp::new(U256::ZERO), c1: Fp::new(U256::ZERO) },
    };
}

type Fp2Q = Fp2<Secp256k1PMont>;
type Fp6Q = Fp6<Secp256k1PMont>;
type Fp12Q = Fp12<Secp256k1PMont>;

// Seeds a base-field element in Montgomery form via the public mul/R2 API,
// without needing MontBackend's private to_mont helper -- same trick
// field_ops.rs's bench_mont uses.
fn fq(v: u64) -> Fq {
    Fp::new(Backend::mul(U256::from(v), Secp256k1PMont::R2))
}

fn bench_fp2(c: &mut Criterion) {
    let a = Fp2Q::new(fq(0x1234_5678_9abc_def0), fq(0x0fed_cba9_8765_4321));
    let b = Fp2Q::new(fq(0xfedc_ba98_7654_3210), fq(0x1122_3344_5566_7788));

    let mut group = c.benchmark_group("Secp256k1MontFp2");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp6(c: &mut Criterion) {
    let a = Fp6Q::new(Fp2Q::new(fq(1), fq(2)), Fp2Q::new(fq(3), fq(4)), Fp2Q::new(fq(5), fq(6)));
    let b = Fp6Q::new(Fp2Q::new(fq(7), fq(8)), Fp2Q::new(fq(9), fq(10)), Fp2Q::new(fq(11), fq(12)));

    let mut group = c.benchmark_group("Secp256k1MontFp6");
    group.bench_function("add", |bch| bch.iter(|| black_box(a) + black_box(b)));
    group.bench_function("sub", |bch| bch.iter(|| black_box(a) - black_box(b)));
    group.bench_function("mul", |bch| bch.iter(|| black_box(a) * black_box(b)));
    group.bench_function("square", |bch| bch.iter(|| black_box(a).square()));
    group.bench_function("neg", |bch| bch.iter(|| -black_box(a)));
    group.bench_function("inverse", |bch| bch.iter(|| black_box(a).inverse()));
    group.finish();
}

fn bench_fp12(c: &mut Criterion) {
    let a = Fp12Q::new(
        Fp6Q::new(Fp2Q::new(fq(1), fq(2)), Fp2Q::new(fq(3), fq(4)), Fp2Q::new(fq(5), fq(6))),
        Fp6Q::new(Fp2Q::new(fq(7), fq(8)), Fp2Q::new(fq(9), fq(10)), Fp2Q::new(fq(11), fq(12))),
    );
    let b = Fp12Q::new(
        Fp6Q::new(Fp2Q::new(fq(13), fq(14)), Fp2Q::new(fq(15), fq(16)), Fp2Q::new(fq(17), fq(18))),
        Fp6Q::new(Fp2Q::new(fq(19), fq(20)), Fp2Q::new(fq(21), fq(22)), Fp2Q::new(fq(23), fq(24))),
    );

    let mut group = c.benchmark_group("Secp256k1MontFp12");
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
