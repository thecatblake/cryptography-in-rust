use bigint::U256;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use elliptic_curve::{
    AffinePoint as SwAffinePoint, EdwardsAffinePoint, MontgomeryAffinePoint, MontgomeryCurve,
    ShortWeierstrassCurve, TwistedEdwardsCurve,
};
use field::{DefaultBackend, Fp, FpConfig};

// Same 256-bit modulus as field/benches/field_ops.rs (secp256k1 base field
// prime), reused here purely to get a realistic-sized field for timing --
// the curve constants and point coordinates below are arbitrary, not
// secp256k1's own parameters, since Add/Mul never call validate() and so
// don't need a genuine on-curve point to benchmark the arithmetic.
struct P;
impl FpConfig for P {
    const MODULUS: U256 = U256::from_limbs([
        0xffff_fffe_ffff_fc2f,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
        0xffff_ffff_ffff_ffff,
    ]);
}
type F = Fp<DefaultBackend<P>>;

const fn fe(limbs: [u64; 4]) -> F {
    F::new(U256::from_limbs(limbs))
}

// Two arbitrary, distinct, nonzero points -- distinct x so "add" hits the
// chord branch, nonzero y so "double" (self + self) hits the tangent
// branch rather than the vertical-chord (x1==x2, y1==-y2) exceptional case.
const X1: [u64; 4] = [0x1111_1111_2222_2222, 0x3333_3333_4444_4444, 0x5555_5555_6666_6666, 0x0777_7777_8888_8888];
const Y1: [u64; 4] = [0xaaaa_aaaa_bbbb_bbbb, 0xcccc_cccc_dddd_dddd, 0xeeee_eeee_1f1f_1f1f, 0x0101_0101_0202_0202];
const X2: [u64; 4] = [0x2222_2222_3333_3333, 0x4444_4444_5555_5555, 0x6666_6666_7777_7777, 0x0888_8888_9999_9999];
const Y2: [u64; 4] = [0xbbbb_bbbb_cccc_cccc, 0xdddd_dddd_eeee_eeee, 0x1f1f_1f1f_2020_2020, 0x0202_0202_0303_0303];

// Full-width scalar so scalar_mul's double-and-add walks close to all 256
// bits, matching how field_ops.rs benchmarks `pow` with a full-width
// exponent.
const SCALAR: [u64; 4] = [0xfedc_ba98_7654_3210, 0x0123_4567_89ab_cdef, 0x1357_9bdf_2468_ace0, 0x0f0e_0d0c_0b0a_0908];

fn bench_short_weierstrass(c: &mut Criterion) {
    struct Curve;
    impl ShortWeierstrassCurve for Curve {
        type Field = F;
        type Scalar = F;

        const A: F = fe([7, 0, 0, 0]);
        const B: F = fe([11, 0, 0, 0]);
    }
    type Point = SwAffinePoint<Curve>;

    let p1 = Point { x: fe(X1), y: fe(Y1), infinity: false };
    let p2 = Point { x: fe(X2), y: fe(Y2), infinity: false };
    let scalar = U256::from_limbs(SCALAR);

    let mut group = c.benchmark_group("ShortWeierstrass");
    group.bench_function("add", |b| b.iter(|| black_box(p1) + black_box(p2)));
    group.bench_function("double", |b| b.iter(|| black_box(p1) + black_box(p1)));
    group.bench_function("scalar_mul", |b| b.iter(|| black_box(p1) * black_box(scalar)));
    group.finish();
}

fn bench_twisted_edwards(c: &mut Criterion) {
    struct Curve;
    impl TwistedEdwardsCurve for Curve {
        type Field = F;
        type Scalar = F;

        const A: F = fe([2, 0, 0, 0]);
        const D: F = fe([5, 0, 0, 0]);
    }
    type Point = EdwardsAffinePoint<Curve>;

    let p1 = Point { x: fe(X1), y: fe(Y1) };
    let p2 = Point { x: fe(X2), y: fe(Y2) };
    let scalar = U256::from_limbs(SCALAR);

    // No separate "double": twisted Edwards addition is unified, so
    // doubling is the same code path as add (see twisted_edwards.rs) --
    // benchmarking it again under a different name would just duplicate
    // the "add" number.
    let mut group = c.benchmark_group("TwistedEdwards");
    group.bench_function("add", |b| b.iter(|| black_box(p1) + black_box(p2)));
    group.bench_function("scalar_mul", |b| b.iter(|| black_box(p1) * black_box(scalar)));
    group.finish();
}

fn bench_montgomery(c: &mut Criterion) {
    struct Curve;
    impl MontgomeryCurve for Curve {
        type Field = F;
        type Scalar = F;

        const A: F = fe([3, 0, 0, 0]);
        const B: F = fe([1, 0, 0, 0]);
    }
    type Point = MontgomeryAffinePoint<Curve>;

    let p1 = Point { x: fe(X1), y: fe(Y1), infinity: false };
    let p2 = Point { x: fe(X2), y: fe(Y2), infinity: false };
    let scalar = U256::from_limbs(SCALAR);

    let mut group = c.benchmark_group("Montgomery");
    group.bench_function("add", |b| b.iter(|| black_box(p1) + black_box(p2)));
    group.bench_function("double", |b| b.iter(|| black_box(p1) + black_box(p1)));
    group.bench_function("scalar_mul", |b| b.iter(|| black_box(p1) * black_box(scalar)));
    group.finish();
}

criterion_group!(benches, bench_short_weierstrass, bench_twisted_edwards, bench_montgomery);
criterion_main!(benches);
