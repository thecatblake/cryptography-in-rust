use bigint::U256;
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use elliptic_curve::{AffinePoint, JacobianPoint};
use field::{DefaultBackend, Fp, FpConfig};
use secp256k1::{Secp256k1, Secp256k1Field, SECP256K1_P};

// Same secp256k1 base field modulus as Secp256k1Field, but wired through
// the generic, division-based DefaultBackend instead of FieldConfig's
// fold-based fast reduction (see secp256k1::secp256k1_reduce512). Same
// modulus on both sides, so the "Field (fold)" vs. "Field (division)"
// groups below isolate exactly the payoff of the hand-written `mul`.
struct DivisionConfig;
impl FpConfig for DivisionConfig {
    const MODULUS: U256 = SECP256K1_P;
}
type DivisionField = Fp<DefaultBackend<DivisionConfig>>;

// Two arbitrary, distinct, nonzero, full-width values -- not reduced mod p
// (add/mul don't require canonical input to execute correctly-shaped
// arithmetic, and the full 256-bit width matches a real field element's
// typical magnitude). Same convention as
// elliptic-curve/benches/elliptic_curve_ops.rs and field/benches/field_ops.rs.
const A: [u64; 4] = [0x1111_1111_2222_2222, 0x3333_3333_4444_4444, 0x5555_5555_6666_6666, 0x0777_7777_8888_8888];
const B: [u64; 4] = [0xaaaa_aaaa_bbbb_bbbb, 0xcccc_cccc_dddd_dddd, 0xeeee_eeee_1f1f_1f1f, 0x0101_0101_0202_0202];

// A second point's coordinates, distinct from (A, B) -- so "add" hits the
// chord branch and "double" (self + self on (A, B)) hits the tangent
// branch rather than the vertical-chord exceptional case.
const X2: [u64; 4] = [0x2222_2222_3333_3333, 0x4444_4444_5555_5555, 0x6666_6666_7777_7777, 0x0888_8888_9999_9999];
const Y2: [u64; 4] = [0xbbbb_bbbb_cccc_cccc, 0xdddd_dddd_eeee_eeee, 0x1f1f_1f1f_2020_2020, 0x0202_0202_0303_0303];

// Full-width scalar so scalar_mul's double-and-add walks close to all 256
// bits, matching field_ops.rs/elliptic_curve_ops.rs's own full-width scalars.
const SCALAR: [u64; 4] = [0xfedc_ba98_7654_3210, 0x0123_4567_89ab_cdef, 0x1357_9bdf_2468_ace0, 0x0f0e_0d0c_0b0a_0908];

fn bench_field_fold(c: &mut Criterion) {
    let fa = Secp256k1Field::new(U256::from_limbs(A));
    let fb = Secp256k1Field::new(U256::from_limbs(B));

    let mut group = c.benchmark_group("Field (fold reduction)");
    group.bench_function("add", |b| b.iter(|| black_box(fa) + black_box(fb)));
    group.bench_function("sub", |b| b.iter(|| black_box(fa) - black_box(fb)));
    group.bench_function("mul", |b| b.iter(|| black_box(fa) * black_box(fb)));
    group.bench_function("square", |b| b.iter(|| black_box(fa).square()));
    group.bench_function("neg", |b| b.iter(|| -black_box(fa)));
    group.bench_function("inverse", |b| b.iter(|| black_box(fa).inverse()));
    group.finish();
}

fn bench_field_division(c: &mut Criterion) {
    let fa = DivisionField::new(U256::from_limbs(A));
    let fb = DivisionField::new(U256::from_limbs(B));

    let mut group = c.benchmark_group("Field (division reduction)");
    group.bench_function("add", |b| b.iter(|| black_box(fa) + black_box(fb)));
    group.bench_function("sub", |b| b.iter(|| black_box(fa) - black_box(fb)));
    group.bench_function("mul", |b| b.iter(|| black_box(fa) * black_box(fb)));
    group.bench_function("square", |b| b.iter(|| black_box(fa).square()));
    group.bench_function("neg", |b| b.iter(|| -black_box(fa)));
    group.bench_function("inverse", |b| b.iter(|| black_box(fa).inverse()));
    group.finish();
}

// Curve constants (A=0, B=7) and point coordinates below use Secp256k1's
// real curve, but arbitrary (not validated on-curve) points -- add/double/
// scalar_mul never call validate(), so no genuine point is needed to
// benchmark the arithmetic (same convention as elliptic_curve_ops.rs).
fn bench_curve_affine(c: &mut Criterion) {
    type Point = AffinePoint<Secp256k1>;

    let p1 = Point { x: Secp256k1Field::new(U256::from_limbs(A)), y: Secp256k1Field::new(U256::from_limbs(B)), infinity: false };
    let p2 = Point { x: Secp256k1Field::new(U256::from_limbs(X2)), y: Secp256k1Field::new(U256::from_limbs(Y2)), infinity: false };
    let scalar = U256::from_limbs(SCALAR);

    let mut group = c.benchmark_group("ShortWeierstrass (Affine)");
    group.bench_function("add", |b| b.iter(|| black_box(p1) + black_box(p2)));
    group.bench_function("double", |b| b.iter(|| black_box(p1) + black_box(p1)));
    group.bench_function("scalar_mul", |b| b.iter(|| black_box(p1) * black_box(scalar)));
    group.finish();
}

fn bench_curve_jacobian(c: &mut Criterion) {
    type Point = JacobianPoint<Secp256k1>;

    let p1 = Point::from_affine(AffinePoint {
        x: Secp256k1Field::new(U256::from_limbs(A)),
        y: Secp256k1Field::new(U256::from_limbs(B)),
        infinity: false,
    });
    let p2 = Point::from_affine(AffinePoint {
        x: Secp256k1Field::new(U256::from_limbs(X2)),
        y: Secp256k1Field::new(U256::from_limbs(Y2)),
        infinity: false,
    });
    let scalar = U256::from_limbs(SCALAR);

    let mut group = c.benchmark_group("ShortWeierstrass (Jacobian)");
    group.bench_function("add", |b| b.iter(|| black_box(p1) + black_box(p2)));
    group.bench_function("double", |b| b.iter(|| black_box(p1).double()));
    group.bench_function("scalar_mul", |b| b.iter(|| black_box(p1) * black_box(scalar)));
    group.finish();
}

criterion_group!(benches, bench_field_fold, bench_field_division, bench_curve_affine, bench_curve_jacobian);
criterion_main!(benches);
