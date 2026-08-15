# secp256k1

The secp256k1 short Weierstrass curve (`y² = x³ + 7`), built on `elliptic-curve`'s generic `AffinePoint`/`JacobianPoint` and `field-core`'s generic `Fp<B>`.

## What's here

* `Secp256k1`: the curve itself (`ShortWeierstrassCurve` impl, `A = 0`, `B = 7`).
* `Secp256k1Field`: the base field (`Fp<FieldBackend>`), modulus `p = 2^256 - 2^32 - 977`.
* `Secp256k1Scalar`: the scalar field (`Fp<DefaultBackend<ScalarConfig>>`), modulus `n` (the group order).
* `Secp256k1Point`: `JacobianPoint<Secp256k1>` -- the crate's primary point representation. `add`/`double`/`scalar_mul` are all inversion-free here, unlike affine (see Benchmark Results below), so this is the type to reach for unless you specifically need affine's canonical `(x, y)` form.
* `Secp256k1AffinePoint`: `AffinePoint<Secp256k1>`, for interop / serialization / equality checks -- Jacobian's `(X, Y, Z)` isn't unique per point, so converting to affine (`Secp256k1Point::to_affine`) is how you get a canonical, comparable representation back.
* `SECP256K1_P`: the base field modulus, exposed directly.
* `G`: the generator point (`Secp256k1AffinePoint`), the base point of the order-`n` subgroup. Convert via `Secp256k1Point::from_affine(G)` for scalar multiplication.

No serialization or public-key derivation helpers yet — see the root [README](../../README.md)'s "Named Curves" checklist for what's still open.

## The base field's fast reduction

`p = 2^256 - 2^32 - 977` is a Solinas-style prime: `2^256 ≡ C (mod p)` for the small constant `C = 2^32 + 977`. `FieldConfig::mul` (`src/lib.rs`) exploits that directly instead of doing a general `U512 % U256` division: a 512-bit product `x = x_hi·2^256 + x_lo` folds to `x_hi·C + x_lo (mod p)`, and since `C` is only 33 bits, a handful of folds shrinks any product down to canonical range.

The rest of `FieldBackend` isn't hand-written: `FieldConfig` implements `field_core::WideFieldConfig`, supplying only `MODULUS` and `mul`. `add`/`sub`/`neg`/`one`, and the default `square` (`mul(a, a)`), are all derived by `field_core::WideEuclideanBackend` via `U256`'s `WideInt` impl (`bigint`) — which also makes `add`/`sub` correct by construction even though `p` sits close to `U256`'s full 256-bit width (a same-width add can silently overflow there; `WideEuclideanBackend` routes through the widened `U512` type instead). `WideEuclideanBackend` differs from the sibling `WideArithmeticBackend` in exactly one thing: it keeps `FpBackend`'s default GCD-based `inverse` instead of overriding it with Fermat's little theorem, since GCD is far cheaper than Fermat for a modulus this wide (see `gcd_inverse`'s doc comment in `field-core`).

The scalar field (`n`, the group order) has no such power-of-two-minus-small-constant shape, so `Secp256k1Scalar` just uses the generic division-based `DefaultBackend`.

## Benchmark Results

Run via `cargo bench -p secp256k1` (`benches/secp256k1_ops.rs`). Median timings from `criterion`; local-machine numbers, not portable across hardware — useful for relative comparison, not absolute performance claims.

**Base field: fold reduction vs. plain division**

Same modulus (`p`) on both sides — `Field (fold reduction)` is `Secp256k1Field` (`FieldConfig`'s hand-written `mul`); `Field (division reduction)` is the same modulus wired through the generic `DefaultBackend` (`U512 % U256`), to isolate exactly what the fast reduction buys.

| op | Field (fold reduction) | Field (division reduction) |
|---|---|---|
| add | 2.91 ns | 1.62 ns |
| sub | 2.94 ns | 1.85 ns |
| mul | 22.3 ns | 655 ns |
| square | 21.9 ns | 621 ns |
| neg | 1.66 ns | 1.62 ns |
| inverse | 896 ns | 905 ns |

`mul` is ~29x faster with the fold-based reduction — the entire point of hand-writing it. `add`/`sub` are slightly *slower* than the division backend's, not faster: `WideEuclideanBackend` routes them through the widened `U512` type to stay correct when `p` is close to `U256`'s full width (see above), which costs a bit more than a same-width add on 4 limbs instead of 8. `inverse` is essentially identical, since neither backend overrides it — both fall through to the same GCD-based default. `neg` is unaffected either way (it's just `MODULUS - a`, no widening needed since `a <= MODULUS` can't overflow).

**Curve point operations**

`Secp256k1`'s real `A = 0`, `B = 7`, with arbitrary (not validated on-curve) points — `add`/`double`/`scalar_mul` never call `validate()`, so no genuine point is needed to benchmark the arithmetic. `scalar_mul` uses a full-width 256-bit scalar.

| op | Affine | Jacobian |
|---|---|---|
| add | 558 ns | 257 ns |
| double | 1.51 µs | 218 ns |
| scalar_mul | 867 µs | 119 µs |

Jacobian coordinates pay off exactly as the inversion-free design predicts: `add` is ~2.2x faster and `double` is ~6.9x faster than affine, since affine's chord-and-tangent formulas need one field inversion per call (GCD-based `inverse`'s cost varies with the operand's bit pattern -- see the field table above -- but at ~900 ns for an arbitrary value it's on the same order as affine `add`'s entire ~558 ns, so it's a substantial fraction of the cost either way) while Jacobian defers that to a single inversion at `to_affine()`, which this benchmark doesn't even call. `scalar_mul` (double-and-add over a 256-bit scalar) inherits that gap almost directly: ~7.3x faster, since it's dominated by ~256 doublings plus up to 256 adds.
