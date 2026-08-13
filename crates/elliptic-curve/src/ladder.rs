#[cfg(feature = "scalar-mul-ladder-constant")]
use field_core::Field;
#[cfg(feature = "scalar-mul-ladder-variable")]
use field_core::FpRepr;

// Internal helpers shared by the constant-time and variable-time Montgomery
// ladder Mul<R> impls in short_weierstrass.rs, twisted_edwards.rs, and
// montgomery.rs. Not part of the crate's public API (see lib.rs: this
// module isn't re-exported) -- each point type still writes out its own
// `cswap` over its own fields (x/y, x/y/z, x/y/z/t, ...), since there's no
// existing shared "point" trait to hang a single generic cswap off of and
// introducing one just for this would be more machinery than the three
// call sites below warrant.

// Selects `a` when `bit` is false, `b` when `bit` is true, via field
// arithmetic (`a + bit*(b-a)`) rather than branching on which operand to
// keep. The bit -> {0, 1} conversion below still goes through an `if` since
// Field has no branch-free bit-to-element primitive -- true hardware
// constant-time selection needs that primitive (see the README's "Constant-
// Time Verification (subtle, branch-free select)" item), so this only gets
// the *shape* of constant-time right: every iteration of the ladder that
// calls this does the same field operations regardless of the bit, instead
// of skipping an addition the way double-and-add does.
#[cfg(feature = "scalar-mul-ladder-constant")]
pub(crate) fn select<F: Field>(bit: bool, a: F, b: F) -> F {
    let bit = if bit { F::one() } else { F::zero() };
    a + bit * (b - a)
}

// Same idea as `select` above, but for a plain `bool` flag (e.g.
// AffinePoint's `infinity`) that isn't a Field element and so can't go
// through the arithmetic trick. `&`/`|` on bool are eager, not
// short-circuiting like `&&`/`||`, so this stays branch-free the same way
// `select` is.
#[cfg(feature = "scalar-mul-ladder-constant")]
pub(crate) fn select_bool(bit: bool, a: bool, b: bool) -> bool {
    (!bit & a) | (bit & b)
}

// Position of `e`'s highest set bit, plus one (so 0 for e == 0). Runs in
// time proportional to that position, i.e. proportional to the scalar's own
// magnitude -- which is exactly why only the variable-time ladder calls
// this. The constant-time ladder below iterates a fixed `R::BITS` instead,
// specifically to avoid a scalar-dependent loop bound.
#[cfg(feature = "scalar-mul-ladder-variable")]
pub(crate) fn bit_length<R: FpRepr>(mut e: R) -> usize {
    let mut n = 0;

    while e != R::ZERO {
        n += 1;
        e >>= 1;
    }

    n
}
