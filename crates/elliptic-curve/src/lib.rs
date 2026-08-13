// Exactly one scalar-mul-* feature selects, at compile time, which
// algorithm every point type's Mul<R> impl uses -- see each point type's
// two cfg-gated `impl Mul<R>` blocks (e.g. short_weierstrass.rs) and
// ladder.rs for the shared constant-time cswap/select helpers.
#[cfg(not(any(feature = "scalar-mul-variable-time", feature = "scalar-mul-constant-time")))]
compile_error!(
    "exactly one of the `scalar-mul-variable-time` or `scalar-mul-constant-time` features must \
     be enabled -- neither is"
);

#[cfg(all(feature = "scalar-mul-variable-time", feature = "scalar-mul-constant-time"))]
compile_error!(
    "the `scalar-mul-variable-time` and `scalar-mul-constant-time` features are mutually \
     exclusive -- both are enabled"
);

#[cfg(feature = "scalar-mul-constant-time")]
mod ladder;

mod short_weierstrass;
pub use short_weierstrass::{AffinePoint, JacobianPoint, ShortWeierstrassCurve};

mod twisted_edwards;
pub use twisted_edwards::{AffinePoint as EdwardsAffinePoint, ExtendedPoint, TwistedEdwardsCurve};

mod montgomery;
pub use montgomery::{AffinePoint as MontgomeryAffinePoint, MontgomeryCurve};
