// Exactly one scalar-mul-* feature selects, at compile time, which
// algorithm every point type's Mul<R> impl uses -- see each point type's
// three cfg-gated `impl Mul<R>` blocks (e.g. short_weierstrass.rs) and
// ladder.rs for the shared constant-time/variable-time helpers.
#[cfg(not(any(
    feature = "scalar-mul-double-and-add",
    feature = "scalar-mul-ladder-variable",
    feature = "scalar-mul-ladder-constant"
)))]
compile_error!(
    "exactly one of the `scalar-mul-double-and-add`, `scalar-mul-ladder-variable`, or \
     `scalar-mul-ladder-constant` features must be enabled -- none is"
);

#[cfg(any(
    all(feature = "scalar-mul-double-and-add", feature = "scalar-mul-ladder-variable"),
    all(feature = "scalar-mul-double-and-add", feature = "scalar-mul-ladder-constant"),
    all(feature = "scalar-mul-ladder-variable", feature = "scalar-mul-ladder-constant")
))]
compile_error!(
    "the `scalar-mul-double-and-add`, `scalar-mul-ladder-variable`, and \
     `scalar-mul-ladder-constant` features are mutually exclusive -- more than one is enabled"
);

mod ladder;

mod short_weierstrass;
pub use short_weierstrass::{AffinePoint, JacobianPoint, ShortWeierstrassCurve};

mod twisted_edwards;
pub use twisted_edwards::{AffinePoint as EdwardsAffinePoint, ExtendedPoint, TwistedEdwardsCurve};

mod montgomery;
pub use montgomery::{AffinePoint as MontgomeryAffinePoint, MontgomeryCurve};
