use field_core::Field;

// x, y generic over any Field so the same point type works over Fp,
// Fp2, or the small fields (Goldilocks/BabyBear/Mersenne31) alike.
// infinity marks the point at infinity (the group identity); x/y are
// meaningless when it's set, matching the usual affine short-Weierstrass
// convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ECPoint<F: Field> {
    pub x: F,
    pub y: F,
    pub infinity: bool,
}
