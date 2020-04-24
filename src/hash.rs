// Hash coefficients
const A: u64 = 11_927_359_292_700_924_260;
const B: u64 = 6_512_515_406_574_399_413;

/// See https://arxiv.org/pdf/1303.5479v2.pdf for details. This is a 2-independent hash
/// function, aiming for speed.
pub fn hash(value: u32) -> u32 {
    (((A.overflowing_mul(value as u64).0).overflowing_add(B).0) >> 32) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hash() {
        let result = hash(100);
        assert_eq!(result, 48_913_034);
    }
}
