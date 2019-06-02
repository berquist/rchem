use std::f64::consts::PI;

use crate::integrals;

fn fact2(n: isize) -> isize {
    if n <= 0 {
        1
    } else {
        n * fact2(n - 2)
    }
}

#[derive(Clone, Debug)]
pub struct PGTO {
    origin: [f64; 3],
    exponent: f64,
    powers: [u8; 3],
    norm: f64,
}

impl PGTO {
    pub fn new(origin: [f64; 3], exponent: f64, powers: [u8; 3]) -> PGTO {
        let mut ret = PGTO {
            origin,
            exponent,
            powers,
            norm: 0.0,
        };
        ret.norm = ret.normalization();
        ret
    }

    fn order(&self) -> u8 {
        self.powers.iter().sum()
    }

    fn normalization(&self) -> f64 {
        let numer = 2f64.powf(2.0 * self.order() as f64 + 1.5)
            * self.exponent.powf(self.order() as f64 + 1.5);
        let l = self.powers[0] as isize;
        let m = self.powers[1] as isize;
        let n = self.powers[2] as isize;
        let denom: f64 =
            (fact2(2 * l - 1) * fact2(2 * m - 1) * fact2(2 * n - 1)) as f64 * PI.powf(1.5);
        (numer / denom).powf(0.5)
    }
}

// struct CGTO {
//     origin: [f64; 3],
//     exponents: Vec<f64>,
//     powers: [u8; 3],
//     contraction_coefficients: Vec<f64>,
// }

// struct CGTO {
//     primitives: [PGTO],
//     contraction_coefficients: [f64],
// }

pub fn S(a: &PGTO, b: &PGTO) -> f64 {
    let powers = [
        a.powers[0],
        a.powers[1],
        a.powers[2],
        b.powers[0],
        b.powers[1],
        b.powers[2],
    ];
    a.norm * b.norm * integrals::get_overlap(a.exponent, b.exponent, a.origin, b.origin, powers)
}
