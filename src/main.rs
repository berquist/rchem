use std::f64::consts::PI;

// extern crate nalgebra as na;
// use na::U3;
extern crate ndarray as nd;

mod integrals;

fn main() {
    let za = 1.8;
    let zb = 2.8;
    let ra = [0.0, 0.0, 0.0];
    let rb = [0.5, 0.8, -0.2];
    println!(
        "{}",
        integrals::get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    );

    // let v = na::Vector3::new(1, 2, 3);
    // let w: na::Vector3<f64> = na::Vector3::new_random();
    // // let x = na::Vector3<f64>::new_random();
    // println!("{}", v);
    // println!("{}", w);
    // // let a: na::MatrixMN<f64, U3, U3> = na::Matrix<_, _, _, _>>::new_random(3, 3);
    // let a = na::DMatrix::<f64>::new_random(3, 3);
    // println!("{}", a);
    // let b = na::Matrix3x3::<f64>::new_random();
    // println!("{}", b);

    // let a = nd::arr2(&[[3.0, 1.0, 1.0],
    //                    [1.0, 3.0, 1.0],
    //                    [1.0, 1.0, 3.0]]);
    // println!("{}", a];
    let coords = nd::arr2(&[
        [0.000000000000, -0.143225816552, 0.000000000000],
        [1.638036840407, 1.136548822547, 0.000000000000],
        [-1.638036840407, 1.136548822547, 0.000000000000],
    ]);
    println!("{}", coords);

    // let a = PGTO::new(
    //     [0.000000000000, -0.143225816552, 0.000000000000],
    //     1.0,
    //     [0, 0, 0],
    // );
    let a = PGTO::new(
        [coords[[0, 0]], coords[[0, 1]], coords[[0, 2]]],
        1.0,
        [0, 0, 0],
    );
    println!("{:?}", a);
    let b = a.clone();
    println!("{:?}", b);
    println!("{}", S(a, b));
}

fn fact2(n: isize) -> isize {
    if n <= 0 {
        1
    } else {
        n * fact2(n - 2)
    }
}

#[derive(Clone, Debug)]
struct PGTO {
    origin: [f64; 3],
    exponent: f64,
    powers: [u8; 3],
    norm: f64,
}

impl PGTO {
    fn new(origin: [f64; 3], exponent: f64, powers: [u8; 3]) -> PGTO {
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

struct CGTO {
    origin: [f64; 3],
    exponents: Vec<f64>,
    powers: [u8; 3],
    contraction_coefficients: Vec<f64>,
}

// struct CGTO {
//     primitives: [PGTO],
//     contraction_coefficients: [f64],
// }

fn S(a: PGTO, b: PGTO) -> f64 {
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
