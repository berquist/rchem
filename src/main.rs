// extern crate nalgebra as na;
// use na::U3;
extern crate ndarray as nd;

mod basis;
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

    let a = basis::PGTO::new(
        [coords[[0, 0]], coords[[0, 1]], coords[[0, 2]]],
        1.0,
        [0, 0, 0],
    );
    println!("{:?}", a);
    println!("{}", basis::S(&a, &a));
}
