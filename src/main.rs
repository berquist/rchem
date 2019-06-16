#![allow(non_snake_case, unused_variables)]

use std::f64;

use chemfiles;
use ndarray::{Array, Axis, Ix1, Ix2, Slice};
use ndarray_linalg::*;

mod basis;
mod integrals;
mod shell;

fn main() {
    // http://www.patorjk.com/software/taag/#p=display&f=3D%20Diagonal&t=rchem
    let logo = r#"
                   ,---,                       ____
                 ,--.' |                     ,'  , `.
  __  ,-.        |  |  :                  ,-+-,.' _ |
,' ,'/ /|        :  :  :               ,-+-. ;   , ||
'  | |' | ,---.  :  |  |,--.   ,---.  ,--.'|'   |  ||
|  |   ,'/     \ |  :  '   |  /     \|   |  ,', |  |,
'  :  / /    / ' |  |   /' : /    /  |   | /  | |--'
|  | ' .    ' /  '  :  | | |.    ' / |   : |  | ,
;  : | '   ; :__ |  |  ' | :'   ;   /|   : |  |/
|  , ; '   | '.'||  :  :_:,''   |  / |   | |`-'
 ---'  |   :    :|  | ,'    |   :    |   ;/
        \   \  / `--''       \   \  /'---'
         `----'               `----'
"#;
    println!("{}", logo);

    let mut trajectory = chemfiles::Trajectory::open("water_crawford.xyz", 'r').unwrap();
    let mut frame = chemfiles::Frame::new();
    trajectory.read(&mut frame).unwrap();
    let natom = frame.size();
    let atomcoords = frame.positions();
    let atomnos: Vec<_> = (0..natom).map(|i| frame.atom(i).atomic_number()).collect();
    let basis_set = basis::Basis::new(&atomnos, &atomcoords, "STO-3G");
    println!("{:#?}", basis_set);

    let S = basis::S(&basis_set);

    let (overlap_eigvals, overlap_eigvecs) = S.eigh(UPLO::Upper).unwrap();
    let overlap_eigvals_inv_sqrt = vec_to_diag_mat(&overlap_eigvals)
        .inv()
        .unwrap()
        .mapv(f64::sqrt);
    let symm_orthog = overlap_eigvecs
        .dot(&overlap_eigvals_inv_sqrt)
        .dot(&overlap_eigvecs.t());

    let T = basis::T(&basis_set);
    let V = basis::V(&basis_set, &atomcoords, &atomnos);
    let H = T + V;
    let F_prime = symm_orthog.t().dot(&H).dot(&symm_orthog);
    let (eps_vec, C_prime) = F_prime.eigh(UPLO::Upper).unwrap();
    let C = symm_orthog.dot(&C_prime);
    let D = build_density(&C, 5);

    let e_elec_new = calc_elec_energy(&D, &H, &H);
    println!("{}", e_elec_new);

    let thresh_e = 1.0e-15;
    let thresh_d = 1.0e-10;
    let max_iterations: u64 = 1024;
    let mut iteration = 0;

    while iteration < max_iterations {
        iteration += 1;
    }
}

fn vec_to_diag_mat(vec: &Array<f64, Ix1>) -> Array<f64, Ix2> {
    let dim = vec.shape()[0];
    let mut mat = Array::zeros((dim, dim));
    for i in 0..dim {
        mat[[i, i]] = vec[i];
    }
    mat
}

fn build_density(C: &Array<f64, Ix2>, nocc: usize) -> Array<f64, Ix2> {
    C.slice_axis(Axis(1), Slice::from(..nocc))
        .dot(&C.slice_axis(Axis(1), Slice::from(..nocc)).t())
}

fn calc_elec_energy(D: &Array<f64, Ix2>, H: &Array<f64, Ix2>, F: &Array<f64, Ix2>) -> f64 {
    ((H + F) * D).sum()
}
