#![allow(non_snake_case, unused_variables)]

use std::f64;

use ndarray::{Array, Axis, Ix1, Ix2, Slice};
use ndarray_linalg::*;

use rchem::basis;

/// The main driver. Currently just an example of computing the HF/STO-2G energy of a water molecule.
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
    let basis_set = basis::Basis::new(&atomnos, &atomcoords, "STO-2G");
    // println!("{:#?}", basis_set);

    // let I = basis::build_I(&basis_set);

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
    let mut D = build_density(&C, 5);

    let mut e_elec_new = calc_elec_energy(&D, &H, &H);

    let thresh_e = 1.0e-11;
    // let thresh_d = 1.0e-10;
    let max_iterations = 1024;
    let mut iteration = 0;

    while iteration < max_iterations {
        let F = build_fock(&D, &H, &basis_set);
        let F_prime = symm_orthog.t().dot(&F).dot(&symm_orthog);
        let (eps_vec, C_prime) = F_prime.eigh(UPLO::Upper).unwrap();
        let C = symm_orthog.dot(&C_prime);
        // let D_old = D.clone();
        D = build_density(&C, 5);
        let e_elec_old = e_elec_new;
        e_elec_new = calc_elec_energy(&D, &H, &F);
        let e_total = e_elec_new;
        let delta_e = e_elec_new - e_elec_old;
        println!("{:4} {:20.12} {:20.12}", iteration, e_elec_new, delta_e);
        if delta_e.abs() < thresh_e {
            println!("Convergence achieved!");
            break;
        }
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

fn build_fock(
    D: &Array<f64, Ix2>,
    H: &Array<f64, Ix2>,
    basis_set: &basis::Basis,
) -> Array<f64, Ix2> {
    let (J, K) = basis::JK_direct(&basis_set, &D);
    (2.0 * J - K) + H
}

// fn build_fock_inmem(
//     D: &Array<f64, Ix2>,
//     H: &Array<f64, Ix2>,
//     I: &Array<f64, Ix4>,
// ) -> Array<f64, Ix2> {
//     let (J, K) = basis::JK_inmem(&I, &D);
//     (2.0 * J - K) + H
// }
