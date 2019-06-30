#![allow(non_snake_case)]

use std::collections::HashMap as Map;
use std::collections::HashSet as Set;
use std::f64::consts::PI;

use cpython::{PyDict, Python};
use ndarray::{Array, Axis, Ix2, Ix4};
use serde::{Deserialize, Deserializer};
use serde_json;

use crate::integrals;
use crate::shell;

#[derive(Debug, Deserialize)]
struct BSEResult {
    name: String,
    description: String,
    // TODO how to make this optional?
    // notes: String,
    elements: Map<u8, BSEElement>,
}

#[derive(Debug, Deserialize)]
struct BSEElement {
    electron_shells: Vec<BSEElectronShell>,
}

#[derive(Debug, Deserialize)]
enum BSEFunctionType {
    #[serde(rename = "gto")]
    GTO,
    #[serde(rename = "gto_spherical")]
    GTOSpherical,
    #[serde(rename = "gto_cartesian")]
    GTOCartesian,
    #[serde(rename = "sto")]
    STO,
}

#[derive(Debug, Deserialize)]
enum BSERegion {
    #[serde(rename = "")]
    All,
    #[serde(rename = "valence")]
    Valence,
    #[serde(rename = "polarization")]
    Polarization,
    #[serde(rename = "core")]
    Core,
    #[serde(rename = "tight")]
    Tight,
    #[serde(rename = "diffuse")]
    Diffuse,
}

fn vec_strings_to_f64(v: &Vec<String>) -> Vec<f64> {
    v.iter().map(|x| x.parse::<f64>().unwrap()).collect()
}

fn deserialize_vec_vec_string_to_vec_vec_f64<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<Vec<f64>>, D::Error> {
    let v: Vec<Vec<String>> = Deserialize::deserialize(deserializer)?;
    Ok(v.iter().map(|x| vec_strings_to_f64(&x)).collect())
}

fn deserialize_vec_string_to_vec_f64<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<Vec<f64>, D::Error> {
    let v: Vec<String> = Deserialize::deserialize(deserializer)?;
    Ok(vec_strings_to_f64(&v))
}

#[derive(Debug, Deserialize)]
struct BSEElectronShell {
    angular_momentum: Vec<usize>,
    // both the coefficients and exponents are stored as strings in the Basis
    // Set Exchange in order to maintain scientific notation
    #[serde(deserialize_with = "deserialize_vec_vec_string_to_vec_vec_f64")]
    coefficients: Vec<Vec<f64>>,
    #[serde(deserialize_with = "deserialize_vec_string_to_vec_f64")]
    exponents: Vec<f64>,
    function_type: BSEFunctionType,
    region: BSERegion,
}

fn get_bse_json(py: Python, basis_set_name: &str, elements: &Vec<u64>) -> BSEResult {
    let locals = PyDict::new(py);
    locals
        .set_item(py, "bse", py.import("basis_set_exchange").unwrap())
        .unwrap();
    // Copy over the elements
    let unique_elements: Set<u64> = elements.into_iter().map(|&x| x).collect();
    let unique_elements: Vec<String> = unique_elements.iter().map(|x| x.to_string()).collect();
    let unique_elements = unique_elements.join(", ");
    let call = format!(
        "bse.get_basis(\"{}\", elements=[{}], fmt=\"json\")",
        basis_set_name, unique_elements
    );
    let jsonstr: String = py
        .eval(&call, None, Some(&locals))
        .unwrap()
        .extract(py)
        .unwrap();
    serde_json::from_str(&jsonstr).unwrap()
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
    powers: [usize; 3],
    exponent: f64,
    norm: f64,
}

impl PGTO {
    fn new(origin: [f64; 3], powers: [usize; 3], exponent: f64) -> PGTO {
        let mut ret = PGTO {
            origin: origin,
            powers: powers,
            exponent: exponent,
            norm: 0.0,
        };
        ret.norm = ret.normalization();
        ret
    }

    fn order(&self) -> usize {
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

#[derive(Debug)]
struct CGTO {
    origin: [f64; 3],
    powers: [usize; 3],
    coefs: Vec<f64>,
    primitives: Vec<PGTO>,
}

impl CGTO {
    fn from_pgtos(pgtos: &Vec<PGTO>, coefs: &Vec<f64>) -> CGTO {
        // TODO Ok and Err?
        assert!(pgtos.len() > 0);
        CGTO {
            origin: pgtos[0].origin,
            powers: pgtos[0].powers,
            // exponents: pgtos.iter().map(|pgto| pgto.exponent).collect(),
            // norms: pgtos.iter().map(|pgto| pgto.norm).collect(),
            coefs: coefs.clone(),
            primitives: pgtos.clone(),
        }
    }
}

#[derive(Debug)]
pub struct Basis {
    name: String,
    cgtos: Vec<CGTO>,
}

impl Basis {
    pub fn new(atomnos: &Vec<u64>, all_atomcoords: &[[f64; 3]], basis_set_name: &str) -> Basis {
        let mut cgtos: Vec<CGTO> = Vec::new();
        let gil = Python::acquire_gil();
        let bseresult = get_bse_json(gil.python(), basis_set_name, &atomnos);
        for (i, &atomno) in atomnos.iter().enumerate() {
            let atomcoords = all_atomcoords[i];
            let element = &bseresult.elements[&(atomno as u8)];
            for shell in &element.electron_shells {
                for angular_momentum in &shell.angular_momentum {
                    assert_eq!(
                        shell.exponents.len(),
                        shell.coefficients[*angular_momentum].len()
                    );
                    for powers in shell::get_ijk_list(*angular_momentum) {
                        let pgtos: Vec<_> = shell
                            .exponents
                            .iter()
                            .map(|exponent| PGTO::new(atomcoords.clone(), powers, *exponent))
                            .collect();
                        let cgto = CGTO::from_pgtos(&pgtos, &shell.coefficients[*angular_momentum]);
                        cgtos.push(cgto);
                    }
                }
            }
        }
        Basis {
            name: basis_set_name.to_string(),
            cgtos: cgtos,
        }
    }
}

fn overlap_pgto(a: &PGTO, b: &PGTO) -> f64 {
    let powers = [
        a.powers[0],
        a.powers[1],
        a.powers[2],
        b.powers[0],
        b.powers[1],
        b.powers[2],
    ];
    a.norm * b.norm * integrals::get_overlap(a.exponent, b.exponent, &a.origin, &b.origin, &powers)
}

fn overlap_cgto_left(a: &CGTO, b: &PGTO) -> f64 {
    a.primitives
        .iter()
        .zip(&a.coefs)
        .map(|(pa, ca)| ca * overlap_pgto(&pa, &b))
        .sum()
}

pub fn S(basis_set: &Basis) -> Array<f64, Ix2> {
    let dim = basis_set.cgtos.len();
    let mut mat: Array<f64, _> = Array::zeros((dim, dim));
    for mu in 0..dim {
        let a = &basis_set.cgtos[mu];
        for nu in 0..mu + 1 {
            let b = &basis_set.cgtos[nu];
            mat[[mu, nu]] = b
                .primitives
                .iter()
                .zip(&b.coefs)
                .map(|(pb, cb)| cb * overlap_cgto_left(&a, &pb))
                .sum();
            mat[[nu, mu]] = mat[[mu, nu]];
        }
    }
    mat
}

fn kinetic_pgto(a: &PGTO, b: &PGTO) -> f64 {
    let powers = [
        a.powers[0],
        a.powers[1],
        a.powers[2],
        b.powers[0],
        b.powers[1],
        b.powers[2],
    ];
    a.norm * b.norm * integrals::get_kinetic(a.exponent, b.exponent, &a.origin, &b.origin, &powers)
}

fn kinetic_cgto_left(a: &CGTO, b: &PGTO) -> f64 {
    a.primitives
        .iter()
        .zip(&a.coefs)
        .map(|(pa, ca)| ca * kinetic_pgto(&pa, &b))
        .sum()
}

pub fn T(basis_set: &Basis) -> Array<f64, Ix2> {
    let dim = basis_set.cgtos.len();
    let mut mat: Array<f64, _> = Array::zeros((dim, dim));
    for mu in 0..dim {
        let a = &basis_set.cgtos[mu];
        for nu in 0..mu + 1 {
            let b = &basis_set.cgtos[nu];
            mat[[mu, nu]] = b
                .primitives
                .iter()
                .zip(&b.coefs)
                .map(|(pb, cb)| cb * kinetic_cgto_left(&a, &pb))
                .sum();
            mat[[nu, mu]] = mat[[mu, nu]];
        }
    }
    mat
}

fn nuclear_pgto(a: &PGTO, b: &PGTO, atomcoords: &[f64; 3]) -> f64 {
    let powers = [
        a.powers[0],
        a.powers[1],
        a.powers[2],
        b.powers[0],
        b.powers[1],
        b.powers[2],
    ];
    a.norm
        * b.norm
        * integrals::get_nuclear(
            a.exponent, b.exponent, &a.origin, &b.origin, atomcoords, &powers,
        )
}

fn nuclear_cgto_left(a: &CGTO, b: &PGTO, atomcoords: &[f64; 3]) -> f64 {
    a.primitives
        .iter()
        .zip(&a.coefs)
        .map(|(pa, ca)| ca * nuclear_pgto(&pa, &b, atomcoords))
        .sum()
}

pub fn V(basis_set: &Basis, atomcoords: &[[f64; 3]], atomnos: &Vec<u64>) -> Array<f64, Ix2> {
    let dim = basis_set.cgtos.len();
    let natoms = atomcoords.len();
    let mut mat: Array<f64, _> = Array::zeros((dim, dim, natoms));
    for (c, single_atomcoords) in atomcoords.iter().enumerate() {
        let atomno = atomnos[c] as f64;
        for mu in 0..dim {
            let a = &basis_set.cgtos[mu];
            for nu in 0..mu + 1 {
                let b = &basis_set.cgtos[nu];
                let elem: f64 = b
                    .primitives
                    .iter()
                    .zip(&b.coefs)
                    .map(|(pb, cb)| cb * nuclear_cgto_left(&a, &pb, single_atomcoords))
                    .sum();
                mat[[mu, nu, c]] = atomno * elem;
                mat[[nu, mu, c]] = mat[[mu, nu, c]];
            }
        }
    }
    mat.sum_axis(Axis(2))
}

fn coulomb_pgto(a: &PGTO, b: &PGTO, c: &PGTO, d: &PGTO) -> f64 {
    let powers = [
        a.powers[0],
        a.powers[1],
        a.powers[2],
        b.powers[0],
        b.powers[1],
        b.powers[2],
        c.powers[0],
        c.powers[1],
        c.powers[2],
        d.powers[0],
        d.powers[1],
        d.powers[2],
    ];
    a.norm
        * b.norm
        * c.norm
        * d.norm
        * integrals::get_coulomb(
            a.exponent, b.exponent, c.exponent, d.exponent, &a.origin, &b.origin, &c.origin,
            &d.origin, &powers,
        )
}

pub fn JK_direct(basis_set: &Basis, D: &Array<f64, Ix2>) -> (Array<f64, Ix2>, Array<f64, Ix2>) {
    let dim = basis_set.cgtos.len();
    let mut J: Array<f64, _> = Array::zeros((dim, dim));
    let mut K: Array<f64, _> = Array::zeros((dim, dim));
    for mu in 0..dim {
        let a = &basis_set.cgtos[mu];
        for nu in 0..dim {
            let b = &basis_set.cgtos[nu];
            for lambda in 0..dim {
                let c = &basis_set.cgtos[lambda];
                for sigma in 0..dim {
                    let d = &basis_set.cgtos[sigma];
                    let mut j_contr = 0.0;
                    let mut k_contr = 0.0;
                    for (pa, ca) in a.primitives.iter().zip(&a.coefs) {
                        for (pb, cb) in b.primitives.iter().zip(&b.coefs) {
                            for (pc, cc) in c.primitives.iter().zip(&c.coefs) {
                                for (pd, cd) in d.primitives.iter().zip(&d.coefs) {
                                    j_contr += ca
                                        * cb
                                        * cc
                                        * cd
                                        * coulomb_pgto(&pa, &pb, &pc, &pd)
                                        * D[[lambda, sigma]];
                                    k_contr += ca
                                        * cb
                                        * cc
                                        * cd
                                        * coulomb_pgto(&pa, &pc, &pb, &pd)
                                        * D[[lambda, sigma]];
                                }
                            }
                        }
                    }
                    J[[mu, nu]] += j_contr;
                    K[[mu, nu]] += k_contr;
                }
            }
        }
    }
    (J, K)
}

pub fn build_I(basis_set: &Basis) -> Array<f64, Ix4> {
    let dim = basis_set.cgtos.len();
    let mut I: Array<f64, _> = Array::zeros((dim, dim, dim, dim));
    for mu in 0..dim {
        let a = &basis_set.cgtos[mu];
        for nu in 0..mu + 1 {
            let b = &basis_set.cgtos[nu];
            for lambda in 0..dim {
                let c = &basis_set.cgtos[lambda];
                for sigma in 0..lambda + 1 {
                    let d = &basis_set.cgtos[sigma];
                    let mut val = 0.0;
                    for (pa, ca) in a.primitives.iter().zip(&a.coefs) {
                        for (pb, cb) in b.primitives.iter().zip(&b.coefs) {
                            for (pc, cc) in c.primitives.iter().zip(&c.coefs) {
                                for (pd, cd) in d.primitives.iter().zip(&d.coefs) {
                                    val +=
                                        ca * cb * cc * cd * coulomb_pgto(&pa, &pb, &pc, &pd);
                                }
                            }
                        }
                    }
                    I[[mu, nu, lambda, sigma]] = val;
                    I[[nu, mu, lambda, sigma]] = val;
                    I[[mu, nu, sigma, lambda]] = val;
                    I[[nu, mu, sigma, lambda]] = val;
                }
            }
        }
    }
    I
}

pub fn JK_inmem(I: &Array<f64, Ix4>, D: &Array<f64, Ix2>) -> (Array<f64, Ix2>, Array<f64, Ix2>) {
    let dim = I.shape()[0];
    let mut J: Array<f64, _> = Array::zeros((dim, dim));
    let mut K: Array<f64, _> = Array::zeros((dim, dim));
    for mu in 0..dim {
        for nu in 0..mu + 1 {
            let mut j_contr = 0.0;
            let mut k_contr = 0.0;
            for lambda in 0..dim {
                for sigma in 0..dim {
                    j_contr += I[[mu, nu, lambda, sigma]] * D[[lambda, sigma]];
                    k_contr += I[[mu, lambda, nu, sigma]] * D[[lambda, sigma]];
                }
            }
            J[[mu, nu]] = j_contr;
            K[[mu, nu]] = k_contr;
            // TODO can't do J += J.t()?
            J[[nu, mu]] = j_contr;
            K[[nu, mu]] = k_contr;
        }
    }
    (J, K)
}
