use std::collections::HashMap as Map;
use std::collections::HashSet as Set;
use std::f64::consts::PI;

use cpython::{PyDict, Python};
use serde::{Deserialize, Deserializer};
use serde_json;

use crate::integrals;

#[derive(Debug, Deserialize)]
pub struct BSEResult {
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
    angular_momentum: Vec<u8>,
    // both the coefficients and exponents are stored as strings in the Basis
    // Set Exchange in order to maintain scientific notation
    #[serde(deserialize_with = "deserialize_vec_vec_string_to_vec_vec_f64")]
    coefficients: Vec<Vec<f64>>,
    #[serde(deserialize_with = "deserialize_vec_string_to_vec_f64")]
    exponents: Vec<f64>,
    function_type: BSEFunctionType,
    region: BSERegion,
}

pub fn get_bse_json(py: Python, basis_set_name: &str, elements: Vec<u8>) -> BSEResult {
    let locals = PyDict::new(py);
    locals
        .set_item(py, "bse", py.import("basis_set_exchange").unwrap())
        .unwrap();
    let unique_elements: Set<u8> = elements.into_iter().collect();
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
    let tmp: serde_json::Value = serde_json::from_str(&jsonstr).unwrap();
    // println!("{:#?}", tmp);
    // println!(
    //     "{:#?}",
    //     tmp["elements"]["1"]["electron_shells"][0]["exponents"]
    // );
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
