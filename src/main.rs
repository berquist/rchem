use std::collections::HashMap as Map;
use std::collections::HashSet as Set;

use cpython::{PyDict, Python};
use ndarray as nd;
use serde::{Deserialize, Deserializer};
use serde_json;

mod basis;
mod integrals;

fn main() {
    let za = 1.8;
    let zb = 2.8;
    let ra = [0.0, 0.0, 0.0];
    let rb = [0.5, 0.8, -0.2];
    // println!(
    //     "{}",
    //     integrals::get_overlap(za, zb, ra, rb, [2, 1, 0, 1, 1, 0])
    // );

    let coords = nd::arr2(&[
        [0.000000000000, -0.143225816552, 0.000000000000],
        [1.638036840407, 1.136548822547, 0.000000000000],
        [-1.638036840407, 1.136548822547, 0.000000000000],
    ]);
    // println!("{}", coords);

    let a = basis::PGTO::new(
        [coords[[0, 0]], coords[[0, 1]], coords[[0, 2]]],
        1.0,
        [0, 0, 0],
    );
    // println!("{:?}", a);
    // println!("{}", basis::S(&a, &a));

    let gil = Python::acquire_gil();
    let elements = vec![8, 1, 1];
    let bseresult = get_bse_json(gil.python(), "STO-3G", elements);
    println!("{:#?}", bseresult);
}

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

fn get_bse_json(py: Python, basis_set_name: &str, elements: Vec<u8>) -> BSEResult {
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
