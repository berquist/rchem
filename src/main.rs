extern crate ndarray as nd;
use cpython::{PyDict, PyResult, Python};
use serde_json;

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

    let gil = Python::acquire_gil();
    let json = get_bse_json(gil.python());
    let v: serde_json::Value = serde_json::from_str(&json).unwrap();
    println!("{:#?}", v);
}

// TODO how to pass the basis set name?
// TODO how to pass the element numbers?
fn get_bse_json(py: Python) -> String {
    let locals = PyDict::new(py);
    locals
        .set_item(py, "bse", py.import("basis_set_exchange").unwrap())
        .unwrap();
    return py
        .eval(
            "bse.get_basis(\"STO-3G\", elements=[8, 1], fmt=\"json\")",
            None,
            Some(&locals),
        )
        .unwrap()
        .extract(py)
        .unwrap();
}
