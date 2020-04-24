// extern crate bindgen;
use cmake;

// use std::env;
// use std::path::PathBuf;

fn main() {
    // println!("cargo:rustc-link-lib=cint");
    // println!("cargo:rerun-if-changed=wrapper.h");
    // let bindings = bindgen::Builder::default()
    //     .header("wrapper.h")
    //     .parse_callbacks(Box::new(bindgen::CargoCallbacks))
    //     .generate()
    //     .expect("Unable to generate bindings");
    // let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    // bindings
    //     .write_to_file(out_path.join("bindings.rs"))
    //     .expect("Couldn't write bindings!");

    let dst = cmake::Config::new("libpyquante2").build();
    println!("cargo:rustc-link-search={}", dst.display());
    println!("cargo:rustc-link-lib=static=pyquante2");
}
