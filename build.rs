use bindgen;
use cmake;
use std::env;
use std::path::PathBuf;

fn main() {
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());

    let dst = cmake::Config::new("libpyquante2").build();
    println!("cargo:rustc-link-search={}", dst.display());
    println!("cargo:rustc-link-lib=static=pyquante2");
    let bindings_libpyquante2 = bindgen::Builder::default()
        .header("libpyquante2/cints.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate libpyquante2 bindings");
    bindings_libpyquante2
        .write_to_file(out_dir.join("bindings_libpyquante2.rs"))
        .expect("Couldn't write libpyquante2 bindings!");

    let simint_base = "/home/eric/data/opt/apps/simint/0.7-g9.3.0-avxfma";
    println!("cargo:rustc-link-search={}/lib", simint_base);
    println!("cargo:rustc-link-lib=static=simint");
    let bindings_simint = bindgen::Builder::default()
        .header(format!("{}/include/simint/simint.h", simint_base))
        .whitelist_var("SIMINT_.*")
        .whitelist_function("simint_.*")
        .whitelist_type("simint_.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate simint bindings");
    bindings_simint
        .write_to_file(out_dir.join("bindings_simint.rs"))
        .expect("Couldn't write simint bindings!");
}
