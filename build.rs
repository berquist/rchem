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
        .expect("Unable to generate bindings");
    bindings_libpyquante2
        .write_to_file(out_dir.join("bindings_libpyquante2.rs"))
        .expect("Couldn't write bindings!");
}
