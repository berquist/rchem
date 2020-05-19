[![Build Status](https://travis-ci.com/berquist/rchem.svg?branch=master)](https://travis-ci.com/berquist/rchem)
[![License](https://img.shields.io/github/license/berquist/rchem.svg)](LICENSE)
[![Issues](https://img.shields.io/github/issues/berquist/rchem.svg)](https://github.com/berquist/rchem/issues)
<!-- [![crates.io](https://img.shields.io/crates/v/rchem.svg)](https://crates.io/crates/rchem) -->
<!-- [![Docs.rs](https://docs.rs/rchem/badge.svg)](https://docs.rs/rchem) -->

# rchem

_Ab initio_ quantum chemistry in Rust.

## Dependencies

- Rust toolchain (stable)
- CMake and a C compiler (for `pyquante2` molecular integrals)
- Python and [`basis-set-exchange`](https://pypi.org/project/basis-set-exchange/) (for reading basis sets)

## Running

To run the simple restricted Hartree-Fock program,
``` sh
$ cargo run
```
To run benchmarks for the different molecular integral implementations,
``` sh
$ cargo bench
```
To view the library documentation,
``` sh
$ cargo doc
```
