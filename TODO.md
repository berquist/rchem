# quantum chemistry calculation workflow

Requirements:
- start with molecular coordinates
- calculate quantum mechanical wavefunction


1. Input molecular coordinates from file
2. Gather basis set from file and place on atoms
3. Calculate overlap (S), kinetic energy (T), nuclear-electron attraction (V) integrals
4. Calculate initial guess by diagonalizing H = T + V
5. Do SCF iterations: integral direct TEIs, solve HF equations...
6. ...
7. Profit!

# Plan post-2019-04-20

1. Find crate containing gamma function: `rgsl`?
2. Find crate containing n-dimensional array crate w/ eigendecomposition
3. What format should input files be in (Conf, TOML, JSON? not YAML)
4. Read in JSON-formatted basis set definitions from EMSL
5. How to store basis sets? Can we wrap https://github.com/MolSSI-BSE/basis_set_exchange (Python)?
