# CTSI2025-PBPK-in-Julia

## Instructions

### Installing the packages

- Clone the repository
- Start a Julia REPL and make sure you're in the project directory (the one that has the `Project.toml` and `Manifest.toml` files)
- In the REPL, type 
  ```
  using Pkg
  Pkg.activate(".")
  Pkg.instantiate()
  ```
  This will activate the package environment and install all packages.

### Run Pluto notebooks

- In the REPL, type
  
  ```
  using Pluto
  Pluto.run()
  ```
  This will load the Pluto package and open a server to run the Pluto notebooks.

- Open a notebook (`nb-julia-basics.jl` or `nb-example-voriPBPK.jl`), click on `Run this notebook` and wait for the magic to happen.
