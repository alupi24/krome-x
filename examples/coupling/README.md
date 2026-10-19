# Calling KROME-X from your code

The same problem in each language KROME-X can be called from: the cell
of the quick start (the main README), primordial gas at n_H = 1e4 cm^-3
and 100 K, with atomic and H2 cooling and the chemical and compression
heating, followed in time from 1 yr to 1e6 yr in steps 10 % longer each
time.

| folder | language | API |
|---|---|---|
| `c/` | C | the C interface (`krome.h`, `krome_net.h`): `krome_solve_cells` |
| `cpp/` | C++ | the header-only C++ API (`krome.hpp`): `solve` |
| `fortran/` | Fortran | the Fortran binding (modules `krome`, `krome_net`): `krome_solve` |
| `python/` | Python | Jupyter notebooks: `quickstart.ipynb` (this cell), `collapse.ipynb` (a one-zone collapse, with and without the continuum cooling), `networks.ipynb` (documenting and converting networks) |

Each compiled example writes `abundances.dat`: the time [yr], Tgas [K]
and the abundances n_i / n_H of every species.

## Build and run

Generate and build the KROME-X project once, from the KROME-X folder
(`options.toml` holds the options):

    python3 -m kromex examples/coupling/options.toml --set project.outdir=build_coupling
    make -C build_coupling lib fortran shared

Then each example, in its folder (`KROME` = the project folder, default
`../../../build_coupling`):

    cd examples/coupling/c       && make && ./example
    cd examples/coupling/cpp     && make && ./example
    cd examples/coupling/fortran && make && ./example

and compare them:

    python3 examples/coupling/compare.py

The C and Fortran examples call the same library and agree to the
printed digits; the C++ example compiles the solver itself and agrees
to ~1e-6 or better (different floating-point contractions). The
notebooks generate their own projects (`jupyter notebook
examples/coupling/python/`).

## From here to a hydro code

A hydro code makes the same call once per step, over all its cells at
once: `krome_solve_cells` (C), `solve_cells` (C++), `krome_solve` with
`n(nspec, ncell)` (Fortran), `net.solve` with `n[ncell, nspec]`
(Python), with per-cell parameters (here the free-fall time `tff`) as
rows of `par`. The manual (docs/manual.md, sections 5 and 6) covers the
options, the per-cell inputs, the GPU path and the failure flags.
