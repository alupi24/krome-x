# Coupling KROME-X to a hydro code

Four small, complete drivers, one per calling style of the usual host codes.
They all solve **the same problem**, so their outputs agree to the solver
tolerance.


| directory | language / API | modelled on | data layout it starts from |
|---|---|---|---|
| `fortran_ramses/` | Fortran, `krome_solve_x` (mass fractions) | RAMSES `cooling_fine` | `uold(ncell, nvar)`, species as passive scalars rho X, processed in vectors of `nvector` cells |
| `c_enzo/` | C, `krome_solve_cells` | ENZO | one field per quantity (`species[s][cell]` mass densities), code units |
| `cpp_gizmo/` | C++ header-only API, `solve_cells` with per-particle `Params` | GIZMO / GADGET | array of particle structs, per-particle dt |
| `python_post/` | Python `net.solve`, `net.thermal` | post-processing a snapshot | numpy arrays over cells |

## The problem

16 cells with n_H = 1 ... 1e4 cm^-3, 1e4 K, 1 % ionized primordial gas
(network `primordial`, H/He/D chemistry with atomic and H2 cooling,
composition-dependent adiabatic index `thermal.gamma = "EXACT"`). 20 hydro steps of
dt = 3e11 s, operator split:

1. **hydro** (a stand-in): a uniform compression rho -> 1.02 rho; the
   species are advected with the gas, and the internal energy gets the
   adiabatic PdV work, e -> e 1.02^(gamma-1). (No `thermal.heating = ["COMPRESS"]`:
   the compression heating is the hydro code's.)
2. **chemistry and cooling**, one KROME call over all the cells:
   - Tgas from the specific internal energy with KROME's mean molecular
     weight and adiabatic index, T = (gamma - 1) mu m_p e / k_B;
   - integrate chemistry and temperature over dt (`evolve_T`);
   - the internal energy from the new Tgas with the NEW mu and gamma (H2
     formation changes both).

The output `result.dat` has, per cell: the index, the initial n_H, the final
Tgas, n(H2)/n_H and n(e-)/n_H.

## Build and run

Generate and build the KROME project once, from the KROME-X root
(`examples/coupling/coupling.toml` holds the options):

    python3 -m kromex examples/coupling/coupling.toml --set project.outdir=build_coupling
    make -C build_coupling lib fortran shared

Then each example (`KROME` = the project directory, default
`../../../build_coupling`):

    make -C examples/coupling/fortran_ramses && examples/coupling/fortran_ramses/coupling
    make -C examples/coupling/c_enzo         && examples/coupling/c_enzo/coupling
    make -C examples/coupling/cpp_gizmo      && examples/coupling/cpp_gizmo/coupling
    python3 examples/coupling/python_post/coupling.py build_coupling/libkrome.dylib   # .so on Linux

(each writes `result.dat` in the current directory).

## What a real coupling needs, and where it is in the drivers

- **One options handle, created once** (`krome_opts_new` / `SolverOpts`):
  tolerances, `evolve_T`, the solver path. Physics is fixed when the project
  is generated.
- **Many cells per call.** `krome_solve_cells` / `krome_solve_x` /
  `solve_cells` take `n[ncell][nspec]` (one cell's species contiguous; in
  Fortran `n(nspec, ncell)`), a per-cell `dt`, and parameter rows (one row
  broadcast with `par_stride = 0`, or one row per cell). With the library
  compiled with `-fopenmp`, the cells run in parallel. The host code's own
  layout is gathered into this buffer and scattered back: see `gather` /
  `scatter` in `c_enzo`, and the vector loop in `fortran_ramses`.
- **Units.** KROME works in cgs: number densities [cm^-3] (or mass
  fractions plus rho [g cm^-3] with the `_x` calls), Tgas [K], dt [s]. The
  ENZO example converts from code units on the way in and back on the way
  out.
- **Energy <-> temperature.** Tgas goes in, and the internal energy is
  rebuilt from the returned Tgas with `krome_get_mu` and `krome_get_gamma`.
  With `thermal.gamma = "CONST"` the index is the `gamma_ad` parameter.
- **Failures.** A cell that fails is left untouched and flagged
  (`flags[c]`, `krome_flag_message`). The drivers retry it with 2, then 4
  substeps; a host code can also shorten its hydro step. A failed solve
  writes `KROME_ERROR_REPORT`, the full input of the cell, ready to replay
  (`krome_opts_set_error_file`: one file per MPI rank).
- **Per-cell physics inputs** are parameters (`NETWORK.md` of the project
  lists them): a local radiation field, cosmic-ray rate, metallicity,
  column densities, ... The C++ example keeps a `Params` per particle.
- **Conservation.** The solver conserves elements and charge to its
  tolerance. The `conserve` option projects the result back onto the input
  totals exactly, for long runs that need it.
