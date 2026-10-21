---
name: kromex
description: "Set up KROME-X chemistry and microphysics for a simulation or a one-zone calculation, from a chemical network to the calls in the host code. Use when choosing or writing a network (.knet, KIDA, UMIST), configuring the physics (cooling, heating, radiation, dust), generating and building a KROME-X project, calling it from C, C++, Fortran or Python, coupling it to a hydrodynamics code (many cells per call, GPU), or diagnosing failed solves."
---

# KROME-X: from a network to a simulation

KROME-X generates a chemistry solver for one network and one choice of
physics, then the host code calls it on its cells. Work through the
steps below in order. Paths are relative to the KROME-X folder; the full
reference is `docs/manual.md`, and the sections named here are its
sections.

## Rules that are easy to get wrong

- **Physics is fixed at generation; numerics are runtime.** Coolers,
  heaters, shielding, photon bins, dust and the adiabatic index are
  chosen when the project is generated. Changing them means
  regenerating. At run time only the solver options (tolerances,
  `evolve_T`, the linear solver, the CPU/GPU path), the per-cell
  parameters, the photon field and the dust bins change.
- **Units are cgs:** number densities in cm^-3, temperatures in K,
  times in s, rho in g cm^-3, photon fluxes in photons cm^-2 s^-1 per
  bin, bin energies in eV, heating and cooling in erg cm^-3 s^-1.
- **Never hard-code species or parameter indices.** Use the generated
  names: `IDX_H2`, `KROME_IDX_H2`, `krome_idx_H2`, `net.idx("H2")`;
  `KROME_PAR_tff`, `krome_par_tff`, `par={"tff": ...}`. The species
  order is listed in the project's `NETWORK.md`.
- **Required parameters have no default.** A parameter marked REQUIRED
  in `NETWORK.md` (e.g. `tff` with the COMPRESS heating) must be set
  for every cell, otherwise the solve returns -105 and nothing is
  integrated.
- **A failed cell is left untouched.** Its inputs are returned
  unchanged and its flag is negative. Always check the return value or
  the per-cell flags, and print them with `flag_message` /
  `krome_flag_message` / `net.flag_message` (section 10).
- **Do not edit generated files.** Only `Makefile` and the examples are
  yours; everything else is rewritten when the project is regenerated.
  The project records its configuration in `krome.toml`, and
  `python3 -m kromex krome.toml` regenerates it.

## 1. Choose or write the network

- **Shipped networks:** `networks/*.knet`, used by name
  (`network = "primordial"`). `docs/network-format.md` describes the
  format for a new one.
- **From the databases:** `kida2krome FILE out.knet` (KIDA) and
  `umist2krome FILE out.knet` (UMIST RATE22) (manual section 2).
- **Check it before anything else:**
  `python3 -m kromex --check config.toml` (or `kp.check()`) compares
  the generated rates, right-hand side and Jacobian with an
  independent evaluation, and flags NaN or suspicious rates and species
  that are never formed or destroyed. Exit status 1 on a failure.
- **Read it:** `python3 -m kromex_doc NETWORK` writes HTML pages with
  every reaction, its rate plotted against temperature, the species,
  graphs, and possibly missing reactions (manual section 2).

## 2. Configure the physics

The configuration is a TOML file or a Python object, with sections
`project`, `solver`, `rates`, `thermal`, `radiation`, `dust`, `surface`
(manual section 3; `python3 -m kromex --template FILE` writes every
option with its description).

```toml
[project]
network = "primordial"
outdir = "chem"
bindings = ["fortran", "python"]   # C and C++ are always available

[thermal]
cooling = ["ATOMIC", "H2"]
heating = ["CHEM"]
```

```python
from kromex import Preprocessor
kp = Preprocessor(network="primordial")
kp.thermal.cooling = ["ATOMIC", "H2"]
kp.thermal.heating = ["CHEM"]
```

- Every option is validated when set; a misspelled name is an error
  that suggests the right one.
- **In a hydro code, leave out `COMPRESS` heating**: the hydro step
  already does the PdV work. It is for one-zone collapse models.
- Choose only the processes the problem needs: each one costs time in
  every right-hand-side evaluation.

## 3. Generate and build

```sh
python3 -m kromex config.toml            # or kp.run("chem")
make -C chem lib                         # libkrome.a: the C interface
make -C chem fortran                     # Fortran modules (if bound)
make -C chem shared                      # libkrome.so/.dylib for Python
make -C chem check                       # this build's rates vs reference
```

- From Python, `net = kp.build()` builds the shared library and loads
  it.
- Large networks (hundreds of species): build with an optimized LAPACK
  (`make LAPACK=-lopenblas`, MKL, or `"-framework Accelerate"`), never
  the reference Netlib BLAS. Very sparse networks use the generated
  sparse LU automatically.
- GPU: `make lib KOKKOS=<kokkos prefix> CXX=<the Kokkos compiler>`, e.g.
  `nvcc_wrapper` for CUDA (manual section 6.5).
- Generation needs no network access, and the data tables are compiled
  into the code: compute nodes need nothing beyond the executable.

## 4. Call it from the host code

`examples/coupling/` has the same problem in C (`c/`), C++ (`cpp/`),
Fortran (`fortran/`) and Python (`python/`), each with a Makefile.
Start from the one in the host code's language. Minimal calls:

```c
/* C: krome.h, krome_net.h, link libkrome.a and the C++ runtime */
double par[KROME_NPAR];
krome_params_default(par);                  /* REQUIRED ones are NaN */
krome_opts* h = krome_opts_new();
krome_opts_set(h, "evolve_T", 1);
int r = krome_solve_cells(ncell, n, Tgas, dt, par, 0, NULL, flags, h);
```

```fortran
! Fortran: use krome; use krome_net; n(krome_nspec, ncell)
call krome_params_default(par)
h = krome_opts_new()
call krome_set(h, "evolve_T", 1d0)
ierr = krome_solve(n, Tgas, dt, h, par=par, flags=flags)
```

```python
import krome                       # in the project folder (or kp.build())
net = krome.load()
n = net.densities({"H": 1e4, "H+": 1.0, "e-": 1.0, "He": 775.0})
res = net.solve(n, Tgas=100.0, dt=3.15e10, opts=net.options(evolve_T=True))
```

C++ uses `solve` / `solve_cells` from `krome.hpp` with `Params` and
`SolverOpts` (manual section 6.1).

**Coupling to a hydrodynamics code** (operator splitting):

1. After the hydro step, gather the cells and make **one call for all of
   them** (`krome_solve_cells`, `solve_cells`, `krome_solve` with
   `n(nspec, ncell)`, `net.solve` with `n[ncell, nspec]`), each cell
   over its own `dt`.
2. Convert the code's variables: species mass densities or passive
   scalars to number densities (or use the mass-fraction calls,
   `krome_solve_cells_x` / `solve_x`, with the cell's rho); the
   temperature from the internal energy with KROME-X's mean molecular
   weight and adiabatic index (`get_mu`, `get_gamma` and their C,
   Fortran and Python forms, manual section 7).
3. With `evolve_T`, rebuild the internal energy from the returned Tgas
   with the *new* mu and gamma (chemistry changes both).
4. Cell-dependent quantities (dust-to-gas ratio, cosmic-ray rate, Av,
   free-fall time) are per-cell parameters: one `par` row per cell.
5. Handle failed cells: their state is unchanged; retry them with
   substeps or relaxed tolerances, and log the flag. Each failure also
   writes `KROME_ERROR_REPORT`; in MPI runs give each rank its own file
   (`error_file` option).
6. Create the options handle once and reuse it.

**Many cells and GPUs** (manual section 6.5): with OpenMP
(`CXXFLAGS += -fopenmp`) the cells of a call run in parallel on the CPU.
In a Kokkos build, initialize the device once (`Kokkos::initialize`,
or `krome_device_initialize()` in C and Fortran, `net.device_initialize()`
in Python; finalize before exiting): until then every call runs on the
CPU. Calls of at least `min_device` cells (default 10^4) then go to the
GPU, smaller ones stay on the CPU; a GPU pays off only with many cells
per call. Calls larger than the GPU memory are split into chunks
automatically (`device_memory`).

## 5. Validate before production

- `make check` on the target machine (compiler and flags included).
- A one-zone run of a representative cell, compared with expectations
  or with a known result: the shipped tests (`python3 -m kromex
  --tests`, `--test NAME`) cover collapses, clouds, photon-dominated
  regions, dust and the early Universe.
- Inspect the chemistry and thermal balance with the queries (manual
  section 7): `get_thermal` for heating and cooling per process,
  `get_fluxes` and `print_best_flux` for the dominant reactions.
- For the data behind the physics and what to cite, see
  `data/data_sources.md`.

## Diagnosing a failed solve

| flag | meaning | what to do |
|---|---|---|
| -105 | a REQUIRED parameter is unset | set it for every cell (`NETWORK.md`) |
| -103 | photo reactions without photon bins or a field | set the bins and `Phi` (manual section 8) |
| -102 | dust cooling without dust bins | set the dust bins (manual section 9) |
| -1 | step budget exhausted | a stiff or very long step: raise `max_steps`, or split `dt` |
| -3, -4 | repeated error-test or Newton failures | check the inputs (negative or NaN densities, extreme T), the network (`--check`), then loosen `rtol`/`atol_n` |
| -110 | Tgas went non-positive | the cooling outruns the heating: check the thermal processes and units |

Every flag is listed in manual section 10, and `KROME_ERROR_REPORT`
records the failing cell's full input state.
