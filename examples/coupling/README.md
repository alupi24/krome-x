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
| `gpu/fortran/` | CUDA Fortran | a GPU-resident host code in the style of mini-ramses: `krome_solve_cells_x_device` on device arrays (below) |
| `gpu/cuda/` | CUDA C++ | the same with a flat array of cells |
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

## GPU-resident host codes (`gpu/`)

A host code whose fields live in GPU memory calls KROME-X on its device
arrays (manual section 6.8): the data never leave the GPU. The two GPU
examples keep a grid of 32768 cells on the GPU (n_H from 1 to 1e6 cm^-3,
one of them the cell above) and, each step: gather the cells into
contiguous arrays in a kernel (mass fractions from the passive scalars
rho x, Tgas from the internal energy), restore their element composition
(consistent advection, `krome_conserve_cells_x_device`), integrate them
(`krome_solve_cells_x_device`, one call), and scatter the results back in
a kernel (rho x, and the internal energy from the new Tgas and
composition). They write the quick start's cell to `abundances.dat`, so
`compare.py` checks them against the CPU examples (1e-4: the device's
arithmetic).

They need an NVIDIA GPU, the CUDA toolkit (and the NVIDIA HPC SDK for
CUDA Fortran), and a Kokkos installation for CUDA (`KOKKOS`). The project
is built with Kokkos, its Fortran modules by nvfortran:

    python3 -m kromex examples/coupling/options.toml --set project.outdir=build_coupling_gpu
    make -C build_coupling_gpu lib fortran KOKKOS=$KOKKOS \
      CXX=$KOKKOS/bin/nvcc_wrapper FC=nvfortran \
      CXXFLAGS="-O2 -std=c++17 -arch=sm_80 --expt-extended-lambda --expt-relaxed-constexpr"

(`-arch=sm_80` / `-gpu=cc80` in the examples' Makefiles: an A100; add
`-fopenmp` to CXXFLAGS when Kokkos runs OpenMP on the host). Then

    cd examples/coupling/gpu/fortran && make KOKKOS=$KOKKOS && ./example
    cd examples/coupling/gpu/cuda    && make KOKKOS=$KOKKOS && ./example
    python3 examples/coupling/compare.py

With a module system, load the GCC that built Kokkos after the HPC SDK
(its libstdc++ must come first). The coupling of a GPU-resident code,
step by step: `skills/kromex-coupling-gpu/SKILL.md`.

## From here to a hydro code

A hydro code makes the same call once per step, over all its cells at
once: `krome_solve_cells` (C), `solve_cells` (C++), `krome_solve` with
`n(nspec, ncell)` (Fortran), `net.solve` with `n[ncell, nspec]`
(Python), with per-cell parameters (here the free-fall time `tff`) as
rows of `par`. The manual (docs/manual.md, sections 5 and 6) covers the
options, the per-cell inputs, the GPU path and the failure flags.
