---
name: kromex-coupling-gpu
description: "Couple KROME-X chemistry and microphysics to a hydrodynamics code whose fields live in GPU memory (a GPU-resident host code: CUDA C++, CUDA Fortran, HIP, OpenACC or Kokkos, e.g. mini-ramses), so that the chemistry runs on the device data without copies through host memory. Use when the user asks to add KROME-X to such a code, or to check or debug that coupling: device arrays and layouts, the gather and scatter kernels, units and energy, device selection under MPI, memory, failed cells, the build with Kokkos and the GPU compilers, tests. For CPU codes, or GPU codes that keep their fields in host memory, use kromex-coupling instead."
---

# Coupling KROME-X to a GPU-resident hydrodynamics code

The host code is the hydrodynamics code the user names (ask if they did
not). It keeps its state on the GPU between steps; KROME-X integrates the
chemistry of its cells there, on device arrays, with its device-resident
calls (manual section 6.8). Paths are relative to the KROME-X folder; the
manual is `docs/manual.md`. Read `skills/kromex/SKILL.md` first (making a
KROME-X project), and section 2 of `skills/kromex-coupling/SKILL.md` for
the physics of any coupling (units, energy, adiabatic index, radiation,
cosmic rays, dust, time step): this skill adds what is particular to
device data and does not repeat it.

**Ask before assuming.** Units, what each variable stores, which energy
is evolved, how the arrays are laid out in device memory, which GPU and
which stream: when the host code's source and documentation do not
settle one, ask the user, say what you found and what is unclear. Never
guess a unit, a layout or a device convention.

## 0. Is this the right skill?

Find where the host code's fields are during its time step:
- in **GPU memory**, copied to the host only for output or communication
  (CUDA / HIP kernels, CUDA Fortran `device` arrays, OpenACC `present`
  data, Kokkos views in a device space): this skill;
- in **host memory** (a CPU code; a code that offloads only some
  kernels and copies back each step): `skills/kromex-coupling`, whose
  host-array calls move the cells to the GPU themselves.

## 1. Read the host code

Besides the physics items of the general coupling skill, find:
- its **GPU programming model** (CUDA C++, CUDA Fortran, HIP, OpenACC,
  Kokkos) and the **GPU runtime** (CUDA or HIP): KROME-X must be built
  for the same one;
- the **layout of the state in device memory**: one array per variable
  (`u(ncell, nvar)`), per cell (`u(nvar, ncell)`), or blocked (octs of
  8 cells, `uold(8, nvar, noct)`, patches with ghost zones); where the
  species (passive scalars) are and what they store (rho x, x, n);
- how it **lists the cells to update** (leaf cells, active cells of a
  level, a patch interior);
- its **kernels' conventions**: block sizes, streams (the default
  stream or its own), where it synchronizes;
- **MPI and GPUs**: how a rank selects its GPU (e.g.
  `cudaSetDevice(mod(rank, ngpu))`);
- the **device memory** it uses and how much is left for KROME-X;
- its **build**: compilers (nvcc, nvfortran, hipcc, ...), CUDA or ROCm
  version, the GPU architecture (`sm_80`, `cc80`, `gfx90a`), MPI wrappers.

## 2. The coupling sheet

Fill this with the physics lines of the general sheet and show it to the
user before implementing:

| item | host code | KROME-X side |
|---|---|---|
| GPU model and runtime | | Kokkos built for CUDA or HIP (the same runtime) |
| GPU architecture | | Kokkos' arch flag (`-arch=sm_80`, ...) |
| GPU of each MPI rank | | `krome_device_initialize_id(id)` (or the host's own Kokkos) |
| cells per call | | all of a rank's cells to update, in one call |
| species storage | | mass fractions x (with rho) or number densities n |
| array layout | | species-major `x(ncell, nspec)` (default) or cell-major |
| energy variable | | Tgas [K], converted in the host's kernels |
| per-cell parameters | | rows `par(npar, ncell)` on the device, or one row |
| photon field, Td, grain radii | | rows `[ncell][bins]` on the device, optional |
| consistent advection | | reference element totals, one row or per cell |
| stream | | KROME-X fences the device before and after the call |
| failed cells | | flags on the device: 0, +1 recovered on the host, < 0 failed |
| device memory for KROME-X | | its buffers: `device_memory` (default 80 % of free) |
| build | | `libkrome.a` by nvcc_wrapper / hipcc, the host's Fortran modules |

## 3. Build KROME-X for the device

1. **Kokkos** for the host code's GPU: an installation with the CUDA (or
   HIP) backend and the right architecture, built with a GCC whose
   libstdc++ the whole build will use. Kokkos is KROME-X's dependency,
   not the host code's: a host code that does not use Kokkos only links
   `libkokkoscore`.
2. **The project** (`skills/kromex/SKILL.md`), then the library with the
   Kokkos compiler, and the Fortran modules with the host's Fortran
   compiler when it is Fortran:

       make lib fortran KOKKOS=$KOKKOS CXX=$KOKKOS/bin/nvcc_wrapper FC=nvfortran \
         CXXFLAGS="-O2 -std=c++17 -arch=sm_80 --expt-extended-lambda --expt-relaxed-constexpr"

   (add `-fopenmp` when Kokkos runs OpenMP on the host).
3. **The host code's link line**: `libkrome.a` (and, from Fortran,
   `krome.o krome_net.o`), `-L$KOKKOS/lib64 -lkokkoscore`, `-lstdc++
   -ldl`, `-lgomp` with Kokkos' OpenMP host space; the GPU runtime comes
   with the host's GPU compiler. The Makefiles of `examples/coupling/gpu/`
   show the lines.
4. **Pitfalls** (seen on a cluster with modules): an HPC SDK module can
   put an older GCC first, and every link then fails on `GLIBCXX_*`
   symbols: load the GCC that built Kokkos last (or point the linker at
   its libstdc++). The CUDA versions of Kokkos' build and of the Fortran
   compiler's CUDA should be compatible. Keep the KROME-X project path in
   one build variable.

## 4. The coupling step

Start from the example in the host's language, `examples/coupling/gpu/`
(`fortran/`: CUDA Fortran, octs `uold(8, nvar, noct)` and a leaf list,
in the style of mini-ramses; `cuda/`: CUDA C++, one array per variable).
Both run the steps below and check themselves against the CPU examples.

**Once, at start-up:**
- `krome_device_initialize_id(rank % ngpu)`, the GPU the rank already
  uses (a Kokkos host code initializes Kokkos itself and skips this);
- the options handle (`krome_opts_new`, `krome_set`), kept for the run;
- KROME-X's argument arrays on the device, sized for the most cells a
  call will have: species `x(ncell, nspec)` (species-major; a leading
  dimension `ld >= ncell` allowed), `rho`, `Tgas`, `dt`, the parameter
  rows, `flags`, the optional rows (Phi, Td, adust);
- the reference element totals for the consistent advection
  (`krome_get_element_totals_x` of the initial composition: one row, or
  one per cell, e.g. scaled with the metallicity), copied to the device.

**Each hydro step** (after the hydro update, operator-split):
1. **Gather kernel** (the host's own, one thread per cell): from the
   host's cells into KROME-X's arrays, converting on the way: code units
   to cgs; the passive scalars to x (or n); the internal energy (the
   total minus the kinetic and magnetic parts) to Tgas = (gamma - 1) rho
   e / (k_B n_tot), with gamma and n_tot (the species' n summed) of
   KROME-X's network; the per-cell parameters (e.g. the free-fall time,
   the column density, the cosmic-ray rate) into their rows.
2. **Consistent advection**: `krome_conserve_cells_x_device` (or `_device`
   for n) back to the reference composition and neutrality. Its flags:
   -112 an element without carriers (corrected as far as possible), +2 a
   species changed by more than 50 % (a warning).
3. **The chemistry**: `krome_solve_cells_x_device` (or
   `krome_solve_cells_device`), every cell over its own dt, one call.
4. **Scatter kernel**: back into the host's cells for the cells with
   flag >= 0: x to the passive scalars, and the internal energy from the
   new Tgas and the new n_tot (the chemistry changes the number of
   particles, e.g. when H2 forms), the kinetic and magnetic parts added
   back; code units again.

**Device pointers.** CUDA C++/HIP: the raw pointers. CUDA Fortran:
`c_loc` of the `device` arrays (as mini-ramses passes them to CUB; not
`c_devloc`). OpenACC: `c_loc` inside `!$acc host_data use_device(...)`.
Kokkos: `view.data()` of views in the device space (the host code and
KROME-X must use the same Kokkos installation).

**Synchronization.** The calls fence the device before and after: the
gather kernel is complete when KROME-X reads, and the scatter kernel sees
the results. A host code on its own streams needs nothing more; one that
overlaps other work must not touch these arrays during the call.

**Memory and size.** One call over all the cells to update: small calls
are latency-bound (on an A100 the device beats a CPU node from a few
thousand cells per call). KROME-X allocates its own buffers per cell (a
few kB for a small network, MB for hundreds of species) and splits a
call into chunks that fit 80 % of the free memory; when the host code
allocates more GPU memory later, bound KROME-X with `device_memory`
(bytes) or free its cache (`device_cache` 0).

## 5. Failed cells

A cell the device cannot integrate is re-solved on the host from its
inputs: flag +1 (`krome_flag_host_recovered`), its results written as
any other's. A cell that fails there too (e.g. out of steps,
`max_steps`) keeps its input, flag < 0; the summary (`device`, `host`,
`failed`) counts them. The host code decides: the scatter kernel skips
them; then retry them (compact their indices and call KROME-X again on
those cells with a smaller dt, or a larger `max_steps`), or accept them
unchanged for this step and log it. Never use a failed cell's output.

## 6. Test the coupling

1. **The examples on the target machine**: build and run
   `examples/coupling/gpu/` and `compare.py` (the GPU examples against the
   CPU ones, 1e-4): it checks the toolchain, the link and the device.
2. **One cell through the host's kernels**: a cell of known density,
   composition and temperature in the host's units, gathered, solved and
   scattered, against the same cell solved by KROME-X's host interface
   (roundoff, ~1e-5): it checks the units and the energy conversion.
3. **Many cells**: the host's real cells through the coupling against the
   host-array call on the same cells (`krome_solve_cells_x`): device
   against host to ~1e-5 - 1e-4; the same cells failing.
4. **Conservation**: element totals, charge and mass before and after
   the step (the device-side consistent advection keeps them).
5. **A Sod shock tube with cooling, run in the host code** (the general
   coupling skill, section 6, item 3).
6. **Cost**: the time of the KROME-X call against the hydro step, the
   device memory, the failed cells per step.

Report what each test compared, with the numbers, and anything that did
not match.
