---
name: kromex-coupling
description: Couple KROME-X chemistry and microphysics to a hydrodynamics code (the host code, named by the user, e.g. RAMSES, ENZO, GIZMO, AREPO, FLASH, Athena++, PLUTO, Phantom, or a custom code). Use when the user asks to add KROME-X to a simulation code, to replace its cooling or chemistry module with KROME-X, or to check or debug an existing KROME-X coupling: units, passive scalars, energy and temperature, adiabatic index, radiation, cosmic rays, dust, time-stepping, CPU/GPU, build system and tests.
---

# Coupling KROME-X to a hydrodynamics code

The host code is the hydrodynamics code the user names (if they did not
name it, ask). This skill turns KROME-X into the host code's chemistry
and cooling step. Read the general skill `skills/kromex/SKILL.md` first:
it covers generating and building a KROME-X project. Paths are relative
to the KROME-X folder; the manual is `docs/manual.md`.

**Ask before assuming.** A coupling is only as correct as its units and
conventions. Whenever the host code's choice is not clear from its
source and documentation (units, what a variable stores, which energy
is evolved, how the time step is chosen), ask the user, and say what you
found and what is unclear. Never guess a unit or a convention silently.

## 1. Read the host code, then fill the coupling sheet

Before writing any code, find in the host code's source and
documentation:
- its **unit system** (code units and their conversion to cgs);
- its **hydro variables** and **equation of state**;
- how it stores **extra fields** (passive scalars, tracers, particle
  properties);
- where its existing **cooling or chemistry step** is called, if any,
  and how it is **operator-split** from the hydro step;
- how it chooses the **time step** (CFL, cooling-time limiters,
  subcycling, individual time steps);
- its **parallelization** (MPI domains, OpenMP, GPU-resident data);
- its **build system** (Makefile, CMake, configure options) and
  language.

Then write down this sheet, and show it to the user before
implementing. Every line needs an answer, either found in the code or
given by the user:

| item | host code | KROME-X side |
|---|---|---|
| units of density, velocity, energy, time, length | | cgs (manual section 1) |
| species storage | | number densities n [cm^-3], or mass fractions x with rho |
| energy variable | | Tgas [K] |
| adiabatic index | | `thermal.gamma` |
| mean molecular weight | | `get_mu` |
| radiation / photochemistry | | none, Av/G0-based, or photon bins |
| cosmic rays | | a network parameter |
| dust | | `d2g`, or binned dust |
| dust temperature needed by the host? | | `Tdust` parameter, or computed per bin |
| time step and splitting | | one call per cell per hydro step |
| hardware | | CPU (OpenMP), GPU (Kokkos) |
| build system and language | | C, C++, Fortran or Python interface |

## 2. The items, one by one

**Units (the most common source of error).** KROME-X works in cgs: n in
cm^-3, rho in g cm^-3, Tgas in K, times in s, heating and cooling in
erg cm^-3 s^-1. Convert every quantity at the interface, in one place,
with the host code's own conversion factors (never retyped constants).
Comoving cosmological codes: convert to physical (proper) density and
time, including the scale factor, before calling KROME-X. Check the
conversion with a cell of known physical density and temperature.

**Species as passive scalars.** Each KROME-X species becomes one
advected field of the host code. Find out what the host advects:
- **mass fractions** (the usual choice): call the mass-fraction
  interface (`krome_solve_cells_x` in C, `solve_cells_x` in C++,
  `krome_solve_x` in Fortran, `net.solve_x` in Python) with the cell's
  rho; `renormalize_x` rescales the result to the input sum when the
  host expects fractions summing to one;
- **partial mass densities** (rho x): divide by rho first, or convert to
  number densities with the species masses;
- **number densities**: use the number-density interface directly.

Species indices and masses come from the generated project
(`krome_net.h`, `krome_net.f90`, `krome.hpp`, `NETWORK.md`), never
hard-coded. Advection errors break element conservation: KROME-X's
`conserve` option (manual section 5) can restore the charge, or the
element totals and charge, after each call.

**Energy and temperature.** KROME-X evolves Tgas (`evolve_T`). Find out
which energy the host evolves: specific internal energy e [erg g^-1],
internal energy density rho e [erg cm^-3], or total energy (internal +
kinetic + magnetic). Then:
1. take the internal part (subtract the kinetic and magnetic energy
   from a total energy);
2. Tgas = (gamma - 1) mu m_p e / k_B, with gamma and mu from KROME-X for
   the cell's composition (`get_gamma`, `get_mu`, manual section 7) and
   KROME-X's constants (`krome::constants`, `KROME_CONST_*`,
   `krome_constants`, manual section 7);
3. after the call, rebuild e from the new Tgas with the new gamma and mu
   (chemistry changes both, e.g. when H2 forms), and add the kinetic and
   magnetic energy back.

**Adiabatic index.** It must be the same in the host's equation of
state and in KROME-X. With a constant gamma in the host, configure
`thermal.gamma = "CONST"` and set the parameter `gamma_ad` to the host's
value. With `"EXACT"` or `"POPOVAS"` (gamma from the composition and
temperature), the host's hydro solver must accept a variable gamma per
cell; if it cannot, ask the user which compromise to take.

**Mean molecular weight.** Use KROME-X's mu for every conversion between
energy and temperature, including the host's own temperature outputs
and cooling-time estimates, so that the two codes agree.

**Photochemistry.** Ask which case applies:
- **none**: no photo reactions in the network;
- **Av- or G0-based** (a field attenuated by extinction): set the
  parameters `G0` and `Av` per cell; column densities for shielding come
  from `radiation.num2col` (with `"length"`, the per-cell parameter
  `Lcol`, e.g. the cell size or a Sobolev or shielding length computed
  by the host);
- **photon bins** (the host does radiative transfer, or a given
  spectrum): configure the bins and pass each cell's field `Phi` per
  bin (manual section 8). With a radiative-transfer host, decide whether
  the field is static during the chemistry step or co-evolved with it
  (`radiation.evolve`, which conserves the photons absorbed by the
  chemistry).

**Cosmic rays.** The ionization rate is a network parameter (its name is
in `NETWORK.md`, e.g. `crate`): constant, or set per cell when the host
provides it (cosmic-ray transport, attenuation with column density).
`CR` heating uses it too.

**Dust.** Ask which description the host has:
- a **dust-to-gas ratio** (`d2g`, per cell if it varies) with one
  effective grain population;
- **size-binned dust** (manual section 9): the bins are set on the
  options, and the per-cell temperatures (and radii with `dust.evolve`)
  are passed and returned as rows; if grain sizes evolve, they are extra
  fields the host must store and advect.

**Dust temperature.** Find out whether the host needs it (dust emission,
radiative transfer, outputs). If not, the `Tdust` parameter is enough
for the rates. If yes, use binned dust: KROME-X computes each bin's
temperature from the radiative balance (or evolves it, `evolve_Td`) and
returns it per cell.

**Time-step consistency.** KROME-X integrates each cell over the host's
time step dt, with its own adaptive internal steps: do not subcycle it
yourself. But the host may need to limit its time step by the cooling
time where cooling is fast, t_cool = (rho e) / |heating - cooling|
(`get_thermal`, manual section 7), or the operator splitting loses
accuracy. Ask whether the host already has such a limiter. With
individual or hierarchical time steps, each cell passes its own dt in
the same call.

**Failed cells.** A failed cell is returned unchanged with a negative
flag. Retry it with smaller substeps, and log it; never use its output.
Each failure writes `KROME_ERROR_REPORT` with the cell's full input; in
MPI runs give each rank its own file (`error_file`).

**Hardware.**
- **CPU with MPI:** each rank calls KROME-X on its own cells.
- **CPU with OpenMP:** either KROME-X parallelizes the call (the project
  built with `-fopenmp`, called outside the host's parallel regions), or
  each host thread calls KROME-X on its own cells (the project built
  without OpenMP); not both, which would nest the parallelism.
- **GPU:** build the project with Kokkos (`make KOKKOS=...`, manual
  section 6.5) and initialize the device once. Calls of at least
  `min_device` cells (default 10^4) go to the GPU. KROME-X takes
  host-memory arrays and copies them to the device itself: a GPU-resident
  host code must copy its fields back to host memory for the call, or
  ask the user how they want to handle it.
- In every case, **gather many cells per call** (a whole patch, block or
  rank domain), never one call per cell.

## 3. Generate the KROME-X project

Choose the network and the physics with the user, from the sheet: the
coolers and heaters (leave out `COMPRESS`: the hydro already does the
PdV work), the radiation, the dust, the adiabatic index. Then generate
it (`skills/kromex/SKILL.md`, steps 1-3), and check the network
(`python3 -m kromex --check`). The project is regenerated whenever the
physics changes; the host only links it.

## 4. Implement the coupling step

One routine, called where the host's cooling or chemistry step is (or
after the hydro update, operator-split):

1. gather the cells (a patch, a block, the rank's cells) into contiguous
   arrays;
2. convert: densities or mass fractions, rho, Tgas from the energy,
   per-cell parameters (`par` rows), and the field, dust and Td rows if
   used;
3. one KROME-X call over all of them, each with its own dt;
4. check the flags; retry the failed cells with substeps;
5. convert back: species fields, internal energy from the new Tgas (new
   gamma and mu), dust fields;
6. create the options handle once at start-up and reuse it.

Start from the example in the host's language in `examples/coupling/`
(C, C++, Fortran, Python): it shows the calls, the parameters and the
error handling.

## 5. Modify the build

- **C or Fortran host:** add the project's include path, link
  `libkrome.a` and the C++ standard library (`-lc++` on macOS,
  `-lstdc++` on Linux); a Fortran host also compiles against the
  modules (`krome.mod`, `krome_net.mod`) and links `krome.o` and
  `krome_net.o` (`make lib fortran` in the project). The Makefiles in
  `examples/coupling/` show the lines.
- **C++ host:** include `krome.hpp` (header-only) and link the
  project's `net_<name>_data.o`.
- Add `-fopenmp` for OpenMP, a LAPACK for large networks
  (`LAPACK=...`), the Kokkos flags for GPUs.
- Keep the KROME-X project path as one build variable, so that a
  regenerated project is picked up without editing the host's build.

## 6. Test the coupling

1. **One cell through the host's routine** (a static or one-zone setup):
   the result must match the same cell solved by KROME-X alone (a
   standalone call with the same inputs) to the solver tolerance. This
   checks units and conversions.
2. **Conservation:** total mass, element totals and charge before and
   after the step (`get_element_totals`); mass fractions summing to one.
3. **A Sod shock tube with cooling, run in the host code with KROME-X
   coupled** (the end-to-end test, as in the KROME paper, Grassi et al.
   2014, Sect. 4.5):
   - first **without cooling** (chemistry on, `evolve_T` off): the
     hydrodynamics must match the analytic Sod solution and the run
     without KROME-X, since the chemistry only follows the gas; the
     species are advected with the flow and stay positive and
     conserved;
   - then **with cooling** (`evolve_T` on): the post-shock gas cools,
     its density rises, and the shock and the abundances (e.g. the
     ionization and H2 behind the shock) change accordingly; check that
     the energy lost matches the cooling KROME-X reports
     (`get_thermal`).

   For the expected behaviour, KROME-X ships the KROME paper's version
   of this test, a 1D spherical shock of 100 shells (Bodenheimer et al.
   2006), without and with cooling: `python3 -m kromex --test shock1D`
   and `--test shock1Dcool`, with their own hydro driver and the legacy
   profiles. These are a reference to compare with, not a substitute for
   the run in the host code.
4. **A run of the user's target setup**, small: check the failure flags,
   the time spent in KROME-X against the hydro, and the temperature
   floor and ceiling behaviour.

Report what each test compared, with the numbers, and anything that did
not match.
