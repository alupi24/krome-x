# KROME-X user manual

This manual covers every generation option, every runtime option and
parameter, every return flag, and how to call KROME from C++, C,
Fortran and Python. The declarations themselves are documented in the
headers:
- `krome.hpp`: the C++ API;
- `krome.h`: the C API;
- `krome.f90`: the Fortran module;
- `krome.py`: the Python module.

Other references:
- the `.knet` network format: [network-format.md](network-format.md).

Contents:

1. [Workflow](#1-workflow)
2. [Generating a project](#2-generating-a-project)
3. [Generation options](#3-generation-options)
4. [Runtime parameters](#4-runtime-parameters)
5. [Solver options](#5-solver-options)
6. [Calling KROME](#6-calling-krome)
7. [Queries](#7-queries)
8. [Photon bins and radiation](#8-photon-bins-and-radiation)
9. [Binned dust](#9-binned-dust)
10. [Return flags and the error report](#10-return-flags-and-the-error-report)
11. [Migrating from legacy KROME](#11-migrating-from-legacy-krome)

## 1. Workflow

1. **Write or choose a network**: a `.knet` file in `networks/`. All
   legacy networks are already converted. A legacy `react_*` file can
   be given directly, and the generator converts it.
2. **Configure and generate a project for it**, choosing the physics:
   coolers, heaters, shielding, adiabatic index, evolved blocks. The
   physics is fixed in the generated code. The configuration is a Python
   object (`kromex.Preprocessor`) or a TOML file, and every project
   records it as `krome.toml`.
3. **Build the project** (`make`) and call KROME from your code. At
   runtime you set only:
   - the solver numerics;
   - `evolve_T`: fixed-Tgas chemistry, or Tgas evolved;
   - the per-cell parameters;
   - the photon field and the dust bins.

Units are cgs throughout:

| quantity | unit |
|---|---|
| number densities n | cm^-3 |
| temperatures | K |
| times | s |
| rho | g cm^-3 |
| photon flux per bin Phi | photons cm^-2 s^-1 |
| bin energies | eV |
| heating and cooling rates | erg cm^-3 s^-1 |

## 2. Generating a project

A project is generated from a configuration: the network and the options
of section 3, grouped in sections (`project`, `rates`, `thermal`,
`radiation`, `dust`, `surface`, `solver`). It can be written in Python or
in a TOML file, and run from either.

**Installation.** `pip install .` in the KROME-X folder installs the
package `kromex` (with numpy) and the command `kromex`; `pip install -e .`
does the same in editable mode, reading the folder itself. Without
installation, the package works from the KROME-X folder (`v2/`): Python
started there (interactively, in a notebook, or `python3 -m kromex`)
finds it, and a script elsewhere adds the folder itself,
`sys.path.insert(0, "/path/to/krome-x")`, or is run with
`PYTHONPATH=/path/to/krome-x`. Python 3.8 or newer. An installed package
extracts the ZCIE tables (the `Z_CIE` cooling) into the user cache,
`~/.cache/kromex` (or `$XDG_CACHE_HOME/kromex`), on first use.

### From Python

```python
from kromex import Preprocessor

kp = Preprocessor(network="primordial")   # a shipped network, or a path
kp.thermal.cooling = ["ATOMIC", "H2"]
kp.thermal.heating = ["CHEM", "COMPRESS"]
kp.thermal.gamma = "EXACT"
kp.run("myproj")          # generate the project
net = kp.build()          # make shared, and the loaded Network (section 6.4)
```

- **Options** are attributes of their section (`kp.thermal.cooling`), or
  `kp.set("thermal.cooling", [...])` for a name built in a script.
  `kp.get(name)` reads one, `Preprocessor.help("thermal.gamma")`
  describes it, and `print(kp)` lists those that differ from the
  defaults.
- **Every value is checked when it is set.** An unknown section or key
  is an error that names the closest ones ("did you mean
  thermal.cooling?"), so a misspelled option can never be ignored; a
  value outside its choices, or of the wrong type, is an error listing
  the valid ones. Choices are matched regardless of case.
- **Files:** `Preprocessor.from_file("run.toml")` reads a configuration,
  `kp.save("run.toml")` writes one, `kp.to_toml()` returns its text;
  `kp.to_dict()` and `Preprocessor.from_dict(d)` convert to and from a
  dictionary `{section: {key: value}}`; `Preprocessor.template()` is the
  commented template of every option.
- **The shipped tests:** `Preprocessor.from_test("collapse")` (section
  6.6), `Preprocessor.tests()` lists them.
- **`kp.check()`** checks the network with these options and writes
  nothing (below); **`kp.run(outdir)`** generates the project (default:
  `project.outdir`, else `generated/`, or `build_NAME` for a test) and
  returns its folder; **`kp.build()`** compiles its shared library (adding
  the Python binding if the project has none) and returns the loaded
  `Network`.
- Errors are exceptions: `kromex.ConfigError` (an invalid option),
  `kromex.GenerationError` (a network or a combination of options the
  generator rejects, with the reason).

### From a TOML file

```toml
# myrun.toml
[project]
network = "networks/primordial.knet"
outdir = "myproj"
bindings = ["python"]

[thermal]
cooling = ["ATOMIC", "H2"]
heating = ["CHEM", "COMPRESS"]
gamma = "EXACT"
```

```sh
python3 -m kromex myrun.toml
python3 -m kromex myrun.toml --set thermal.gamma=CONST --set project.outdir=other
```

- **The command line** runs a file, or a shipped test:

  | command | does |
  |---|---|
  | `python3 -m kromex FILE.toml` | generate the project of the file (installed: `kromex FILE.toml`, and likewise below) |
  | `--set SECTION.KEY=VALUE` | override an option (repeatable); VALUE in TOML syntax (`1e-6`, `true`, `["H2", "Z"]`), a bare word is a string (`--set thermal.gamma=EXACT`) |
  | `python3 -m kromex --test NAME` | a shipped test (section 6.6); `--tests` lists them |
  | `python3 -m kromex --check FILE.toml` | check the network, write nothing (below) |
  | `python3 -m kromex --template FILE` | write the commented template |
  | `python3 -m kromex --version` | the KROME-X version (`v2/VERSION` plus the git revision, e.g. `2.0.0+g04bb53f`; the generated code carries it: C++ `krome::KROME_VERSION`, C `krome_version()`, Fortran `krome_version_f()`, Python `net.version`) |

  There are no per-option flags: the options are the file's, changed with
  `--set`.
- **The template:** `options_template.toml`, in the root of `v2/`, lists
  every option commented out, with its description, default and an
  example value: copy it and uncomment what you need
  (`python3 -m kromex --template FILE` writes the same file).
- **Paths** in a file (the network, coolant and tolerance files) are
  relative to the file's folder.
- **A legacy `react_*` network** can be given as `project.network`: it
  is converted into the project as `<name>.knet` with
  `kromex/convert_legacy.py`, and the project is generated from that
  file. Legacy option files are not read: section 11 maps each legacy
  option to its key.

**The record.** Every project writes `krome.toml`: every option it was
generated with (defaults included), the network as an absolute path, and
a `[meta]` table with the KROME-X version and the SHA-256 of the network
file. `python3 -m kromex myproj/krome.toml`, or
`Preprocessor.from_file("myproj/krome.toml").run()`, regenerates the
project in place; a network that changed since is reported.

### Checking a network

```sh
python3 -m kromex --check myrun.toml
```

or `kp.check()`, checks a network with the options it will be generated
with, and writes nothing. It reports:

- **The generated code** against an independent evaluation of the
  network (a Python interpreter of the parsed network): the rate
  coefficients at nine temperatures from 12 K to 8e7 K, the right-hand
  side from the stoichiometry, the analytic Jacobian against finite
  differences, the freeze directives. A difference is a failure of the
  code generator: please report it. Photo, surface and tabulated rates
  are computed at run time and not evaluated independently.
- **The rates**, at the same temperatures inside each reaction's T range
  (the solver evaluates them from T_CMB to 1e9 K):
  - a rate that is NaN or infinite: a failure;
  - a negative rate, or a two-body rate above 1e-5 cm^3 s^-1 or a
    three-body rate above 1e-24 cm^6 s^-1: a warning. These are often
    fits used beyond their validity (declare their T range), or a typo
    in an exponent or the units.
- **The network:** species formed but never destroyed, destroyed but
  never formed, or in no reaction; the same reaction written more than
  once with the same rate (counted twice). Reactions written more than
  once with different rates are summed as channels (a note).

Parameters are set to a probe value (1.3e-17, Tdust 10 K): a problem
with a rate that depends on a parameter is a warning only. The exit
status is 1 on a failure, so the check can run in a script; `kp.check()`
returns False.

### Documenting a network

```sh
python3 -m kromex doc NETWORK [-o DIR] [--range Av=0.1,20] [--set crflux=1.3e-17]
```

(`kromex doc` once installed; in Python `from kromex.doc import
document`) writes the documentation of a network from its `.knet` file
alone, without generating or building it (legacy `docmake`). NETWORK is
a file or a shipped network by name; DIR defaults to `doc_<network>`.
Open `DIR/index.html` in a browser:

- **Species:** mass, charge, composition, phase, the formation enthalpy
  at 298.15 K when the thermochemical tables have it, and the reactions
  forming and destroying each.
- **Reactions:** the rate as a formula, its line in the network file, its
  T range and reference, and k against T (and against the `--range`
  quantities), with the other pieces of a piecewise rate. The plot also
  shows KROME-X's rate outside the range (`--outside`, as the generation
  option `rates.outside`; default: zero below, the bound's value above)
  and, faint, the formula extrapolated.
- **Photo reactions:** the cross section sigma(E) from the source the
  generator reads (the `xsec` expression, the Verner et al. 1996 fits of
  `xsec=auto`, the Leiden and SWRI files in `data/database`), its energy
  window, the shielding factor, and the unshielded rate in the Draine
  (1978) field. A missing data file is flagged: the generated code gives
  that reaction a zero rate.
- **Piecewise rates:** the jump of each piecewise rate at its bounds, in
  red above 0.1%: pieces that do not meet are usually fits joined at the
  wrong temperature, or a typo.
- **Missing reactions and branches:** reactant sets the network never
  uses, and products its reactant sets could give but do not (same
  elements and charge, one or two species of the same phase), each with
  its reaction enthalpy (negative: exothermic) when the tables have
  every species. A lead for the author, not a verdict. Above
  `--max-missing` species per phase (150) only the branches are listed.
- **Graphs:** one per element, from each reactant to each product with
  it (`graphs/*.dot`; rendered by Graphviz when installed, else in the
  browser).
- **LaTeX:** `latex/<network>_table.tex`, the reaction table for a paper
  (longtable: number, reaction, rate coefficient, T range, the `ref`
  citation), and `latex/<network>_reactions.tex`, a document around it
  with the definitions of the vars written by name. Run LaTeX two or
  three times: longtable settles its column widths over the runs.

The rates are evaluated in Python from the network's expressions, at the
parameters' defaults unless `--set` gives a value (a rate proportional to
a parameter whose default is 0, such as a cosmic-ray rate, is zero until
set). Surface rates, rates from the database, tables and extern
functions are computed by the generated code only: their reactions are
documented without plots (the photo reactions: their cross sections). The formulas and the graphs
without Graphviz are drawn by MathJax and viz.js, loaded from a CDN: the
site needs an internet connection for them.

### Networks from KIDA and UMIST

```sh
python3 tools/kida2krome.py gas_reactions_kida.uva.2024.in kida2024.knet
python3 tools/umist2krome.py rate22_final.rates rate22.knet
```

convert a KIDA network (the kida.uva files) or a UMIST Database for
Astrochemistry ratefile (RATE22, RATE12) to a `.knet`, with the
databases' rate formulas (the docstring of each tool cites them):

- **Two-body rates:** Kooij, alpha (T/300)^beta exp(-gamma/T); KIDA's
  ion-polar formulas (ionpol1, ionpol2) as written in KIDA.
- **Cosmic rays:** alpha times the parameter `crflux`, the ionisation
  rate in the database's convention (KIDA: zeta of H2, default 1.3e-17
  s^-1; UMIST: CP alpha zeta/zeta0 and CR alpha (T/300)^beta gamma /
  (1 - omega) zeta/zeta0, zeta0 = 1.36e-17 s^-1, omega 0.5 with
  `--albedo`; default zeta0).
- **UV photorates:** alpha exp(-gamma Av) times `uv_scale`, the field in
  units of the database's (Draine 1978 for KIDA), with `Av` the implicit
  parameter.
- **Temperature ranges:** as in the database. A rate given over several
  ranges is a piecewise rate, its gaps closed at the start of the next
  range; a rate measured at one temperature has no range. With the
  default `rates.outside = "asym"` a rate is zero below its range:
  `rates.outside = "clamp"` keeps it for colder gas.
- **Duplicates:** KIDA lists alternative values for the same reaction
  and range (the user chooses one): the best recommended is kept. UMIST
  lists sums of fits as separate lines, kept (summed); the same rate
  written twice is kept once. The others are comments in the output.
- **Species:** isomer prefixes become labels (`l-C3H` is `C3H_l`); the
  pseudo-species (cosmic rays, photons) leave the reaction. Formulas the
  network cannot express (KIDA's grain charging and surface H) are
  comments.

Each reaction carries a comment with its database number (UMIST: its
DOI as `ref`). The tools flag entries whose parameters do not fit their
formula (KIDA 2024 has a few, e.g. ionpol parameters under the Kooij
formula, giving k = 1 cm^3 s^-1): check them before use.
`kromex --check` then tests the network; `kromex doc` documents it.

### The project directory

| file | what | regenerated |
|---|---|---|
| `krome.hpp` | the C++ API (the only header a C++ caller includes) | yes |
| `krome.h`, `krome_net.h` | the C API; the network constants (`KROME_NSPEC`, `KROME_IDX_<sp>`, `KROME_PAR_<p>`) | yes |
| `krome.f90`, `krome_net.f90` | Fortran modules `krome`, `krome_net` (`project.bindings` with `fortran`) | yes |
| `krome.py` | Python module (`project.bindings` with `python`, or added by `Preprocessor.build()`) | yes |
| `krome_network.hpp`, `net_<name>*.hpp`, `net_<name>_data.cpp` | the generated network | yes |
| `core/` | the implementation headers (internal) | yes |
| `core/cvode/` | the integrator (CVODE's BDF, ported) and the SUNDIALS license (internal) | yes |
| `krome.cpp` | the C interface implementation (`libkrome.a`) | yes |
| `NETWORK.md` | this network: options used, species, parameters, reactions | yes |
| `krome.toml` | the configuration this project was generated with (section 2): `python3 -m kromex krome.toml` regenerates it | yes |
| `krome_check.cpp` | `make check`: reference rate coefficients of this network | yes |
| `Makefile`, `example.cpp`, `example.f90`, `example.py` | build and examples, yours to edit (for a shipped test: its driver instead of the examples, section 6.6) | never overwritten |

Makefile targets:

| target | builds |
|---|---|
| `make` | the library and the C++ example (for a shipped test: its driver), plus `shared` with the Python binding |
| `make lib` | `libkrome.a` (C interface + network) |
| `make example` | the C++ example |
| `make fortran` | the Fortran modules |
| `make example_f90` | the Fortran example |
| `make shared` | `libkrome.so` / `.dylib` for Python (in every project) |
| `make check` | every target above, then `krome_check`: the rate coefficients of this build (your compiler and flags) against reference values computed when the project was generated, by an independent evaluation of the network; exit status 1 on a difference above 1e-8 |

Makefile variables:

| variable | meaning |
|---|---|
| `CXX`, `FC` | compilers |
| `CXXFLAGS` | default `-O3 -std=c++17 -ffp-contract=fast -fopenmp-simd`; add `-fopenmp` for cells in parallel; bitwise-reproducible builds: `-O2 -std=c++17 -ffp-contract=off` |
| `LAPACK="-framework Accelerate"` / `LAPACK=-lopenblas` | LAPACK dense LU (`linsol`), recommended for large networks |
| `KOKKOS=<prefix>` | the device backend (section 6.5); with the matching `CXX` |

**Which LAPACK.** The LAPACK LU (`linsol` LAPACK, or auto from 64
equations) is only as fast as the BLAS behind it. The integrator's own
dense LU (`linsol` dense, the dense LU of SUNDIALS, ported) uses no BLAS: plain
loops, the same everywhere.
- **Use an optimized BLAS:** OpenBLAS or MKL on Linux, Accelerate on macOS.
  One 655x655 factorization: the integrator's own LU 22.3 ms, OpenBLAS (1 thread)
  4.9 ms, Accelerate 1.7 ms.
- **Avoid the reference (Netlib) BLAS/LAPACK.** It is the default
  `-lblas -llapack` of many Linux systems (Debian/Ubuntu alternatives,
  the conda-forge `netlib` variant). It is unoptimized and can be slower
  than the integrator's own LU on large networks. A Linux test with it lost KROME-X's
  advantage over legacy; with OpenBLAS it recovered the speedups of the
  benchmark.
- **The library linked** through `LAPACK=` is the one called. `ldd` (Linux)
  or `otool -L` (macOS) on the executable shows which. With OpenBLAS, `OPENBLAS_VERBOSE=2` prints the
  core type it selected when it loads.
- **Threads:** many cells under OpenMP need `OPENBLAS_NUM_THREADS=1`
  (MKL: `MKL_NUM_THREADS=1`), so that each cell's factorization stays
  serial.
- **No rebuild needed to compare:** `linsol` is a runtime option (C key
  `linsol` 1 / 2 / 3, Python `linsol="dense"` / `"lapack"` / `"sparse"`).

**The sparse LU (`linsol` sparse).** Large chemical networks have sparse
Newton matrices. The generator orders the species to limit the fill
(minimum degree) and writes the rows of L + U with their fill
(`sp_*` tables in the Jacobian header; the header's comment gives the fill
and the flops of one factorization). The integrator factorizes them with
static pivots, as legacy KROME's DLSODES did. Tgas and the other non-species
equations form a dense border, eliminated through a small dense Schur
complement; a zero pivot is a recoverable failure (a smaller step). It
needs no library, and it runs in the device's per-thread kernel.
- **auto** takes it when one factorization needs at most a tenth of the
  dense flops, the border has at most 4 equations and the pattern is
  exact (`sparse_auto_`). The fill decides, not the size. One core, fixed T:
  DarkCloud (Semenov, 655 species, 11% filled) 2.2 s against 4.8 s with
  LAPACK (Accelerate) and 35 s with the dense LU; cloud (452) 0.37 s against
  0.69 s; ions (55, 8%) 3x faster; the 9-34-species benchmark networks,
  45-87% filled, are as fast or up to 1.5x slower, and stay dense.
- **Exactness:** the pattern holds every entry of the analytic Jacobian
  except those added by difference quotients (photo shielding, J_LW rates,
  secondary ionization, dust exchange), which the sparse LU drops: an
  approximate Newton matrix (the same solution, possibly more Newton
  iterations). auto avoids it there; `linsol` sparse forces it.
- **Jacobian:** it works with the analytic and the FD Jacobian alike (it
  reads the saved Jacobian; the FD one still costs one rhs call per
  column).

## 3. Generation options

The options of a configuration (section 2), by section: in a TOML file,
`key = value` under `[section]`; in Python, `kp.section.key = value`;
with `--set`, `section.key=value`. `options_template.toml` lists them
all, with their defaults. None of them can be changed at runtime (the
`[solver]` ones are only the defaults the runtime starts from).

### [project]

| key | meaning |
|---|---|
| `network` | the `.knet` network (or a legacy `react_*` file): a path, or the name of a shipped network (`"primordial"`, from `networks/`); required |
| `outdir` | the project directory (default `generated`, or `build_NAME` for a shipped test) |
| `bindings` | the optional bindings over the C interface: `["fortran"]` (`krome.f90`, `krome_net.f90`, `example.f90`, `make fortran`), `["python"]` (`krome.py`, `example.py`, `make shared`), or both; default none |
| `single` | `true`: the network as one `net_<name>.hpp` instead of the split files |
| `runtime` | `false`: only the generated network files, no runtime copy, Makefile or examples (default `true`) |

### [solver]

| key | meaning |
|---|---|
| `rtol` | the project's default relative tolerance (`SolverOpts::rtol` starts from it; default 1e-8; legacy `-RTOL`) |
| `atol` | the project's default absolute tolerance of the species [cm^-3] (`SolverOpts::atol_n`; default 1e-20; legacy `-ATOL`) |
| `atol_file` | per-species default absolute tolerances: `SPECIES ATOL` per line (the species by name or legacy index name, `idx_Hj`; Fortran reals; `#` comments; species not in the network are skipped with a warning); the other species keep `atol`. The starting `SolverOpts::atol_species` (legacy `-customATOL`) |
| `rtol_file` | the same for the relative tolerances (`SolverOpts::rtol_species`, legacy `-customRTOL`) |

These are defaults only: the caller can still change both at runtime.

### [rates]

| key | meaning |
|---|---|
| `outside` | a rate outside its declared T range: `"asym"` (default) zero below and clamped above; `"clamp"` both ends; `"zero"` both ends (legacy) |
| `rhs` | the mass-action rhs, fluxes and analytic Jacobian as generated statements (`"unrolled"`) or as loops over index tables (`"loop"`, legacy `-iRHS`): the same operations in the same order, much smaller code (large networks compile in a fraction of the time, and fit a GPU's instruction cache better). `"auto"` (default): loops above 500 reactions |
| `reactant_multiplicity` | `true`: 1/c! for c identical reactants (nuclear networks; legacy `-nuclearMult`) |
| `screening` | `true`: electron screening of nuclear rates (legacy `stars_screen`) |
| `check_reverse` | `true`: write `reverse_check.txt`, each reaction's reverse (or its absence) and dH(300 K) |

**Cosmological recombination.** The rate functions `recfast_H_rec(n(X))`,
`recfast_H_ion(n(X))`, `recfast_He_rec(n(He), n(H))` and
`recfast_He_ion(n(He), n(H))` are the H (and D) and He I recombination of
RECFAST 1.5.2 (Seager, Sasselov & Scott 1999, 2000; Wong, Moss & Scott
2008; Rubino-Martin et al. 2010) as rate coefficients: the case-B
recombination and the CMB photoionization from n = 2, each times the
Peebles / Sobolev escape factor, with RECFAST's fudge factors and the He I
triplets. They read `zredshift`, `hubble0` and `omega0`, with H(z) =
hubble0 (1 + z) sqrt(1 + omega0 z), the matter-dominated form of the
`EXPANSION` cooling. `networks/earlyUniverse.knet` uses them
(the shipped test `earlyUniverse`).

### [thermal]

`cooling` and `heating` are lists of process names (the legacy names),
e.g. `cooling = ["ATOMIC", "H2"]`. Thermal evolution as a whole is
switched at runtime (`evolve_T`). The processes themselves are fixed at
generation.

| cooling | process |
|---|---|
| `ATOMIC` | H and He atomic cooling: collisional excitation and ionization (Cen 1992), and the recombination cooling of the network's own recombination rates, n_e n_X kB T (3/2 k + T dk/dT) |
| `H2` | H2 rovibrational line cooling |
| `HD` | HD line cooling (uses Tcmb) |
| `CIE` | H2 collision-induced emission (uses Tcmb) |
| `CONT` | continuum cooling of primordial gas |
| `COMPTON` | Compton cooling on the CMB (uses zredshift) |
| `FF` | free-free (bremsstrahlung) |
| `EXPANSION` | cosmological expansion (REQUIRES hubble0, omega0) |
| `CHEM`, `DISS` | chemical cooling: the energy of dissociation reactions (the formation heating is the `CHEM` heating) |
| `Z` | metal fine-structure lines from level populations (coolZ.dat: C, Si, Fe, O, S and the ions C+, O+, Si+, Fe+; select with `metals`; C, C+ and O with current collision data from LAMDA, Lique et al. 2018, Klos et al. 2021, Yan & Babb 2023, including H2 and He as colliders; Si + H and S I from Yan & Babb 2023) |
| `Z_CIE`, `Z_CIENOUV` | Shen+13 CIE metal cooling tables (with and without UV background), always with their heating; scaled by `total_Z` |
| `Z_EXTENDED` | `Z` below 1e4 K merged with `Z_CIE` above |
| `CO`, `OH`, `H2O`, `HCN` | Neufeld & Kaufman molecular cooling tables (columns from `radiation.num2col`) |
| `DUST` | gas-dust thermal exchange with the effective grain population (d2g, Tdust) |
| `DUSTTAB` | tabulated dust cooling (`dust.tables`) |
| `DUSTBINS` | gas-dust exchange of the binned dust (section 9) |
| `DH` | reaction enthalpies at Tgas from the thermo30 NASA polynomials (alternative to CHEM) |

| heating | process |
|---|---|
| `COMPRESS` | compressional heating n kB T / tff (REQUIRES tff) |
| `VISCOUS` | viscous heating (REQUIRES visc_nu, visc_omega) |
| `PHOTO` | photoheating by the photon bins |
| `CR` | cosmic-ray heating |
| `CHEM` | chemical heating (H2 formation, including on dust; Solomon pumping) |
| `DH` | reaction enthalpies (with the `DH` cooling) |
| `PHOTOAV` | H2 photodissociation + UV pumping heating (not with CHEM) |
| `PHOTODUST` | dust photoelectric heating, Bakes & Tielens 1994 (G0, Av), or Weingartner & Draine 2001 with `photoelectric = "WD01"` |
| `PHOTODUSTNET` | the same, net of grain recombination cooling (G0, Av) |
| `NUCLEAR` | nuclear energy release, sum of flux x qeff |

| key | meaning |
|---|---|
| `cooling` | the cooling processes of the first table above (and the coolants of `coolant_files`) |
| `heating` | the heating processes of the second table above |
| `gamma` | adiabatic index: `"CONST"` = the parameter `gamma_ad` (default); `"EXACT"` = from composition and Tgas, every diatomic with its exact rotation and vibration (Omukai & Nishi 1998, Boley+07); `"POPOVAS"` = H2 and CO from the Popovas & Jørgensen (2016) ro-vibrational partition function tables. With H2 both use `opratio`, or the network's `H2_ortho` / `H2_para` abundances when it has them. The legacy approximations FULL, VIB, ROT and REDUCED were removed in 2.0.0 |
| `photoelectric` | the model of `PHOTODUST` / `PHOTODUSTNET`: `"BT94"`, Bakes & Tielens 1994 (default, legacy), or `"WD01"`, Weingartner & Draine 2001 eq. 44 (heating) and eq. 45 (grain recombination cooling), scaled by d2g / 0.01 |
| `wd01` | the row of WD01 Tables 2-3, `[RV, bC, case, field]` (default `["3.1", "6.0", "A", "ISRF"]`: Milky Way dust in the solar-neighbourhood ISRF; `field` B0 = a 3e4 K blackbody; `bC` the carbon in the smallest grains, 1e-5 per H) |
| `floor` | the coolers evaluated as cool(Tgas) - cool(Tfloor) (parameter `Tfloor`): `true`, all configured line coolers; a list, only those (legacy `-useIndividualFloor`); default `false` |
| `h2_opacity` | optically thick suppression of the H2 line cooling: `"RIPAMONTI"` or `"OMUKAI"` |
| `metals` | the metals of `Z` / `Z_EXTENDED` (C, Si, Fe, O, S, C+, O+, Si+, Fe+; default all present) |
| `photo_pumping` | `true`: radiative pumping of the metal-line levels by the photon field (legacy `-usePhotoInduced`) |
| `cmb_pumping` | `true`: the CMB (Tcmb) in the metal-line level balance |
| `coolant_files` | coolant files in the coolZ.dat format (legacy `-coolFile`: LAMDA conversions by `tools/lamda2krome.py`, e.g. `data/cooling/coolO2.dat`; CHIANTI conversions by `tools/chianti2krome.py`), added to the embedded coolZ.dat. Select their coolants by name in `cooling` (the legacy names: Roman numerals for atoms and ions, e.g. `FeII`; the species for molecules, e.g. `O2`), with `FILE` (all), or with `metals` (a metal a file defines comes from the file). Each is its own cooling term (`get_thermal` names it), optically thin, with the pumping and floor options of `Z` |
| `levels` | keep levels 0..L of every metal / file coolant, with their lines (legacy `-coolLevels`) |
| `escape` | the optical depth of the line coolants (the coolZ.dat metals and the `coolant_files` coolants): `"LVG"`, each line's radiative rates times its Sobolev escape probability (below); default optically thin |

Shipped coolant files (data/cooling/):
- `coolO2.dat`: O2 from LAMDA (36 of its 48 levels reachable).
- `coolFeII_16.dat`: Fe+ on its 16 levels below 14000 cm^-1 (a6D, a4F,
  a4D, a4P: the 26, 35, 51, 87 um, 5.3, 17.9 um and 1.26, 1.64 um lines),
  instead of the 5 a6D levels of coolZ.dat:
  `coolant_files = ["data/cooling/coolFeII_16.dat"]` with the `Z` cooling
  (an ion a file defines replaces its coolZ.dat block). Same sources as the coolZ.dat Fe+
  (CHIANTI electron rates, NIST A-values; H collisions within a6D only).
  It matters above ~2000 K (x1.26 at 4000 K in a test cell) and costs
  ~2.1 us per call against 0.17 us.

Converters (both write the coolZ.dat format):
- `python3 tools/lamda2krome.py o2.dat coolO2.dat`: a LAMDA molecular data
  file.
- `python3 tools/chianti2krome.py $XUVTOP coolC.dat c o_3 [--max-level 20]
  [--tmin 10 --tmax 1e9 --per-decade 4]`: ions (`o_3` = O++) or whole
  elements (`c`) of a CHIANTI v8+ database (the directory with the element
  directories). It writes the levels (`.elvlc`), the lines (`.wgfa`), and
  the electron (`.scups`) and proton (`.psplups`) de-excitation rates,
  descaled as CHIANTI does (Burgess & Tully 1992, natural cubic spline)
  and tabulated as log10 q on a log10 T grid. The proton data are
  excitation rate coefficients (the files' units): their de-excitation
  rates follow by detailed balance, held constant outside the fits'
  stated validity range. Legacy `tools/chianti.py` (Python 2) read the
  pre-v8 formats and converted the proton rates as collision strengths
  (about 1e-7 of their value).

**Optically thick lines** (`escape = "LVG"`, new in KROME-X: legacy
declared escape-probability arrays that nothing used). Every radiative term
of a line (A, and the pumping of `photo_pumping` / `cmb_pumping`) is
multiplied by its escape probability beta = (1 - e^-tau) / tau, in the
level balance and in the cooling (trapped photons are absorbed again
locally). The optical depth is the Sobolev one,
tau = A lambda^3 g_u / (8 pi g_l) (N / (b n)) (n_l - n_u g_l / g_u), with
N the species column of the `radiation.num2col` estimator (e.g. `"length"`
with the per-cell parameter `Lcol`: N = n Lcol, the Sobolev length) and b =
sqrt(2 kB T / m + vturb^2), `vturb` a per-cell parameter [cm s^-1], default
0. Masers (tau <= 0) keep beta = 1. The level solve is iterated, from the
optically thin solution as in RADEX: beta from under-relaxed populations
(0.3 old + 0.7 new after 5 plain steps) until
they change by less than 1e-10 of the species density, at most 100
iterations (realistic columns converge in 10-30). The populations queries
and the cooling terms include it.

Coolant-file physics, and how it differs from legacy:
- **Line energy:** the cooling of every line uses the level energies. Legacy
  used hc/lambda when a wavelength was given, and 0 for a non-positive one.
  In `data/cooling/coolChianti.dat` (legacy `tools/coolChianti.dat`) that silenced the 2s two-photon decays, which
  do radiate the level energy, as a continuum.
- **Levels:** the levels no transition connects to the ground are dropped
  (zero population), and all the others are solved. Legacy cut every level
  above the last one excited directly from the ground, including levels fed
  by cascades.
- **Rates:** evaluated directly at Tgas. Legacy interpolated them linearly
  on a 100-point log-T table, which overestimates the steep low-T rates.

### [radiation]

| key | meaning |
|---|---|
| `evolve` | `true`: co-evolve the photon field with the chemistry (as ln Phi, written back); default: the caller's field is held fixed |
| `shielding` | H2 self-shielding model of `krome_fshield()`: `"WG11"`, Wolcott-Green+ 2011 (default since 2.0.0), `"DB96"`, Draine & Bertoldi 1996, `"R14"`, Richings+ 2014, or `"NONE"` (f = 1, legacy's default); the H2 column is that of `H2`, or `H2_ortho` + `H2_para` |
| `habing_dust_shielding` | `true`: dust attenuation of the Habing flux from the photon bins |
| `secondary_ionization` | secondary ionization by fast photoelectrons on these species, e.g. `["H", "H2", "He"]` (Quadri+25) |
| `num2col` | the column-density estimator N(n) of the molecular cooling, the H2 opacity, the line-coolant escape probability, the shielding and the opacity query: `"default"` 1.87e21 (n / 1e3)^(2/3) [cm^-2], `"jeans"` / `"jeans40"` half a Jeans length (Tgas capped at 40 K), `"length"` n `Lcol` with the per-cell parameter `Lcol` [cm], required (e.g. the cell size, or a Sobolev or shielding length computed by the caller) |

### [dust]

| key | meaning |
|---|---|
| `tables` | embed the legacy dust tables (Tdust, dust cooling, `krate_H2dust_table()`): `"HM2012"` or `"MC"`, 2D for that field; `"3D"` with the per-cell Av; `"PHOTO3D"` with Av from the photon bins |
| `evolve_Td` | `true`: per-bin dust temperatures as solver state (needs `DUSTBINS`); default: radiative-balance bisection |
| `evolve` | the grain radii of the bins as solver state (section 9), by the processes listed: `["GROWTH", "SPUTTER", "EVAP"]`, plus `"KELVIN"` for the curvature correction to EVAP |
| `gas_total` | `true`: each gas/surface pair's surface slot holds gas + surface (a stability option) |

### [surface]

Grain-surface chemistry.

| key | meaning |
|---|---|
| `binding_energies` | binding energies: `"KROME"`, `data/Ebare_ice.dat` (E_diff = 2/3 E_bind), or `"SEMENOV"`, the OSU set of Semenov+10 (E_diff = 0.77 E_bind) |
| `nu` | attempt frequency: `"HASEGAWA"` (default) or `"LEGACY"` |
| `competition` | `false`: two-body surface rates without reaction-diffusion competition (legacy); default `true` |
| `tunnelling` | `false`: two-body surface reactions with a barrier by thermal hopping only, exp(-Ea/Tdust) (legacy Semenov framework); default `true`, the larger of that and tunnelling through a 1 A barrier |
| `reactive_desorption` | gas branch of two-body surface reactions: `"NONE"`, `"FIXED"` or `"MINISSALE"` |
| `bins` | per-bin surface species for NB dust bins (runtime `DustBins` with NB bins); default 0 |
| `cr_desorption` | cosmic-ray spot-heating desorption driven by the CR-rate parameter named here (e.g. `"crate"`) |

The grain population itself is part of the network. By default it is the
MRN distribution holding rho·d2g, or a runtime `DustBins` population. A
network can declare one grain size instead with `surface_grains single(n
= ..., a2 = ...)`, the legacy Semenov framework's grains. A reaction can
give its own sticking with `surface(stick=...)`. See
[network-format.md](network-format.md) section 5.3.

## 4. Runtime parameters

Each cell carries its own parameters:
- C++: the `Params` struct (`p.tff`);
- C and Fortran: a `par` row, indices `KROME_PAR_<name>` /
  `krome_par_<name>`;
- Python: `par={"tff": ...}`.

A network has only the parameters it uses. `NETWORK.md` lists them, with
their defaults. A parameter marked REQUIRED has no default: the physics
chosen at generation needs it. Until it is set, the solver rejects the
cell with `FLAG_PARAM_UNSET` (-105).

| parameter | default | meaning |
|---|---|---|
| `Tcmb` | 2.725 | CMB temperature [K]; also the rates' lower temperature cap: rates are evaluated at min(max(Tgas, Tcmb), 1e9 K), the heating and cooling at min(max(Tgas, 2.73 K), 1e9 K) (the legacy caps; below Tcmb the cooling keeps falling with Tgas, and `thermal.floor` stops the radiative cooling at Tfloor) |
| `gamma_ad` | 5/3 | adiabatic index (`thermal.gamma = "CONST"`) |
| `tff` | REQUIRED | compression timescale [s] (`thermal.heating = ["COMPRESS"]`) |
| `zredshift` | 0 | redshift |
| `hubble0` | REQUIRED | Hubble constant H0 [s^-1] (`EXPANSION`, `elec_recomb_ST93()`, the `recfast_*()` rates) |
| `omega0` | REQUIRED | matter density parameter (same) |
| `visc_nu` | REQUIRED | kinematic viscosity [cm^2 s^-1] (`thermal.heating = ["VISCOUS"]`) |
| `visc_omega` | REQUIRED | orbital frequency [s^-1] of the sheared flow (same) |
| `Tdust` | 10 | dust temperature [K] |
| `G0` | 0 | local radiation field [Habing] |
| `Av` | 0 | visual extinction [mag] |
| `d2g` | 0.01 | dust-to-gas mass ratio |
| `opratio` | 3 | H2 ortho:para ratio, splitting total `H2` where a process needs ortho and para (metal-line colliders, `thermal.gamma = "EXACT"` / `"POPOVAS"`); absent when the network has `H2_ortho` / `H2_para` (their abundances are used) |
| `total_Z` | 1 | metallicity [Z_sun] of the ZCIE tables |
| `clump` | 1 | clumping factor of the Jura H2-on-dust rate |
| `Tfloor` | 2.725 | cooling-floor temperature [K] (`thermal.floor`) |

A network can also declare its own parameters, such as `crate` for the
cosmic-ray rate. They appear in `NETWORK.md` as "network parameter".

## 5. Solver options

These are runtime numerics.
- C++: fields of `SolverOpts` (and `CellsOpts`, `EquilOpts`).
- C: `krome_opts_set(h, key, value)` on a handle from
  `krome_opts_new()`.
- Fortran: `krome_set(h, key, value)`.
- Python: `net.options(key=value, ...)`, or `opts[key] = value`.

| C++ field | key | default | meaning |
|---|---|---|---|
| `SolverOpts::evolve_T` | `evolve_T` | false | evolve Tgas with the chemistry (the configured thermal processes); false: fixed Tgas |
| `SolverOpts::jac` | `jac` | `JAC_AUTO` (0) | Jacobian: `JAC_AUTO` = analytic at fixed Tgas, difference quotients with Tgas evolved; `JAC_ANALYTIC` (1): the analytic chemistry Jacobian, with Tgas evolved the Tgas column and the thermal rows by difference quotients (faster for large networks, about 1.6 times on 655 species, but less robust at very small atol or with evolving dust); `JAC_FD` (2): difference quotients throughout (legacy `MF=222`) |
| `SolverOpts::linsol` | `linsol` | `LINSOL_AUTO` (0) | LU of the Newton systems: auto (the sparse LU where it saves 10x the flops, section 2; else LAPACK from `LAPACK_MIN_NEQ` = 64 equations when built with it), `LINSOL_DENSE` (1, the integrator's own: the dense LU of SUNDIALS, ported), `LINSOL_LAPACK` (2), `LINSOL_SPARSE` (3, static pivots in the generated order) |
| `SolverOpts::rtol` | `rtol` | 1e-8, or codegen `solver.rtol` | relative tolerance |
| `SolverOpts::atol_n` | `atol_n` | 1e-20, or codegen `solver.atol` | absolute tolerance of the species [cm^-3] |
| `SolverOpts::atol_species`, `rtol_species` | C `krome_opts_set_atol_species`, `krome_opts_set_rtol_species` (getters `krome_opts_get_atol_species`, `krome_opts_get_rtol_species`); Fortran `krome_set_atol_species`, `krome_set_rtol_species` (`krome_get_atol_species`, `krome_get_rtol_species`); Python `opts.set_atol_species`, `set_rtol_species` (an array or a dict `{species: value}`; properties `opts.atol_species`, `opts.rtol_species`) | none, or codegen `solver.atol_file` / `solver.rtol_file` | per-species tolerances: NSPEC values replacing `atol_n` / `rtol` for the species; a negative entry (and every species without an array) takes the scalar. The C++ arrays are the caller's; the C, Fortran and Python setters copy them into the handle, NULL / no argument / `None` clears them, and the getters return the tolerance each species gets. With a per-species rtol the error weights use a per-slot rtol (CVODE's formula), on the host and on the device |
| `SolverOpts::atol_T` | `atol_T` | 1e-6 | absolute tolerance of Tgas and dust temperatures [K] |
| `SolverOpts::atol_u` | `atol_u` | 1e-10 | absolute tolerance of ln Phi (`radiation.evolve`) |
| `SolverOpts::max_steps` | `max_steps` | 100000 | step budget; exhausting it is final (`FLAG_TOO_MUCH_WORK`) |
| `SolverOpts::max_err_fails` | `max_err_fails` | 7 | error-test failures per step |
| `SolverOpts::max_conv_fails` | `max_conv_fails` | 10 | Newton convergence failures per step |
| `SolverOpts::max_restarts` | `max_restarts` | 5 | restarts at the reached state after error-test or convergence failures, each from a tenth of the step size the integrator had reduced to |
| `SolverOpts::clip_negative` | `clip_negative` | true | n = max(n, 0) after a solve |
| `SolverOpts::renormalize_x` | `renormalize_x` | false | mass-fraction solves (`solve_x`, `krome_solve_cells_x`, `net.solve_x`): rescale the returned fractions to their input sum (legacy `-useX` always did) |
| `SolverOpts::conserve` | `conserve` | 0 | 0 off; 1 electrons restore the input charge; 2 the minimal correction restoring the input element totals and charge; 3 the electrons neutralize the ions, also while a `freeze` holds e- fixed (legacy `-computeElectrons`) |
| `SolverOpts::error_report` | `error_report` | true | write the error report on a failure (section 10) |
| `SolverOpts::error_file` | `krome_opts_set_error_file` / `krome_set_error_file` / `opts["error_file"]` | `KROME_ERROR_REPORT` | the report file |
| `SolverOpts::photo`, `SolverOpts::Phi` | `krome_opts_set_photo`, the `Phi` argument | none | photon bins and field (section 8) |
| `SolverOpts::dust` | `krome_opts_set_dustbins` | none | binned dust (section 9) |
| `CellsOpts::path` | `path` | `PATH_AUTO` (0) | many cells on the CPU or the device (section 6.5): `PATH_AUTO`, `PATH_PERCELL` (1, the CPU), `PATH_DEVICE` (2, the Kokkos backend) |
| `CellsOpts::device_cache` | `device_cache` | true | keep the device buffers between calls (sized for the largest call so far; memory only, every call writes its own inputs); `kk::release_device_cache()` frees them, `Kokkos::finalize` does |
| `CellsOpts::max_device_cells` | `max_device_cells` | 0 | the device path integrates at most this many cells at a time (larger calls in chunks on the same buffers); 0: all at once |
| `CellsOpts::min_interleave` | `min_interleave` | -1 | the device path's memory layout: interleaved across cells (coalesced on a GPU) for calls of at least this many cells, contiguous per cell below; -1: contiguous, and on CUDA interleaved above the cells the device holds resident at once (SMs x 256; 27 648 on an A100), with the launch bounds measured best for each (contiguous 32,1, interleaved 128,8). The results agree to roundoff (bitwise when built without FMA contraction) |
| `CellsOpts::min_device` | `min_device` | 10000 | `PATH_AUTO`: calls from this many cells go to the device, when the build has one (section 6.5) |
| `CellsOpts::Phi` | `krome_cells::Phi` (C), `phi=` / `Phi=` | none | per-cell fields `Phi[ncell][nbins]`, read (and with `radiation.evolve` updated). Without them every cell starts from `SolverOpts::Phi`; with `radiation.evolve` a single-cell call updates it, a call of several cells does not return the evolved fields |
| `CellsOpts::flags` | the `flags` argument | none | per-cell return flags |
| `CellsOpts::stats` | (C++ only) | none | per-cell `SolveStats[ncell]` |
| `CellsOpts::Td` | `krome_cells::Td` (C), Fortran `td=`, Python `Td=` | none | per-cell dust temperatures `Td[ncell][ndust]`: each cell starts from its row and the row receives its final Td (binned dust) |
| `CellsOpts::adust` | `krome_cells::adust` (C), `adust=` | none | per-cell grain radii `adust[ncell][ndust]` (`dust.evolve`) |
| `EquilOpts::t0` | `eq_t0` | 100 yr | first output time of the equilibrium solver [s] |
| `EquilOpts::t_max` | `eq_tmax` | 5e8 yr | give up after this time [s] |
| `EquilOpts::factor` | `eq_factor` | 3 | output times t0 x factor^k |
| `EquilOpts::ertol` | `eq_rtol` | 1e-5 | converged when the relative change between outputs is below this |
| `EquilOpts::eatol` | `eq_atol` | 1e-12 | ... with n floored at eatol x ntot |

**Charge in the Newton iteration.** The integrator is CVODE's BDF with
one addition: the Newton iteration keeps the net charge exactly. An
exact Jacobian conserves it, but difference quotients (`JAC_AUTO` with
Tgas evolved, `JAC_FD`, and the difference-quotient entries of
`JAC_ANALYTIC`) break it at roundoff over tiny increments, and the drift,
conserved from then on, can outgrow the electrons of a collapsing cloud.
The electron component of each Newton correction is therefore set so
that the correction conserves the charge; the Jacobian is untouched, so
the solution is the same. The correction stays the solver's own where
the electrons carry less than 1e-6 of the summed |charge x density| of
the other species, or would become negative (an input with more net
charge than electrons). It is always on, with every Jacobian, in a
network with e- and no freezable species, and off with `dust.gas_total`;
`conserve` (below) instead corrects the state after a solve.

Python also accepts names for the enumerated keys:
- `jac`: `"auto"`, `"analytic"`, `"fd"`;
- `linsol`: `"auto"`, `"dense"`, `"lapack"`, `"sparse"`;
- `path`: `"auto"`, `"percell"`, `"device"`.

A key that does not exist is rejected: `krome_opts_set` returns -1, and
Python raises. That includes old process switches, because the physics
is fixed at generation.

## 6. Calling KROME

### 6.1 C++

```cpp
#include "krome.hpp"
using namespace krome;

double n[NSPEC] = {};                  // [cm^-3], order: NETWORK.md
n[IDX_H] = 1e4; n[IDX_Hj] = n[IDX_E] = 1.0; n[IDX_He] = 775.0;
double Tgas = 100.0;
Params p;  p.tff = 1e14;               // REQUIRED by thermal.heating = ["COMPRESS"]
SolverOpts o;  o.evolve_T = true;
int flag = solve(n, Tgas, p, 1e10, o); // 0: n, Tgas updated
if (flag != FLAG_OK) printf("%s\n", flag_message(flag));
```

Link the network data (no external library; add `-DKROME_USE_LAPACK` and
a LAPACK for the LAPACK LU):

```sh
c++ -O2 -std=c++17 -I myproj my.cpp myproj/net_<name>_data.o
```

| function | purpose |
|---|---|
| `solve(n, Tgas, p, t_end, o, stats)` | one cell, from 0 to t_end; `SolveStats` optional |
| `solve_x(x, rho, Tgas, p, t_end, o, stats)` | the same with mass fractions |
| `solve_cells(ncell, n, Tgas, dt, p or params, o, co)` | many cells, each over its own dt |
| `solve_cells_x(ncell, x, rho, Tgas, dt, p or params, o, co)` | the same with mass fractions |
| `solve_equilibrium(n, Tgas, p, o, e, t_reached)` | to steady state |
| `flag_message(flag)` | a flag's description |

`krome.hpp` can be included in any number of translation units. Every
declaration it lists is the supported interface. Everything in `core/`
is internal.

### 6.2 C

```c
#include "krome.h"
#include "krome_net.h"
double par[KROME_NPAR];
krome_params_default(par);            /* REQUIRED ones are NaN */
par[KROME_PAR_tff] = 1e14;
krome_opts* h = krome_opts_new();
krome_opts_set(h, "evolve_T", 1);
int r = krome_solve_cells(ncell, n, Tgas, dt, par, 0 /* one row */,
                          NULL /* no per-cell extras */, flags, h);
if (r) printf("%s\n", krome_flag_message(r));
krome_opts_free(h);
```

The optional per-cell arrays go in a `krome_cells` struct (NULL: none).
Its members are `Phi` (photon bins), `Td` (dust bins) and `adust`
(`dust.evolve`); unused members stay NULL, and new per-cell fields are
added there:

```c
krome_cells ex = {0};
ex.Phi = Phi;    /* [ncell][nbins] */
krome_solve_cells(ncell, n, Tgas, dt, par, 0, &ex, flags, h);
```

Link `libkrome.a` and the C++ runtime (`-lc++` or
`-lstdc++`).

Array layouts:

| array | layout |
|---|---|
| species | `n[ncell][KROME_NSPEC]` |
| parameters | `par[ncell][KROME_NPAR]` (`par_stride` = KROME_NPAR), or one row for every cell (`par_stride` 0) |
| photon field | `Phi[ncell][nbins]` |
| flags | `flags[ncell]` (may be NULL) |

The C functions:

| group | functions |
|---|---|
| solving | `krome_solve_cells`, `krome_solve_cells_x`, `krome_equilibrium_cells`, `krome_flag_message` |
| options | `krome_opts_new`, `krome_opts_free`, `krome_opts_set`, `krome_opts_get`, `krome_opts_set_error_file`, `krome_opts_set_photo`, `krome_opts_set_photo_broadening`, `krome_opts_nbins`, `krome_opts_set_dustbins`, `krome_opts_ndust`, `krome_opts_set_dust_Td`, `krome_opts_get_dust_Td` |
| device | `krome_has_device`, `krome_device_initialize`, `krome_device_finalize` |
| network | `krome_version`, `krome_network_name`, `krome_nspec`, `krome_nrea`, `krome_npar`, `krome_species_name`, `krome_species_index`, `krome_param_name`, `krome_param_index`, `krome_params_default`, `krome_param_required_by`, `krome_species_mass`, `krome_species_charge`, `krome_nelem`, `krome_element_name`, `krome_species_natoms`, `krome_reaction_name` |
| queries | section 7 |
| radiation | section 8 |

### 6.3 Fortran

```fortran
use krome
use krome_net
real(8) :: n(krome_nspec, ncell), Tgas(ncell), dt(ncell)
real(8) :: par(krome_npar, ncell)
integer(c_int) :: flags(ncell)
type(c_ptr) :: h
call krome_params_default(par)
par(krome_par_tff, :) = 1d14
h = krome_opts_new()
call krome_set(h, "evolve_T", 1d0)
ierr = krome_solve(n, Tgas, dt, h, par=par, flags=flags)
if (ierr /= krome_flag_ok) print *, trim(krome_flag_message_f(ierr))
call krome_opts_free(h)
```

- **Layout:** the arrays are column-major, `n(nspec, ncell)`, the
  same memory layout as C.
- **Indices** are 1-based: `krome_idx_<sp>`, `krome_par_<name>`.
- **Parameters:** `par` can be one column `par(krome_npar, 1)` for
  every cell.
- **Names:** the Fortran names are the C names. The `krome_get_*`
  queries are functions that return `ierr` and fill output arrays, e.g.
  `ierr = krome_get_rates(n, Tgas, h, k, par)` with `k(nrea, ncell)`.
  Names that return strings end in `_f` (`krome_species_name_f`).
- **Mass fractions:** `krome_solve_x(x, rho, Tgas, dt, h, par, phi,
  flags)`.

The Fortran names that differ from C:

| purpose | Fortran |
|---|---|
| solve | `krome_solve`, `krome_solve_x`, `krome_equilibrium` |
| options | `krome_set(h, key, value)`, `krome_get(h, key)`, `krome_set_error_file`, `krome_set_photo`, `krome_set_photo_broadening`, `krome_nbins`, `krome_set_dustbins`, `krome_ndust`, `krome_set_dust_Td`, `krome_get_dust_Td`, `krome_photo_sigma` |
| counts (the constants are in `krome_net`) | `krome_nspec_c`, `krome_nrea_c`, `krome_npar_c`, `krome_nelem_c`, `krome_nheat_c`, `krome_ncool_c`, `krome_nphoto_c` |
| index of a name (1-based) | `krome_species_idx`, `krome_param_idx` |
| names by 1-based index (and the version, `krome_version_f`) | `krome_species_name_f`, `krome_param_name_f`, `krome_param_required_by_f`, `krome_reaction_name_f`, `krome_element_name_f`, `krome_heating_name_f`, `krome_cooling_name_f`, `krome_flag_message_f` |
| species and photo metadata by 1-based index | `krome_species_mass_f`, `krome_species_charge_f`, `krome_species_natoms_f`, `krome_photo_reaction_f`, `krome_photo_reactant_f`, `krome_photo_xsec_f` |

### 6.4 Python

```python
import krome
net = krome.load()                       # libkrome next to krome.py
n = net.densities({"H": 1e4, "H+": 1.0, "e-": 1.0, "He": 775.0})
res = net.solve(n, Tgas=100.0, dt=1e10, par={"tff": 1e14},
                opts=net.options(evolve_T=True))
res.n, res.Tgas, res.flags               # KromeError on failure
                                         # (check=False to inspect)
```

- **Arrays:** `n` can be `[nspec]` or `[ncell, nspec]`. Tgas, dt and
  each parameter can be a scalar or one value per cell.
- **Solving:** `net.solve`, `net.solve_x`, `net.equilibrium`,
  `net.flag_message`.
- **Setup:** `net.idx`, `net.par_idx`, `net.params_default`,
  `net.param_rows`, `net.densities`, `net.options` (`Options.set`,
  `set_photo`, `set_dustbins`, `nbins`).
- **Device:** `net.has_device`, `net.device_initialize`,
  `net.device_finalize`.
- **Queries and radiation:** the methods listed in sections 7 and 8,
  plus `net.photo_xsec(r, E)` (the cross section of photo reaction r).

### 6.5 Many cells, OpenMP and the device

`solve_cells` / `krome_solve_cells` (and `net.solve` with several cells)
integrate independent cells, on the CPU or on a device.

**Which device: at build time.** A project has a device backend only when
it is built with Kokkos, `make KOKKOS=<prefix>` with the compiler Kokkos
was built with (it defines `KROME_USE_KOKKOS`). The device is the one of
that Kokkos installation:

| Kokkos built for | `CXX` | the device |
|---|---|---|
| CUDA | `$KOKKOS/bin/nvcc_wrapper` | NVIDIA GPUs |
| HIP | `hipcc` | AMD GPUs |
| OpenMP | a host compiler with `-fopenmp` | the CPU threads, one cell per thread |

e.g. `make lib KOKKOS=$HOME/kokkos CXX=$HOME/kokkos/bin/nvcc_wrapper`.
Without `KOKKOS=` the library has no device path, and every call runs on
the CPU (`krome_has_device()` / `net.has_device` says which).

**CPU or device: at run time, call by call.** Every many-cell call goes
through a dispatcher, set by `CellsOpts::path` (the key `path` in C,
Fortran and Python, with the names of section 5):

| `path` | the cells run |
|---|---|
| `PATH_AUTO` (0, default) | on the device when the build has one, Kokkos is initialized **and** the call has at least `min_device` cells (default 10 000); on the CPU otherwise |
| `PATH_PERCELL` (1) | on the CPU: a loop of independent `solve()` calls, parallel over the cells with OpenMP when built with `-fopenmp` |
| `PATH_DEVICE` (2) | on the device, whatever their number (`-104` in a build without one) |

The default threshold follows the measurements of the CPU and GPU
comparison: below a few thousand cells per call a GPU is slower than its
share of CPU cores (`min_device` tunes it). A single-cell `solve()`
always runs on the CPU. A configuration the device backend does not
support runs on the CPU, and cells that fail on the device are solved
again on the CPU, so only cells that really fail report a failure.

**Initialize the device once**, before the first device call:
`Kokkos::initialize` in C++ (or a Kokkos-based host code's own),
`krome_device_initialize()` in C and Fortran, `net.device_initialize()`
in Python; and finalize it before exiting. Before that, `PATH_AUTO`
runs every call on the CPU, and an explicit `PATH_DEVICE` returns
`FLAG_NO_KOKKOS` (-109).

C++:

```cpp
#include <Kokkos_Core.hpp>
#include "krome.hpp"
using namespace krome;

Kokkos::initialize(argc, argv);
{
  CellsOpts co;                 // co.path = PATH_AUTO: the device from
  co.min_device = 5000;         // 5000 cells (default 10000)
  std::vector<int> flags(ncell);
  co.flags = flags.data();
  int f = solve_cells(ncell, n, Tgas, dt, params, o, co);
}
Kokkos::finalize();
```

(`solve_cells_kokkos`, in namespace `krome::kk`, calls the device
backend directly, with no CPU fallback.)

C:

```c
krome_device_initialize();
krome_opts* h = krome_opts_new();
krome_opts_set(h, "path", 2);         /* always the device; 0 = auto */
krome_opts_set(h, "min_device", 5000);   /* the auto threshold */
/* ... krome_solve_cells(..., h) ... */
krome_opts_free(h);
krome_device_finalize();
```

Fortran: the same keys, `call krome_set(h, "path", 2d0)`, and the
functions `krome_device_initialize()`, `krome_device_finalize()`,
`krome_has_device()`.

Python:

```python
net = kp.build()                       # a project built with KOKKOS=...
if net.has_device:
    net.device_initialize()
opts = net.options(path="auto", min_device=5000)
r = net.solve(cells, Tgas, dt, par=par, opts=opts)   # cells[ncell, nspec]
net.device_finalize()                  # before exiting
```

- **What the device supports:** every configuration of the host:
  fixed or evolved Tgas, static or co-evolved photon bins
  (`radiation.evolve`), binned dust (DUSTBINS cooling, the surface hook,
  `surface.bins`), `dust.evolve_Td`, `dust.evolve`,
  `radiation.secondary_ionization`, the dust tables. The per-cell rows
  (`CellsOpts::Phi`, `Td`, `adust`) are read and written as on the host.
- **How it integrates:** the cells of a call are copied to the device
  together, and one kernel runs one thread per cell. Each thread runs
  the host's integrator (`core/cvode/bdf_cell.hpp`) on the host's own rhs
  and Jacobian (the same functions, compiled for the device) to its own
  dt, with its own step sizes and order. The result is bitwise equal to
  the host `solve()` when both are compiled with the same floating-point
  contraction. Work memory per cell, plus the rhs work: with the dense
  (or LAPACK) LU, 15 NEQ + 2 NEQ^2 doubles (the saved Jacobian and the
  Newton matrix); with the sparse LU, 15 NEQ + NEQ^2 doubles plus the
  factors (the nonzeros of L + U, `SPARSE_NNZ` in the Jacobian header)
  and the dense border. The memory, not the speed, limits the cells of
  one call on a large network: DarkCloud (655 species, sparse LU) needs
  about 4 MB per cell, so an A100 (64 GB) holds about 1.5x10^4 cells, and
  `max_device_cells` splits larger calls into chunks. The interleaved kernel (large calls on a GPU) sizes its
  thread-local arrays at compile time for at most
  `KROME_DEVICE_MAX_PHOTO_BINS` co-evolved photon bins (default 32) and
  `KROME_DEVICE_MAX_DUST_BINS` dust bins (default 40); a call with more
  runs the contiguous kernel.
- **Failures:** only the cells that fail on the device are re-solved on
  the host, so only cells that really fail report.

### 6.6 Examples: the shipped tests and hydro coupling

**The shipped tests** are the legacy tests, ported to KROME-X (legacy
`./krome -test=NAME`):

    python3 -m kromex --test collapse          # into build_collapse/
    cd build_collapse && make && ./test        # then: python3 plot.py

or, in Python, `Preprocessor.from_test("collapse").run()`; the test's
options can be changed before running (`--set`, or the attributes).

- **What it writes:** the project with the test's options
  (`examples/tests/NAME/options.toml`, recorded in the project's
  `krome.toml`), the driver `test.cpp` (the C++ API) with the helper
  header `krome_test_utils.hpp`, the plot script `plot.py` (numpy +
  matplotlib, columns picked by name from the output header) and the
  test's data files. The driver is the project's program: `make` builds
  it, and no `example.cpp` is written.
- **The Fortran-interface test** (`--test fortran`): the collapse test
  through the Fortran binding (`test.f90`, helper module
  `krome_test_utils.f90`) and through the C++ API (`test_cpp.cpp`).
  `check_diffs.py` compares the two outputs.
- **Output:** the same quantities as legacy, in files named by their
  content (`abundances.dat`, `cooling.dat`, ...); each test's README
  gives the legacy name of each file.
- **Finding tests:** `python3 -m kromex --tests` (or
  `Preprocessor.tests()`) lists them, and each test's README notes what
  changed from legacy.

**`examples/coupling/`** holds four drivers coupling KROME to a hydro code,
one for each calling style:

- Fortran, RAMSES-like: passive scalars and `krome_solve_x`;
- C, ENZO-like: species fields in code units and `krome_solve_cells`;
- C++, GIZMO-like: particles and `solve_cells` with per-particle `Params`;
- Python: post-processing a snapshot.

All four solve the same operator-split problem: gather the cells, make
one call, then rebuild the internal energy from Tgas with `get_mu` and
`get_gamma`. A failed cell is subcycled. Their README covers units,
energy, failures and per-cell inputs.

## 7. Queries

Queries evaluate the cell as the solver sees it: the rates at the
temperature capped to [Tcmb, 1e9 K], the thermal terms and gamma capped to
[2.73 K, 1e9 K], the photon field and dust of the options. C takes the
`krome_solve_cells` array conventions.

| quantity | C++ | C / Fortran | Python |
|---|---|---|---|
| rate coefficients k[NREA] | `get_rates` | `krome_get_rates` | `rates` |
| fluxes k prod(n) [cm^-3 s^-1] | `get_fluxes` | `krome_get_fluxes` | `fluxes` |
| dust temperature of the rates (legacy `krome_get_table_Tdust`: the `dust.tables` value, else the `Tdust` parameter; -1 without either) | `get_Tdust` | `krome_get_Tdust` | `Tdust` |
| heating, cooling (totals and per term), dT/dt | `get_thermal`, `heating_name`, `cooling_name` | `krome_get_thermal`, `krome_nheat`, `krome_ncool`, `krome_heating_name`, `krome_cooling_name` | `thermal` |
| energy per reaction [MeV] | `get_qeff` | `krome_get_qeff` | `qeff` |
| adiabatic index | `get_gamma` | `krome_get_gamma` | `gamma` |
| mean molecular weight | `get_mu` | `krome_get_mu` | `mu` |
| mass density, total density, H nuclei | `get_rho`, `get_ntot`, `get_Hnuclei` | `krome_get_rho`, `krome_get_ntot`, `krome_get_Hnuclei` | `rho`, `ntot`, `hnuclei` |
| net charge; electrons for neutrality | `get_charge`, `get_electrons_neutral` | `krome_get_charge`, `krome_get_electrons` | `charge`, `electrons` |
| element totals | `get_element_totals` | `krome_get_element_totals` | `element_totals` |
| free-fall time, Jeans length | `get_free_fall_time`, `get_jeans_length` | `krome_get_free_fall_time`, `krome_get_jeans_length` | `free_fall_time`, `jeans_length` |
| free-fall time from a mass density | `get_free_fall_time_rho` | `krome_get_free_fall_time_rho` | `free_fall_time_rho` |
| column density of a number density (the `radiation.num2col` estimator) | `num2col` | `krome_num2col` | `num2col` |
| line coolants (the coolZ.dat metals, then the `thermal.coolant_files` coolants): count, legacy names (CII, FeII, O2), species, levels (energies [K], degeneracies, the file's numbering) | `NCOOLANT`, `coolant_name`, `coolant_species`, `coolant_nlev`, `coolant_level_E`, `coolant_level_g` (`coolant_lev0` offsets) | `krome_ncoolant`, `krome_coolant_name`, `krome_coolant_species`, `krome_coolant_nlev`, `krome_coolant_levels`; Fortran `krome_ncoolant_c`, `krome_coolant_name_f`, `krome_coolant_species_f`, `krome_coolant_nlev_f`, `krome_coolant_levels` | `coolants`, `coolant_species`, `coolant_levels` |
| each coolant's line cooling [erg cm^-3 s^-1] (legacy `krome_cooling<metal>`; the Z and coolant-file terms are their sums, before the `Z_EXTENDED` weight) | `get_coolant_cooling` | `krome_get_coolant_cooling` | `coolant_cooling` |
| a coolant's level populations [cm^-3] (legacy `krome_popcool_dump`; levels no transition connects to the ground are empty) | `get_coolant_populations` | `krome_get_coolant_populations` | `coolant_populations` |
| photoheating per photo reaction [erg cm^-3 s^-1], summing to the PHOTO term (with `radiation.secondary_ionization`, each reaction's heat share; legacy `krome_get_photoBin_heats` gave it per absorber) | `get_photo_heating` | `krome_get_photo_heating` | `photo_heating` |
| the field over the Draine ISRF (legacy `krome_get_ratioFluxDraine`) | `get_ratio_flux_draine` | `krome_get_ratio_flux_draine` | `ratio_flux_draine` |
| n <-> mass fractions | `n2x`, `x2n` | `krome_n2x`, `krome_x2n` | `n2x`, `x2n` |
| optical depth per bin | `get_opacity` | `krome_get_opacity` | `opacity` |
| photon removal per bin [cm^-1] | `get_absorption` | `krome_get_absorption` | `absorption` |
| Habing flux of a field | `PhotoModule::habing` | `krome_get_habing` | `habing` |

**Legacy conveniences** (legacy `krome_user`; one cell, the query
conventions for `par`, `Phi` and the options):

| legacy | C++ | C | Fortran | Python |
|---|---|---|---|---|
| `krome_get_names_header` | `names_header()` | `krome_names_header` | `krome_names_header()` | `names_header()` |
| `krome_print_best_flux` | `print_best_flux(n, T, p, nbest, o, file)` | `krome_print_best_flux` (stdout) | `krome_print_best_flux(n, T, nbest, h, par, phi)` | `print_best_flux`, `best_fluxes` (a list) |
| `krome_explore_flux` | `explore_flux(file, n, T, p, xvar, o)` | `krome_explore_flux(path, ...)` (appends) | `krome_explore_flux(n, T, unit, xvar, h, par, phi)` | `explore_flux(file or path, ...)` |
| `krome_get_info` | `print_info(n, T, file)` | `krome_print_info` (stdout) | `krome_print_info(n, T)` | `print_info` |
| `krome_scale_Z` | `scale_Z(n, Z)` (generated) | `krome_scale_Z(ncell, n, Z)` | `krome_scale_Z(n, Z)`, n(nspec) or n(nspec, ncell) | `scale_Z(n, Z)` (a copy) |
| `krome_popcool_dump` | `print_coolant_populations(file, n, T, p, o)` | `krome_print_coolant_populations(path, ...)` (appends) | `krome_print_coolant_populations(n, T, unit, h, par, phi)` | `print_coolant_populations(file or path, ...)` |
| `krome_dump_cooling` | `dump_cooling(file, n, T, p, o)` | `krome_dump_cooling(path, ...)` | `krome_dump_cooling(n, T, unit, h, par, phi)` | `dump_cooling(file or path, ...)` |
| `krome_dump_rates` | `dump_rates(file, Tmin, Tmax, npts, p, o)` | `krome_dump_rates(path, ...)` | `krome_dump_rates(Tmin, Tmax, npts, unit, h, par, phi)` | `dump_rates(file or path, ...)` |
| `krome_print_phys_variables` | `print_params(p, file)` | `krome_print_params(par)` (stdout) | `krome_print_params(par)` | `print_params(par, file)` |

The formats are legacy's: the populations one row per level (name, level
from 1, Tgas, population, the coolant's total, then a blank line), the
cooling one row (Tgas, total, every term), the rates rows (Tgas, reaction
from 1, k at n = 1e-40) with a blank line per temperature.

`scale_Z` sets the network's neutral atomic metals (species named as an
element, not H or He) to n_H 10^(Z + log10(n/n_H)_sun), with the Asplund+09
photospheric abundances of `data/asplund.dat` and Z = log10(Z/Z_sun),
floored at 1e-40, as legacy.

**Physical constants.** CODATA 2018 (the exact SI values of the 2019
redefinition where they apply) and the IAU astronomical units, from one
table, `kromex/constants.py`: the numbers the solver itself uses. C++
`krome::constants::p_mass` (defined by the generated network header), C
`KROME_CONST_p_mass`, Fortran `use krome_constants` then `p_mass`, Python
`krome.constants.p_mass`. The names are legacy's, plus a few:
`boltzmann_erg` / `_eV` / `_J`, `planck_erg` / `_eV` / `_J`, `hbar_erg`,
`hc_k` (hc/k_B [cm K]), `gravity`, `p_mass`, `e_mass`, `n_mass`, `amu`,
`clight`, `pi`, `eV_to_erg`, `ry_to_eV`, `ry_to_erg`, `seconds_per_year`
(Julian year), `km_to_cm`, `pc_to_cm`, `Mpc_to_cm` (legacy `cm_to_Mpc` =
1 / `Mpc_to_cm`), `stefboltz_erg`, `N_avogadro`, `Rgas_J`, `Rgas_kJ`.
Network expressions can name `boltzmann_erg`, `boltzmann_eV`, `Rgas_kJ`,
`N_avogadro`, `pi`, `planck_eV` and `eV_to_erg`. Legacy KROME used older values (CODATA 2006-2010,
pi = 3.14159265359, 1 Mpc = 3.08e24 cm, a 365-day year): the results move
by ~1e-6, except where the old values were wrong (hc/k in the `thermal.gamma`
temperatures, 0.7 %) or coarse (the Mpc, 0.2 %).

## 8. Photon bins and radiation

A network with photo reactions needs photon bins and a field. The field
is the caller's photon flux per bin, `Phi[nbins]` [cm^-2 s^-1], one
array per cell.
- **Fixed field (default):** the field is held constant during a
  solve.
- **With `radiation.evolve`:** the field is co-evolved and written back.
  The bins only absorb (d ln Phi_b/dt = -c kappa_b), so in a two-moment
  (M1) radiation code the flux of each bin is attenuated by the same
  factor as its photon density: after the call, multiply every flux
  component F_b by Phi_b(after) / Phi_b(before). This is exact, keeps the
  reduced flux |F_b| / (c N_b) unchanged, and needs no flux equations in
  the solver (the rates depend on the photon density only).

Setting up the bins:
- C++:
  ```cpp
  Radiation rad = Radiation::log_bins(B, emin, emax);
  PhotoModule ph;
  ph.init(rad);
  o.photo = &ph;
  o.Phi = phi;
  ```
- C: `krome_opts_set_photo(h, nbins, Elo, Ehi)`, then pass `Phi` to the
  solve.
- Python: `opts.set_photo(Elo, Ehi)`.

Doppler broadening of the tabulated (Leiden/SWRI) cross sections,
legacy `-useBroadening`, is applied once when the bins are built, as in
legacy. It uses a temperature and a turbulent velocity [cm s^-1]:
- C++: `rad.broad_T`, `rad.broad_vturb` before `ph.init(rad)`;
- C: `krome_opts_set_photo_broadening(h, T, vturb)`;
- Fortran: `krome_set_photo_broadening(h, T, vturb)`;
- Python: `opts.set_photo_broadening(T, vturb)`.

Each energy sees the table's mean over a box of width `E sqrt(2 kB T / m
+ vturb^2) / c`. Unlike legacy's widening, this conserves the integrated
cross section.

Building fields:

| | C++ (`spectra::`) | C | Python |
|---|---|---|---|
| bins | `bins_lin`, `bins_log` | `krome_bins_lin`, `krome_bins_log` | `bins_lin`, `bins_log` |
| J <-> Phi (J constant in a bin) | `J_to_phi`, `phi_to_J` | `krome_J_to_phi`, `krome_phi_to_J` | `phi_to_J` |
| blackbody | `blackbody`, `planck_J`, `blackbody_range` | `krome_spectrum_blackbody`, `krome_blackbody_range` | `spectrum_blackbody`, `blackbody_range` |
| Draine 1978 ISRF | `draine` | `krome_spectrum_draine` | `spectrum_draine` |
| J21 power law | `j21` | `krome_spectrum_j21` | `spectrum_j21` |
| Haardt & Madau 2012 | `hm2012` | `krome_spectrum_hm2012` | `spectrum_hm2012` |
| a table J(E) | `table` | `krome_spectrum_table` | `spectrum_table`, `spectrum_file` |
| dust opacity per bin | `dust_kappa_bins`, `dust_kappa_draine` | `krome_dust_kappa_bins`, `krome_dust_kappa_draine` | `dust_kappa`, `dust_kappa_file` |
| G0 and Av of a field | `find_G0_Av` | `krome_find_G0_Av` | `find_G0_Av` |

Every field is the exact integral over the bin. Legacy KROME sampled
the bin midpoint instead.

H2 photodissociation can also come from line data:
`H2 -> H + H [] h2pd() ; pdiss=0.15` sums the X-B lines of
`data/H2pdB.dat` on the photon bins, as legacy
`@next_H2photodissociation` did. See network-format.md section 5.4.

The photo-reaction metadata:

| function | returns |
|---|---|
| `krome_nphoto`, `krome_photo_reaction`, `krome_photo_reactant` | the photo reactions |
| `krome_photo_pdiss` | events per absorbed photon |
| `krome_photo_xsec` | the cross section sigma(E) |
| `krome_opts_photo_sigma` | the bin-averaged cross sections |
| `krome_nlw`, `krome_lw_reaction`, `krome_lw_reactant`, `krome_lw_pdiss`, `krome_opts_lw_weights` | the J_LW rates |

## 9. Binned dust

`DUSTBINS` cooling, `surface.bins` and `dust.evolve_Td` need a binned
grain population: an MRN size distribution per grain type, with Qabs
and Tdust per bin.

Setting it up:
- C++:
  ```cpp
  DustBins d;
  d.nd = 10;
  d.init();
  d.interp_qabs(rad);  // with photon bins
  o.dust = &d;
  ```
- C: `krome_opts_set_dustbins(h, nd, use_Si, use_C, alow, aup, phi)`.
- Python: `opts.set_dustbins(...)`.

### Dust temperature: what KROME does in each mode

| dust description | Tdust | input or computed |
|---|---|---|
| effective (MRN) population, no `dust.tables` | the per-cell parameter `Tdust` (default 10 K), fixed during the solve | **input** |
| effective population, `dust.tables` | interpolated in the legacy tables at every evaluation, from (ntot, Tgas) (plus Av with `3D` / `PHOTO3D`) | **computed**; the parameter is unused |
| binned dust (`DustBins`), default | per bin, the radiative-collisional balance, solved at every evaluation before the rates (see below) | **computed**; always an output, never an input |
| binned dust, `dust.evolve_Td` (Tgas evolved) | per bin, solver state following each grain's energy equation (see below) | **initial condition** in, evolved value out |
| binned dust, `dust.evolve_Td`, Tgas fixed (`evolve_T` off) | the caller's values, fixed (thermal evolution as a whole is off) | **input** |

The balance (default, binned dust) solves, for each bin,

    (intBB(Td) - intBB(Tcmb) - intJ) besc = pre (Tgas - Td):

- the left side is the grain's net radiation: emission minus the CMB
  and the photon field it absorbs;
- the right side is the gas-grain collisional exchange;
- the escape factor `besc` depends on Td through the dust optical
  depth, so the balance is solved to self-consistency (fixed-point
  iteration from Td = Tgas). Legacy took one step of this iteration
  per call, from the previous call's Td.

The balance is a function of the state (Tgas, n, field) only. The rates
(sticking, evaporation, surface reactions, H2 formation on dust) and the
`DUSTBINS` cooling see the same Td, whether Tgas is evolved or not. The
Td you pass in are not used, and the balance at the final state is
written back.

With `dust.evolve_Td`, each bin's Td follows the grain's energy equation:

    C(a, Td) dTd/dt = 4 pi a^2 [pre (Tgas - Td) - (intBB(Td) - intBB(Tcmb) - intJ) besc]

- `C` is the heat capacity of one grain, Draine & Li (2001, eqs. 9-11)
  as published; N (the grain's atoms, in the (N - 2) factor) comes from
  the bulk density and a mean atomic mass.
- Td relaxes to the balance from any start, on the grains' thermal
  timescale.
- The gas cools by the collisional exchange, which conserves energy.
- Legacy's `-usedTdust` instead tracked the balance through dTgas/dt,
  which is valid only when started on it.

Passing and reading the per-bin Td:
- **One cell (C++):** `solve()` reads and writes `dust.Td`.
- **Many cells:** pass a per-cell array `Td[ncell][ndust]`:
  `CellsOpts::Td` in C++, `krome_cells::Td` of `krome_solve_cells` /
  `krome_solve_cells_x` / `krome_equilibrium_cells` in C, `td=` in
  Fortran, `Td=` in Python. Each row receives its cell's final Td; a
  failed cell's row is untouched. With `dust.evolve_Td` the row is also
  the cell's initial condition.
- **The bins' own Td:**
  - C: `krome_opts_set_dust_Td` / `krome_opts_get_dust_Td`;
  - Fortran: `krome_set_dust_Td` / `krome_get_dust_Td`;
  - Python: `opts.set_dust_Td` / `opts.dust_Td`.

  With `dust.evolve_Td` they are every cell's initial condition when no
  per-cell array is given (default 30 K).

### Dust evolution (`dust.evolve`)

`dust.evolve = ["GROWTH", "SPUTTER", "EVAP"]` (legacy `-dustOptions`) makes the
grain radius of every bin solver state. Each bin's grain count stays
fixed; the processes change the radii, and the atoms go to or come from
the gas-phase key element of the grain type: C for carbonaceous grains,
Si for silicates (the network must contain it). Per bin, with m the key
element's mass, rho0 the bulk density (2.25 / 3.13 g cm^-3) and
da = (m / rho0)^(1/3) one atomic layer:

| process | da/dt | from |
|---|---|---|
| `GROWTH` | S(Tg, Td) v n_key m / (4 rho0), v = sqrt(8 kB Tg / (pi m)) | accretion; sticking of Leitch-Devlin & Williams (1985), as legacy |
| `SPUTTER` | -n_p y(Tg) m / (4 rho0) above 1e5 K; n_p = ntot - n_e | thermal sputtering; y: legacy's fit of the Nozawa et al. (2006) yields times the thermal speed |
| `EVAP` | -nu0 da exp(-E0 / Td), nu0 = 1e12 s^-1; E0 = 7.2 eV (carbon, Lenzuni et al. 1995), 4.66 eV (silicate, Nozawa et al. 2006) | sublimation, the Polanyi-Wigner layer removal of Grassi et al. (2017, sec. 2.3): size independent |
| `KELVIN` (with EVAP, optional) | EVAP times exp(2 sigma Omega / (max(a, da) kB Td)), Omega = m / rho0 the volume per key atom, sigma = 1400 (carbon), 800 (silicon) erg cm^-2 (Nozawa et al. 2003) | the curvature (Kelvin) enhancement of a small grain's vapour pressure: smaller grains evaporate faster (x1.3 at 5 nm, x3 at 1 nm, 1500 K); a floored at one layer |

- All rates are multiplied by a^2 / (a^2 + da^2): a grain smaller than one
  layer neither grows nor evaporates, so an evaporated bin stays empty.
- The solver integrates the atoms locked in each bin, so gas plus grains
  is conserved to roundoff. `conserve = 2` is rejected (-106): the
  exchange changes the gas element totals by design.
- Td for these rates: the bins' balance, the `dust.evolve_Td` state, or,
  with `dust.evolve_Td` at fixed Tgas, the caller's Td (legacy's dust test
  set Td = Tgas this way).
- The grain counts are the population's grains per dust mass at its
  reference sizes, times rho_gas x `d2g` of the solve start: `d2g` sets
  the grain number per gas mass, and the dust mass follows the radii.
- The cross sections use the evolved radii (cooling, opacity, the energy
  equation, the surface-chemistry area weights); Qabs stays that of the
  reference sizes, as legacy.
- Dust evolution runs on the host per-cell path and on the device.

Legacy fixes: growth used the proton's thermal speed (3.5 to 5x too
fast); sputtering mixed the initial and current sizes with a fixed
density of 2.3 and counted the electrons; evaporation scaled as 1/a^2
(one atom per grain, far too fast for small grains) and a destroyed bin
was removed inside the rhs with a unit error; `-dustSeed` had no effect.

Passing the radii:
- **One cell (C++):** `solve()` reads and writes `dust.a` (the reference
  sizes after `init()`; 0 = an evaporated bin).
- **Many cells:** `CellsOpts::adust` in C++, `krome_cells::adust` of
  `krome_solve_cells` / `krome_solve_cells_x` / `krome_equilibrium_cells`
  in C, `adust=` in Fortran and Python: `adust[ncell][ndust]`
  [cm], read and updated like `Td`. Hydro codes advect them with the gas.
- **The bins' own radii, counts and types:**
  - C: `krome_opts_set_dust_a`, `krome_opts_get_dust_a`,
    `krome_opts_get_dust_counts` (h, rho_gas, d2g, nd),
    `krome_opts_get_dust_types` (1 carbon, 0 silicate);
  - Fortran: `krome_set_dust_a`, `krome_get_dust_a`,
    `krome_get_dust_counts`, `krome_get_dust_types`;
  - Python: `opts.set_dust_a`, `opts.dust_a`, `opts.dust_counts`,
    `opts.dust_types`.

## 10. Return flags and the error report

Every solver call returns 0 on success.
- **A flag from -1 to -27** means the integration failed after the
  restarts (the integrator's flags; the values are CVODE's). The cell is
  then left untouched: its n, Tgas, Phi, Td and x are unchanged.
- **The KROME flags below** reject the inputs before integrating, except
  `FLAG_NONPOSITIVE_T`. That one fails a solve that would end at a
  non-positive Tgas; a trial state with Tgas <= 0 is also refused during
  the integration, so the integrator shortens its step.
- **Descriptions:** `flag_message` / `krome_flag_message` /
  `krome_flag_message_f` / `net.flag_message` describe any flag.

| C++ | C | Fortran | value | meaning |
|---|---|---|---|---|
| `FLAG_OK` | `KROME_FLAG_OK` | `krome_flag_ok` | 0 | success |
| `FLAG_TOO_MUCH_WORK` | `KROME_FLAG_TOO_MUCH_WORK` | `krome_flag_too_much_work` | -1 | `max_steps` exhausted before t_end |
| `FLAG_TOO_MUCH_ACC` | `KROME_FLAG_TOO_MUCH_ACC` | `krome_flag_too_much_acc` | -2 | tolerances too small for machine precision |
| `FLAG_ERR_FAILURE` | `KROME_FLAG_ERR_FAILURE` | `krome_flag_err_failure` | -3 | repeated error-test failures |
| `FLAG_CONV_FAILURE` | `KROME_FLAG_CONV_FAILURE` | `krome_flag_conv_failure` | -4 | repeated Newton convergence failures |
| `FLAG_LSETUP_FAIL` | `KROME_FLAG_LSETUP_FAIL` | `krome_flag_lsetup_fail` | -6 | Newton matrix setup or factorization failed |
| `FLAG_LSOLVE_FAIL` | `KROME_FLAG_LSOLVE_FAIL` | `krome_flag_lsolve_fail` | -7 | linear solve failed |
| `FLAG_RHSFUNC_FAIL` | `KROME_FLAG_RHSFUNC_FAIL` | `krome_flag_rhsfunc_fail` | -8 | right-hand side failed unrecoverably |
| `FLAG_FIRST_RHSFUNC_ERR` | `KROME_FLAG_FIRST_RHSFUNC_ERR` | `krome_flag_first_rhsfunc_err` | -9 | right-hand side failed at the initial state |
| `FLAG_REPTD_RHSFUNC_ERR` | `KROME_FLAG_REPTD_RHSFUNC_ERR` | `krome_flag_reptd_rhsfunc_err` | -10 | repeated recoverable right-hand side errors |
| `FLAG_UNREC_RHSFUNC_ERR` | `KROME_FLAG_UNREC_RHSFUNC_ERR` | `krome_flag_unrec_rhsfunc_err` | -11 | right-hand side failed unrecoverably in a restart |
| `FLAG_NLS_FAIL` | `KROME_FLAG_NLS_FAIL` | `krome_flag_nls_fail` | -16 | nonlinear solver failure |
| `FLAG_ILL_INPUT` | `KROME_FLAG_ILL_INPUT` | `krome_flag_ill_input` | -22 | illegal input (e.g. error weights not positive) |
| `FLAG_TOO_CLOSE` | `KROME_FLAG_TOO_CLOSE` | `krome_flag_too_close` | -27 | t_end too close to the start |
| `FLAG_NO_DUST` | `KROME_FLAG_NO_DUST` | `krome_flag_no_dust` | -102 | DUSTBINS cooling with Tgas evolved and no dust bins |
| `FLAG_NO_FIELD` | `KROME_FLAG_NO_FIELD` | `krome_flag_no_field` | -103 | photo reactions without bins and field |
| `FLAG_UNSUPPORTED` | `KROME_FLAG_UNSUPPORTED` | `krome_flag_unsupported` | -104 | the requested path cannot do this configuration, or is not built |
| `FLAG_PARAM_UNSET` | `KROME_FLAG_PARAM_UNSET` | `krome_flag_param_unset` | -105 | a REQUIRED parameter is unset (NaN) |
| `FLAG_CONSERVE_GAS_TOTAL` | `KROME_FLAG_CONSERVE_GAS_TOTAL` | `krome_flag_conserve_gas_total` | -106 | conservation with `dust.gas_total` |
| `FLAG_SURFACE_BINS` | `KROME_FLAG_SURFACE_BINS` | `krome_flag_surface_bins` | -107 | `surface.bins = NB` without NB dust bins |
| `FLAG_NO_LAPACK` | `KROME_FLAG_NO_LAPACK` | `krome_flag_no_lapack` | -108 | `linsol` LAPACK in a build without it |
| `FLAG_NO_KOKKOS` | `KROME_FLAG_NO_KOKKOS` | `krome_flag_no_kokkos` | -109 | device path before Kokkos is initialized |
| `FLAG_NONPOSITIVE_T` | `KROME_FLAG_NONPOSITIVE_T` | `krome_flag_nonpositive_t` | -110 | the integration ended at a non-positive or non-finite Tgas; the cell is untouched |

`solve_equilibrium` returns +1 when t_max comes before convergence. The
state at t_max is returned.

### The error report

When an integration fails, KROME writes `KROME_ERROR_REPORT` in the
working directory, like legacy KROME. The file is replaced at each
failure. It contains:
- **The failure:**
  - the flag, its name and its description;
  - in a `solve_cells` / `krome_solve_cells` call, which cell failed;
  - the time reached.
- **The cell properties as passed to KROME:**
  - t_end;
  - the input Tgas (fixed or evolved);
  - rho, ntot, H nuclei, mu and net charge of the input;
  - every parameter;
  - the photon field bin by bin (edges, Phi) with its Habing flux;
  - the dust bins (sizes, input Td);
  - every solver option.
- **The state at the failure:** n, dn/dt and the input n of every
  species.
- **A replay of the failing call:** C/C++ code with 0-based indices
  that sets the input n, Tgas, t_end, all parameters, the options, and
  the photon and dust bins, then calls `solve`. Paste it into a driver
  that includes `krome.hpp`, and it reproduces the failure. Projects
  generated with `project.bindings = ["fortran"]` also get the same call in Fortran, with
  1-based indices, through `krome_solve`.
- **Rates and fluxes:** sorted, with percentages.
- **The integrator statistics:**
  - step sizes, order, counts;
  - failures;
  - the species with the largest local error.

Options:
- **Turn it off:** `error_report = false`.
- **Change the file:** `error_file`. In an MPI run, give each rank its
  own file.

## 11. Migrating from legacy KROME

- **Network files:** all legacy networks are in `networks/` as `.knet`.
  A custom `react_*` file is converted by the generator on the fly (give
  it as `project.network`), or explicitly with
  `python3 kromex/convert_legacy.py react_x -o x.knet --check`.
- **Option files** are not read: a legacy project's options are written
  once as a TOML configuration (section 2). The tables below give, for
  every legacy option, the key that replaces it, or why there is none:
  - **with a key:** the first table;
  - **no generation option:** the option is now a runtime setting, or
    no longer applies;
  - **not ported:** what to use instead.
- **What changed in the calling code:**

  | legacy | KROME-X |
  |---|---|
  | `krome(x, Tgas, dt)` | `solve` (C++), `krome_solve_cells` (C, many cells at once), `krome_solve` (Fortran) |
  | mass fractions (`-useX`) | the `solve_x` family (the solver still integrates number densities) |
  | `krome_idx_*` constants | 1-based in Fortran, 0-based in C/C++ |
  | species order | still alphabetical |
  | `krome_set_user_*` (`user_tff`, `user_crate`, ...) | per-cell parameters (`par` rows, `Params` fields) |

#### With a key

| legacy | KROME-X |
|---|---|
| `-ATOL` | `solver.atol` (the project's default `atol_n`) |
| `-checkReverse` | `rates.check_reverse = true` |
| `-columnDensityMethod` | `radiation.num2col = "jeans"` / `"jeans40"` |
| `-cooling` | `thermal.cooling`, the same names; the per-metal CI, CII, OI, OII, SiI, SiII, FeI, FeII become `Z` in `thermal.cooling` plus `thermal.metals` (e.g. `["C", "C+", "O"]`); H2GP98 was removed in 2.0.0 (use `H2`) |
| `-coolFile` | `thermal.coolant_files`, which ADD to the embedded data/coolZ.dat (legacy replaced it); legacy's data/coolO2.dat is the corrected data/cooling/coolO2.dat |
| `-coolLevels` | `thermal.levels` |
| `-dustOptions` | GROWTH, SPUTTER, EVAP: `dust.evolve`; dT: `dust.evolve_Td = true`; T is the default balance; H2 is the network reaction krate_H2dust |
| `-coolingQuench` | removed in 2.0.0 (an ad hoc switch-off of the cooling below Tcrit; use `thermal.floor`) |
| `-customATOL` | `solver.atol_file` (the same file format; KROME-X warns about species not in the network, which legacy skipped silently) |
| `-customRTOL` | `solver.rtol_file` |
| `-gamma` | `thermal.gamma`: `DEFAULT` = `"CONST"` (the default); EXACT, POPOVAS the same; FULL, VIB, ROT, REDUCED were removed in 2.0.0 (use EXACT); expressions are not ported |
| `-H2opacity` | `thermal.h2_opacity` |
| `-heating` | `thermal.heating`, the same names (XRAY is not one: X-rays are photo reactions) |
| `-interfacePy` | `project.bindings = ["python"]` |
| `-n`, `-network` | `project.network` |
| `-iRHS` | `rates.rhs = "loop"` (the loop rhs; KROME-X keeps the analytic Jacobian, also as loops) |
| `-nuclearMult` | `rates.reactant_multiplicity = true` |
| `-project` | `project.outdir` |
| `-RTOL` | `solver.rtol` (the project's default `rtol`) |
| `-shieldHabingDust` | `radiation.habing_dust_shielding = true` |
| `-shielding` | `radiation.shielding` |
| `-test` | the shipped tests: `python3 -m kromex --test NAME` (section 6.6) |
| `-useCoolFloor` | `thermal.floor = true` |
| `-useIndividualFloor` | `thermal.floor = [TERMS]` |
| `-usePhotoInduced` | `thermal.photo_pumping = true` |
| `-useSemenov` | in `[surface]`: `binding_energies = "SEMENOV"`, `nu = "LEGACY"`, `competition = false`, `tunnelling = false` (and `cr_desorption = "PARAM"` for CR desorption; the single grain size is the network's `surface_grains`, written by the converter) |

#### No generation option (a runtime setting, or not applicable)

| legacy | KROME-X |
|---|---|
| `-dustSeed` | no effect (in legacy either: it was never substituted); grains are not nucleated |
| `-checkConserv` | conservation is checked when the network is parsed; SolverOpts::conserve corrects numerical drift |
| `-checkThermochem` | `rates.check_reverse = true` |
| `-clean` | not applicable (`project.outdir`) |
| `-compact` | the C++ core is header-only; `project.single = true` writes the network as one header |
| `-compiler` | deprecated in legacy; the Makefile's CXX / FC |
| `-computeElectrons` | `SolverOpts::conserve = 3` (C key `conserve`): after each solve the electrons neutralize the ions, also while a `freeze` directive holds e- fixed (legacy computed them in the rhs). Without a freeze every reaction conserves charge (checked at parse time) and conserve = 1 corrects the drift |
| `-conserve` | conservation is a runtime option in KROME-X, not a generation flag: SolverOpts::conserve = 2 (elements and charge; option key conserve) |
| `-conserveE` | charge conservation is a runtime option in KROME-X, not a generation flag: SolverOpts::conserve = 1 (option key conserve) |
| `-conserveLin` | conservation is a runtime option in KROME-X, not a generation flag: SolverOpts::conserve = 2 (the weighted least-squares projection, no LAPACK) |
| `-dry` | not applicable |
| `-enzo` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-fexArgument` | not applicable (KROME-X integrates with CVODE's BDF, ported) |
| `-flash` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-forceMF21` | the Jacobian is chosen at runtime: SolverOpts::jac |
| `-forceMF222` | the Jacobian is chosen at runtime: SolverOpts::jac = JAC_FD |
| `-forceRWORK` | KROME-X uses CVODE's BDF, ported (no work arrays) |
| `-gizmo` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-ierr` | every KROME-X call returns a flag |
| `-interfaceC` | the C interface (krome.h) is always built |
| `-lh` | no Fortran headers are generated |
| `-maxord` | the BDF order is CVODE's (up to 5); not an option |
| `-mergeTlimits` | piecewise rates are one reaction group in KROME-X (reaction_names lists each piece) |
| `-nochargeCheck` | KROME-X always validates the network when parsing it (mass, charge, nucleons): errors cannot be skipped; `python3 -m kromex --check` also reports reactions written twice |
| `-noCheck` | KROME-X always validates the network when parsing it (mass, charge, nucleons): errors cannot be skipped; `python3 -m kromex --check` also reports reactions written twice |
| `-noExample` | the examples are written once and never overwritten |
| `-nomassCheck` | KROME-X always validates the network when parsing it (mass, charge, nucleons): errors cannot be skipped; `python3 -m kromex --check` also reports reactions written twice |
| `-noRecCheck` | KROME-X does not stop on missing recombinations; `python3 -m kromex --check` reports species never destroyed |
| `-noSinkCheck` | KROME-X does not stop on sinks; `python3 -m kromex --check` reports species never destroyed or never formed |
| `-noVerbatimFile` | reaction names come from the network |
| `-pedantic` | deprecated in legacy |
| `-photoBins` | the photon bins are runtime setup in KROME-X (Radiation + PhotoModule::init in C++, krome_opts_set_photo in C); photo reactions come from the network |
| `-photoDustVarAv` | Av is the per-cell parameter Av |
| `-photoDustVarG0` | G0 is the per-cell parameter G0 |
| `-quote` | not ported |
| `-quotelist` | not ported |
| `-ramses` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-ramsesOffset` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-ramsesTH` | KROME-X is a library called by the hydro code through the C/Fortran interface: no code-specific patches are generated |
| `-report` | the error report is a runtime option in KROME-X, not a generation flag: on by default (SolverOpts::error_report, error_file) |
| `-sh` | no Fortran headers are generated |
| `-skipDevTest` | the shipped tests (`--test NAME`) have no development status |
| `-skipDup` | KROME-X always validates the network when parsing it (mass, charge, nucleons): errors cannot be skipped; `python3 -m kromex --check` also reports reactions written twice |
| `-skipJacobian` | the analytic Jacobian is always generated; SolverOpts::jac = JAC_FD selects difference quotients at runtime |
| `-skipODEthermo` | Tgas is evolved or not at runtime: SolverOpts::evolve_T |
| `-source` | not applicable |
| `-unsafe` | generated files are overwritten, user files (Makefile, examples) never |
| `-useBroadening` | the Doppler broadening is set with the photon bins at runtime: Radiation::broad_T / broad_vturb (C krome_opts_set_photo_broadening) |
| `-useChemoToggle` | not ported: the chemistry is always evolved (evolve_T is the one runtime physics switch) |
| `-useDvodeF90` | KROME-X uses CVODE's BDF, ported |
| `-useEquilibrium` | equilibrium is a runtime call in KROME-X: solve_equilibrium / krome_equilibrium_cells / net.equilibrium |
| `-useFileIdx` | not ported: reactions are numbered in network order |
| `-useIERR` | every KROME-X call returns a flag |
| `-useN` | number densities are the default |
| `-usePhIoniz` | obsolete in legacy; photo reactions are in the network |
| `-usePhotoOpacity` | obsolete in legacy |
| `-usePlainIsotopes` | isotopes are written in the network |
| `-useTabs` | tabulated rates were measured and removed (a 7-10 % gain); the implementation is in the git history |
| `-useThermoToggle` | thermal evolution is a runtime option in KROME-X, not a generation flag: SolverOpts::evolve_T (option key evolve_T) |
| `-useX` | mass fractions is a runtime option in KROME-X, not a generation flag: solve_x / solve_cells_x (C++), krome_solve_cells_x (C), krome_solve_x (Fortran), net.solve_x (Python) |
| `-v` | `python3 -m kromex --version` |
| `-ver` | `python3 -m kromex --version` |
| `-verbatimFilename` | reaction names come from the network |
| `-version` | `python3 -m kromex --version` |

#### Not ported

| legacy | KROME-X |
|---|---|
| `-customODE` | custom ODEs are not ported |
| `-dust` | binned dust is runtime setup in KROME-X (DustBins: nd, useSi, useC, then init(); C krome_opts_set_dustbins) with the `DUSTBINS` cooling |
| `-dustTabs` | `dust.tables` (`"HM2012"`, `"MC"`, `"3D"` or `"PHOTO3D"`; the tabulated H2-on-dust rate is the network reaction krate_H2dust_table(); the dust cooling table is the `DUSTTAB` cooling) |
| `-listAutomatics` | not ported |
| `-listSWRI` | not ported |
| `-noTlimits` | `rates.outside = "clamp"` (clamp both ends) instead |
| `-reducer` | not ported: the legacy reducer ranked reactions by flux without checking the reduced network; network reduction is planned as a separate tool that writes and validates a reduced .knet |
| `-reverse` | reverse reactions are written in the network (<->), or added by the converter: kromex/convert_legacy.py --reverse |
| `-stars` | nuclear networks: the `NUCLEAR` heating, `rates.screening = true`, `rates.reactant_multiplicity = true` |
| `-Tlimit` | `rates.outside` (`"asym"`, `"clamp"` or `"zero"`) for the rate temperature ranges |
| `-useAutoNetwork` | not ported |
| `-useCustomCoe` | not ported: rates come from the network |
| `-useDustH2const` | the Jura rate is the network reaction krate_H2dust_Jura() (parameter clump) |
| `-useODEConstant` | not ported (rate_scale = expr in the network scales every rate) |
| `-xsecKernelFunction` | not ported |

#### The legacy tools

| legacy | KROME-X |
|---|---|
| `tools/docmake` | `python3 -m kromex doc NETWORK` (§2, "Documenting a network"): the species and reaction pages, the rate plots, the joints of piecewise rates, the missing reactions and branches, the graphs per element and the LaTeX table, from the `.knet` (KIDA and UMIST files: `tools/kida2krome.py`, `tools/umist2krome.py`, §2). Not ported: the Langevin rates (the polarizability table had no documented source) and the subnetwork filter |

#### The legacy `krome_user` functions

Every public legacy user function (legacy `src/krome_user.f90`, the
public routines of `src/krome.f90`, and the per-network functions legacy
generated) and its KROME-X form. Conventions:
- C++ names are in namespace `krome`; `p` is the `Params` struct, `o`
  the `SolverOpts`; `h` a C/Fortran `krome_opts` handle; `net` a Python
  `Network`, `opts` an `Options`.
- "par X" is the per-cell parameter X: C++ `p.X`, C
  `par[KROME_PAR_X]`, Fortran `par(krome_par_X, :)`, Python
  `par={"X": ...}`.
- "plain" means an ordinary array operation in the caller's code (the
  state, the field `Phi` and the bin edges are the caller's arrays).

**Solve and state**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome(x, Tgas, dt[, xdust, Tdust])` | integrate one cell over dt | `solve(n, Tgas, p, dt, o)`; many cells `solve_cells`; per-bin dust `CellsOpts::Td` / `adust` | `krome_solve_cells` (extras in `krome_cells`) | `krome_solve(n, Tgas, dt, h, par=, phi=, flags=, td=, adust=)` | `net.solve(n, Tgas, dt, par=, Phi=, opts=)`: `Result.n`, `.Tgas`, `.Td` |
| `krome(x, rhogas, Tgas, dt)` (`-useX`) | the same with mass fractions | `solve_x`, `solve_cells_x` | `krome_solve_cells_x` | `krome_solve_x` | `net.solve_x` |
| `krome_equilibrium` | integrate to steady state | `solve_equilibrium` (`EquilOpts`) | `krome_equilibrium_cells` (keys `eq_*`) | `krome_equilibrium` | `net.equilibrium` |
| `krome_init` | load tables, set defaults | nothing to initialize (photon bins `PhotoModule::init`, dust `DustBins::init`) | `krome_opts_new` | `krome_opts_new` | `krome.load()` |
| `krome_thermo_on` / `krome_thermo_off` | thermal evolution on / off | `o.evolve_T` | key `evolve_T` | `krome_set(h, "evolve_T", 1d0)` | `opts["evolve_T"]` |
| `krome_thermo(x, Tgas, dt)` | explicit Euler step of Tgas | `get_thermal` (dT/dt), then Tgas += dT/dt dt | `krome_get_thermal` | `krome_get_thermal` | `net.thermal(...).dTdt` |
| `krome_store` / `krome_restore` | state to / from `krome_dump.dat` (under development in legacy) | plain (the state is the caller's arrays) | plain | plain | plain |
| `krome_dump` (internal) | failure dump | the error report: `o.error_report`, `o.error_file` | keys `error_report`, `krome_opts_set_error_file` | `krome_set_error_file` | `opts["error_report"]`, `opts["error_file"]` |
| `krome_set_mpi_rank` | the rank in failure output | one `o.error_file` per rank | `krome_opts_set_error_file` | `krome_set_error_file` | `opts["error_file"]` |

**Per-cell physical settings** (per-cell parameters in KROME-X)

| legacy | what | KROME-X (every interface) |
|---|---|---|
| `krome_set_Z` | metallicity Z/Z_sun | par `total_Z` |
| `krome_set_dust_to_gas` | dust-to-gas ratio D/D_sun | par `d2g`, a mass ratio: d2g = 0.01 x the legacy value |
| `krome_set_clump` | Jura H2-on-dust clumping factor | par `clump` |
| `krome_set_Tcmb` / `krome_get_Tcmb` | CMB temperature | par `Tcmb` |
| `krome_set_zredshift` / `krome_get_zredshift` | redshift | par `zredshift` |
| `krome_set_orthoParaRatio` / `krome_get_orthoParaRatio` | H2 ortho:para | par `opratio` |
| `krome_set_Tfloor` / `krome_get_Tfloor` | cooling-floor temperature | par `Tfloor` |
| `krome_set_metallicity` / `krome_get_metallicity` | `phys_metallicity` | none: legacy never read the value it set |
| `krome_set_user_<v>` / `krome_get_user_<v>` | the network's `@common: user_*` variables | the network's `param <v>` (`user_` dropped; `user_Tdust` is the built-in par `Tdust`) |
| `krome_set_dust_table_3D_variable` | the third axis (Av) of the 3D dust tables | par `Av` (`dust.tables = "3D"`) |
| `krome_set_J21xray` | the tabulated X-ray J21 | not ported: X-rays go through the photon bins |

**Species and reaction metadata**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_get_names` | species names | `species_names[]` | `krome_species_name(i)` | `krome_species_name_f(i)` | `net.species` |
| `krome_get_index` | index of a species | `IDX_<species>` | `krome_species_index` | `krome_species_idx` | `net.idx` |
| `krome_get_names_header` | the names, space-separated | `names_header()` | `krome_names_header` | `krome_names_header` | `net.names_header()` |
| `krome_get_rnames` | reaction names | `reaction_names[]` | `krome_reaction_name(r)` | `krome_reaction_name_f(r)` | `net.reactions` |
| `krome_get_mass` / `krome_get_imass` | masses [g] / their inverses | `species_mass_g[]` | `krome_species_mass(i)` | `krome_species_mass_f(i)` | `net.masses` |
| `krome_get_charges` | charges | `species_charge[]` | `krome_species_charge` | `krome_species_charge_f` | `net.charges` |
| `krome_get_zatoms` | atomic number per species | the composition: `species_elements` | `krome_species_natoms`, `krome_element_name` | `krome_species_natoms_f` | `net.composition` |
| `krome_get_cooling_names_header` / `krome_get_heating_names_header` | the thermal term names | `cooling_name(i)`, `heating_name(i)` | `krome_ncool`, `krome_cooling_name`, `krome_nheat`, `krome_heating_name` | `krome_cooling_name_f`, `krome_heating_name_f` | `net.cooling_terms`, `net.heating_terms` |
| `krome_get_info` | the species table | `print_info` | `krome_print_info` | `krome_print_info` | `net.print_info` |

**Rates and fluxes**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_get_coef`, `krome_get_coe`, `krome_get_coeT` | rate coefficients | `get_rates` | `krome_get_rates` | `krome_get_rates` | `net.rates` |
| `krome_get_flux` | reaction fluxes | `get_fluxes` | `krome_get_fluxes` | `krome_get_fluxes` | `net.fluxes` |
| `krome_print_best_flux` | the largest fluxes | `print_best_flux` | `krome_print_best_flux` | `krome_print_best_flux` | `net.print_best_flux`, `net.best_fluxes` |
| `krome_explore_flux` | every flux, absolute and relative | `explore_flux` | `krome_explore_flux` | `krome_explore_flux` | `net.explore_flux` |
| `krome_dump_flux` | index and flux to a unit | `explore_flux` (a superset) | `krome_explore_flux` | `krome_explore_flux` | `net.explore_flux` |
| `krome_dump_rates` | k on a Tgas grid to a file | `dump_rates` | `krome_dump_rates` | `krome_dump_rates` | `net.dump_rates` |
| `krome_get_qeff` | nuclear Q per reaction | `get_qeff` | `krome_get_qeff` | `krome_get_qeff` | `net.qeff` |
| `krome_stars_coe` | screened nuclear rates | `get_rates` (`rates.screening`) | `krome_get_rates` | `krome_get_rates` | `net.rates` |
| `krome_stars_energy` | flux x Q per reaction | `get_qeff` x `get_fluxes` [MeV] | the same | the same | `net.qeff` x `net.fluxes` |

**Thermal**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_get_heating` / `krome_get_heating_array` | total / per-term heating | `get_thermal` (`heat`, `heats[]`) | `krome_get_thermal` | `krome_get_thermal` | `net.thermal(...)`: `.heat`, `.heating` |
| `krome_get_cooling` / `krome_get_cooling_array` | total / per-term cooling | `get_thermal` (`cool`, `cools[]`) | `krome_get_thermal` | `krome_get_thermal` | `.cool`, `.cooling` |
| `krome_cooling<metal>` (generated) | the cooling of one coolZ metal | `get_coolant_cooling` | `krome_get_coolant_cooling` | `krome_get_coolant_cooling` | `net.coolant_cooling` |
| `krome_popcool_dump` | the coolants' level populations | `get_coolant_populations`, `print_coolant_populations` | `krome_get_coolant_populations`, `krome_print_coolant_populations` | the same | `net.coolant_populations`, `net.print_coolant_populations` |
| `krome_dump_cooling` | the cooling terms to a unit | `dump_cooling` | `krome_dump_cooling` | `krome_dump_cooling` | `net.dump_cooling` |
| `krome_print_phys_variables` | print the physical variables | `print_params` | `krome_print_params` | `krome_print_params` | `net.print_params` |
| `krome_plot_cooling` | a cooling plot | not ported | not ported | not ported | not ported |

**Gas properties, units and conservation**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_get_gamma` / `krome_get_gamma_x` | adiabatic index (from n / from x) | `get_gamma` (with x: `x2n` first) | `krome_get_gamma` | `krome_get_gamma` | `net.gamma` |
| `krome_get_mu` / `krome_get_mu_x` | mean molecular weight | `get_mu` (with x: `x2n` first) | `krome_get_mu` | `krome_get_mu` | `net.mu` |
| `krome_get_rho` | mass density | `get_rho` | `krome_get_rho` | `krome_get_rho` | `net.rho` |
| `krome_get_Hnuclei` | H nuclei density | `get_Hnuclei` | `krome_get_Hnuclei` | `krome_get_Hnuclei` | `net.hnuclei` |
| `krome_get_electrons` | electrons for neutrality | `get_electrons_neutral` | `krome_get_electrons` | `krome_get_electrons` | `net.electrons` |
| `krome_get_free_fall_time` | free-fall time from n | `get_free_fall_time` | `krome_get_free_fall_time` | `krome_get_free_fall_time` | `net.free_fall_time` |
| `krome_get_free_fall_time_rho` | free-fall time from rho | `get_free_fall_time_rho` | `krome_get_free_fall_time_rho` | `krome_get_free_fall_time_rho` | `net.free_fall_time_rho` |
| `krome_num2col` | the column-density estimator | `num2col` | `krome_num2col` | `krome_num2col` | `net.num2col` |
| `krome_n2x` / `krome_x2n` | number densities <-> mass fractions | `n2x`, `x2n` | `krome_n2x`, `krome_x2n` | `krome_n2x`, `krome_x2n` | `net.n2x`, `net.x2n` |
| `krome_consistent_x` | renormalize x, electrons for neutrality | after a solve: `o.renormalize_x`, `o.conserve = 3` | keys `renormalize_x`, `conserve` | `krome_set` | `opts[...]` |
| `krome_conserve`, `krome_conserveLin_x` | force element and charge conservation | after a solve: `o.conserve = 2` (against the input totals) | key `conserve` | `krome_set(h, "conserve", 2d0)` | `opts["conserve"] = 2` |
| `krome_conserveLinGetRef_x` | the reference element totals | `get_element_totals` | `krome_get_element_totals` | `krome_get_element_totals` | `net.element_totals` |
| `krome_scale_Z` | metals at Asplund+09 x 10^Z | `scale_Z` | `krome_scale_Z` | `krome_scale_Z` | `net.scale_Z` |

**Photon bins and radiation** (the field is the caller's photon flux per
bin `Phi`; legacy's `photoBinJ` converts with `J_to_phi` / `phi_to_J`)

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_calc_photobins`, `krome_set_photobinE_lr`, `krome_set_photobinE_limits` | the bins from their edges | `Radiation{Elo, Ehi}` + `PhotoModule::init` | `krome_opts_set_photo(h, nb, Elo, Ehi)` | `krome_set_photo` | `opts.set_photo(Elo, Ehi)` |
| `krome_set_photobinE_lin` / `_log` | linear / log bins | `spectra::bins_lin`, `Radiation::log_bins` | `krome_bins_lin`, `krome_bins_log` | `krome_bins_lin`, `krome_bins_log` | `net.bins_lin`, `net.bins_log` |
| `krome_set_lineBroadeningVturb` (and the Tgas argument of `set_photobinE_*`) | Doppler broadening | `Radiation::broad_T`, `broad_vturb` | `krome_opts_set_photo_broadening` | `krome_set_photo_broadening` | `opts.set_photo_broadening` |
| `krome_set_photoBinJ` / `krome_get_photoBinJ` | set / get the field | the caller's `Phi` (`o.Phi`); `spectra::J_to_phi`, `phi_to_J` | the `Phi` argument; `krome_J_to_phi`, `krome_phi_to_J` | `phi=`; the same | `Phi=`; `net.J_to_phi`, `net.phi_to_J` |
| `krome_get_photoBinE_left` / `_right` / `_mid` / `_delta` / `_idelta` | bin edges, centres, widths | plain (the caller's Elo, Ehi) | plain | plain | plain |
| `krome_get_photoIntensity` | J of the bin containing E | plain (`phi_to_J` + a bin lookup) | plain | plain | plain |
| `krome_photoBin_scale` / `_scale_array` / `_store` / `_restore` | scale, store, restore the field | plain | plain | plain | plain |
| `krome_opacity_scale_size` | Phi x exp(-tau) over a size | plain: Phi x exp(-`get_absorption` x size) | `krome_get_absorption` | `krome_get_absorption` | `net.absorption` |
| `krome_get_photoBin_rates` | the photo rates | `get_rates` at `photo_reaction[r]` | `krome_get_rates` + `krome_photo_reaction` | `krome_get_rates` + `krome_photo_reaction_f` | `net.rates(...)[..., net.photo_reactions]` |
| `krome_get_photoBin_heats` | photoheating per photo reaction | `get_photo_heating` (per cm^3; divide by the absorber for legacy's per-absorber value) | `krome_get_photo_heating` | `krome_get_photo_heating` | `net.photo_heating` |
| `krome_get_xsec` | cross sections per bin | `PhotoModule::sigma`; pointwise `photo_xsec[r](E)` | `krome_opts_photo_sigma`; `krome_photo_xsec` | `krome_photo_sigma`; `krome_photo_xsec_f` | `opts.photo_sigma()`; `net.photo_xsec` |
| `krome_load_photoBin_file_2col` | J(E) from a file onto the bins | `spectra::table` | `krome_spectrum_table` | `krome_spectrum_table` | `net.spectrum_file`, `net.spectrum_table` |
| `krome_load_photoBin_file` | edges and J from a file | the caller reads the file; `Radiation` + `J_to_phi` | `krome_opts_set_photo` + `krome_J_to_phi` | the same | `opts.set_photo` + `net.J_to_phi` |
| `krome_set_photoBin_HMlog` / `_HMCustom` | the Haardt & Madau field | `spectra::hm2012` (+ bins) | `krome_spectrum_hm2012` | `krome_spectrum_hm2012` | `net.spectrum_hm2012` |
| `krome_set_photoBin_BBlin` / `_BBlog` / `_BBlog_auto` | a blackbody | `spectra::blackbody` (`blackbody_range` for the auto range) | `krome_spectrum_blackbody`, `krome_blackbody_range` | the same | `net.spectrum_blackbody`, `net.blackbody_range` |
| `krome_set_photoBin_draineLin` / `_draineLog` / `_draineCustom` | the Draine ISRF | `spectra::draine` | `krome_spectrum_draine` | `krome_spectrum_draine` | `net.spectrum_draine` |
| `krome_set_photoBin_J21lin` / `_J21log` | a J21 power law | `spectra::j21` | `krome_spectrum_j21` | `krome_spectrum_j21` | `net.spectrum_j21` |
| `krome_get_ratioFluxDraine` | the field over Draine's | `get_ratio_flux_draine` | `krome_get_ratio_flux_draine` | `krome_get_ratio_flux_draine` | `net.ratio_flux_draine` |
| `krome_dump_Jflux` | E and J to a unit | plain | plain | plain | plain |
| `krome_get_opacity` | optical depth per bin | `get_opacity` | `krome_get_opacity` | `krome_get_opacity` | `net.opacity` |
| `krome_get_opacity_size` / `_size_d2g` | optical depth over a size | `get_absorption` x size (dust: `kdust`) | `krome_get_absorption` | `krome_get_absorption` | `net.absorption` |
| `krome_load_opacity_table`, `krome_load_average_kabs` | dust opacity onto the bins | `spectra::dust_kappa_bins`, `dust_kappa_draine` | `krome_dust_kappa_bins`, `krome_dust_kappa_draine` | the same | `net.dust_kappa_file`, `net.dust_kappa` |
| `krome_find_G0_Av` | G0 and Av of a field | `spectra::find_G0_Av` | `krome_find_G0_Av` | `krome_find_G0_Av` | `net.find_G0_Av` |

**Dust**

| legacy | what | C++ | C | Fortran | Python |
|---|---|---|---|---|---|
| `krome_get_table_Tdust`, `krome_get_Tdust` (tables) | Tdust of the rates | `get_Tdust` | `krome_get_Tdust` | `krome_get_Tdust` | `net.Tdust` |
| `krome_init_dust_distribution` | MRN bins | `DustBins{nd, useSi, useC, alow, aup, phi}.init()`, `o.dust`; d2g is par `d2g` | `krome_opts_set_dustbins` | `krome_set_dustbins` | `opts.set_dustbins` |
| `krome_get_dust_distribution` | grains per bin [cm^-3] | `DustBins::xw` x rho x d2g | `krome_opts_get_dust_counts` | `krome_get_dust_counts` | `opts.dust_counts` |
| `krome_set_dust_distribution` | set the grains per bin | none (the MRN weights `DustBins::xw` are fixed at `init`) | none | none | none |
| `krome_scale_dust_distribution` | scale the grain numbers | par `d2g` (the counts scale with it) | par `d2g` | par `d2g` | par `d2g` |
| `krome_get_dust_size` / `krome_set_dust_size` | bin radii | `DustBins::a` | `krome_opts_get_dust_a` / `_set_dust_a` | `krome_get_dust_a` / `krome_set_dust_a` | `opts.dust_a` / `set_dust_a` |
| `krome_set_Tdust`, `krome_set_Tdust_array`, `krome_get_Tdust` (bins) | per-bin Td | `DustBins::Td` (per cell `CellsOpts::Td`); without bins par `Tdust` | `krome_opts_set_dust_Td` / `_get_dust_Td` | `krome_set_dust_Td` / `krome_get_dust_Td` | `opts.set_dust_Td` / `dust_Td` |
| `krome_get_averaged_Tdust` | mean Td | `get_Tdust` with `o.dust` (area-weighted; legacy weighted by number) | `krome_get_Tdust` | `krome_get_Tdust` | `net.Tdust` |
| `krome_set_surface`, `_surface_norm`, `_surface_array`, `krome_get_surface` | the per-bin surface species | plain, on `IDX_X_b<j>_dust` (weights from `krome_opts_get_dust_counts`) | plain | plain | plain |

**Not ported**

| legacy | what | why |
|---|---|---|
| `krome_convert_xmoc`, `krome_return_xmoc`, `krome_set_photobinE_moc` | MOCASSIN coupling | code-specific, as the RAMSES / ENZO patches |
| `krome_get_jacobian` | a Jacobian column | internal to the solver |
| `krome_print_best_flux_frac`, `krome_print_best_flux_spec` | filtered flux prints | `best_fluxes` / `get_fluxes` and the metadata filter them |
| `krome_plot_cooling` | a cooling plot | plotting is left to the caller |
| `krome_reducer` | the network reducer | `-reducer` stops (above) |
