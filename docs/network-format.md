# KROME v2 network format — specification draft (v1)

Status: **draft v2** (Step 1 of the rewrite; decisions from the first review
round are folded in: `clamp` out-of-range default, `Tgas_eV`, electron `e-`,
case-sensitive species, `.knet` extension).
Scope: gas-phase chemistry, photochemistry, X-ray and cosmic-ray chemistry,
freezeout/evaporation surface chemistry, nuclear networks, custom
cooling/heating terms, reverse (thermochemical) reactions (§3.2), and
grain-surface chemistry (§5.3).

Guiding principles:

1. **Human-readable reactions** — prizmo-style arrow syntax, no positional
   columns, no `@format` state.
2. **Rates are a defined expression language, not Fortran** — every rate is
   parsed into an AST with a whitelisted symbol table. This is what makes
   multi-backend code generation (Fortran / C++ / CUDA), automatic Jacobians,
   and CSE possible.
3. **The network file describes chemistry only.** Solver choices (tabulation,
   tolerances), thermal-module selection, and host-code coupling live in the
   build configuration, not in the network. Legacy directives that encoded
   solver behavior (`@noTab*`, `@storeonce`) have no equivalent here on purpose.
4. **Everything is validated at parse time** — unknown symbols, unbalanced
   mass/charge, overlapping temperature ranges are errors before any code is
   generated.

---

## 1. File structure

- Plain text, UTF-8. Recommended extension: `.knet`.
- `#` starts a comment (full-line or trailing). Blank lines are ignored.
- A line ending in `\` continues on the next line (useful for long fits).
- The first non-comment line must be the format header:

  ```
  krome_format 1
  ```

- Every other non-comment line is either a **directive** (`param`, `var`,
  `table`, `extern`, `species`, `ref`, `cooling`, `heating`, and the
  network-level modifiers `rate_scale` and `freeze` of §6.2) or a
  **reaction**.
- Directives are order-independent: a `var` may be used before the line that
  defines it. Cyclic definitions are an error.

## 2. Species

```
species := isotope? unit+ label* phase? charge? | "e-"
unit    := element digits?
isotope := digits            # mass number: 2H, 3He (nuclear networks)
label   := "_" alnum+        # state/site label: O_1D, CH2_1, H_c_dust
phase   := "_ice" | "_dust"  # grain-surface species
charge  := "+"+ | "-"+
```

- Element symbols are **case-sensitive** standard symbols (`H`, `He`, `Si`,
  ...). `D` is deuterium. Examples: `H2`, `HD`, `He++`, `H-`, `CH3OH`,
  `CO_ice`, `3He`, `O_1D` (excited state), `H_c_dust` (chemisorbed site).
  Labeled species are distinct species sharing the core's composition.
- A name that does not parse as an element composition (`FK1` didactic
  species, `M` generic collider) is accepted as an **opaque species**: zero
  reference mass, conservation checks skipped for its reactions, and a
  prominent parser notice — so placeholders work but typos stay visible.
- `e-` is the free electron.
- Photons are **not written**: `H+ + e- -> H` is the radiative recombination
  (legacy `g` product is simply dropped). Cosmic-ray particles and X-ray
  photons are likewise implicit (see §5.2, §6).
- In a reaction, a reactant in curly braces — `H + H + {H} -> 3He` — is
  **non-kinetic**: it is consumed by the stoichiometry (mass/charge balance
  includes it) but its density does not enter the rate product. Used by
  nuclear networks where a lumped reaction consumes more particles than the
  rate-limiting step involves (legacy "curly species").

Species not appearing in any reaction (legacy `@ghost`) are declared
explicitly:

```
species H2, He        # tracked/advected but chemically inert in this network
```

The parser derives mass and charge from the name and checks conservation on
every reaction (overridable per reaction only via an explicit annotation, TBD).

**H2 and its nuclear-spin species.** A network may carry H2 as one
species (`H2`) or as `H2_ortho` and `H2_para` (or all three). The processes
that need H2 read the total, `nH2_(n)` (H2 cooling, CIE, the H2 opacity,
the H2 self-shielding `krome_fshield`, `shield_dust`, the molecular
coolants' H2 collider, the CHEM / PHOTOAV heating), and those that need
ortho and para separately read `nH2o_` / `nH2p_` / `h2_ofrac_`: the
network's own `H2_ortho` / `H2_para`, else `H2` split by the parameter
`opratio` (the metal-line H2or / H2pa colliders, the rotational heat
capacity of `thermal.gamma = "EXACT"` / `POPOVAS`). `opratio` exists only in networks
that need that split. The Solomon and UV dissociation of any H2 form get
their heating with their own reactant.

### 2.1 Species order in the generated code

The generated arrays (`n[NSPEC]`, `IDX_<species>`, `species_names`, the C
`KROME_IDX_*` and Fortran `krome_idx_*` constants) order the species by
name, sorted by code point (so `D, D+, H, H+, H-, H2, ..., He, e-`):
independent of how the network file is written. Code should use the
generated index constants or look names up (`krome_species_index`), not
hard-coded positions. (Legacy KROME used first appearance in the file
followed by a stable sort by charge.)

## 3. Reaction lines

```
reaction   := stoich "->" stoich  trange  rate  annotation*
stoich     := ("{" species "}" | species) ("+" ...)*
trange     := "[" [number] "," [number] "]" | "[]"
annotation := ";" key "=" value
```

Examples:

```
H + e- -> H+ + e- + e-     []              exp(-32.71 + 13.54*lnTe - ...)
H+ + e- -> H               [0, 5.5e3]      3.92e-13 * invTe**0.6353
H+ + e- -> H               [5.5e3, inf]    exp(-28.61 - 0.724*lnTe - ...)
H + H + H -> H2 + H        [0, 3e2]        1.3e-32 * (Tgas/300)**(-0.38)
H + 2H -> 3He              [1e0, 1.1e8]    (1.81e3/T923)*exp(-3.721/T913)  ; qeff=5.493
```

Recognized annotations (all optional): `qeff=<expr>` — energy released per
reaction event [MeV] (consumed by `thermal.heating = ["NUCLEAR"]`, and returned by the
`get_qeff` query);
`pdiss=<number>` in (0, 1] — reaction events per absorbed photon, for
`photo()` and `J_LW` rates only (§5.4; required for `J_LW` rates).
Reactions preceded by a `ref` directive carry that citation as metadata:

```
ref \cite{Janev1987Springer}
H + e- -> H+ + e- + e-     []    exp(...)
```

`ref` attaches to the next reaction line; it has no runtime effect and is
collected into the build report (bibliography of the generated network).

- Reactions have **no index numbers**; the preprocessor assigns stable indices
  and exports a name map for host codes. An optional label `id NNN:` prefix is
  reserved if explicit indexing proves necessary.
- Fields are separated by whitespace (any amount; alignment is cosmetic).
- Up to 3-body reactions in v1 (limit is a parser setting, not a format limit).
- Rate units are the usual cgs conventions: `s^-1`, `cm^3 s^-1`, `cm^6 s^-1`
  by reaction order; densities are `cm^-3`.

### 3.1 Temperature ranges and piecewise rates

- `[Tmin, Tmax]` bounds in K. `0` / `inf` denote unbounded sides; `[]` means
  no restriction.
- A rate valid on `[Tmin, Tmax]` applies for `Tmin <= Tgas < Tmax`
  (lower-inclusive, upper-exclusive). Legacy `.gt.` vs `.ge.` limit operators
  collapse onto this half-open convention.
- The range field always refers to `Tgas`. A gate on any other quantity
  (e.g. a radiation-temperature threshold) is written in the rate itself with
  `where(...)` (§4.1).
- Lines with **identical reactants and products** are grouped as follows:
  - Lines whose ranges **tile** an interval (no overlap, no interior gap) form
    the **pieces of one piecewise rate**. This replaces the legacy
    operator-prefixed limits (`.LE.`, `>`, ...) and `-mergeTlimits`.
  - Several such tiled chains (or unlimited `[]` lines) for the same
    stoichiometry are **separate additive channels** — physically distinct
    processes (gas-phase vs on-dust recombination of `H+ + e- -> H`; direct
    LW photodissociation vs the Solomon process for `H2 -> H + H`; radiative
    vs three-body association of `C+ + O -> CO+`, each its own piecewise
    chain). The parser emits a notice listing the channels so accidental
    duplicates are caught by eye; a piece that neither continues nor overlaps
    any chain (an interior gap) is an error.
  - A **degenerate range** `[x, x]` is legal with a notice: under `clamp` it
    means "always evaluate at T = x" (a rate measured at a single
    temperature). Legacy semantics made such a reaction fire only at exactly
    that temperature — almost certainly never the intent.
- **Outside the declared union** of a piecewise group the default is
  **asymmetric**: below the lower bound the rate is **zero** (single-sided
  lower limits in practice belong to barrier reactions whose true rate
  collapses — clamping would inflate it by orders of magnitude), above the
  upper bound the rate is **clamped** to its value at the bound (fits blow
  up beyond their validity; saturation is the safe continuation). Uniform
  `clamp` and uniform `zero` (legacy KROME behavior) are build-time options,
  globally or per reaction. This default was chosen after measuring an 86%
  HD error from low-T clamping in the primordial regression at 100 K.

### 3.2 Reverse reactions: `<->`

```
N2 + O <-> NO + N        []            3e-10*exp(-38300.0/Tgas)
H + H+ <-> H2+           [10, 6.7e3]   1.85e-23*Tgas**1.8
H + H+ <-> H2+           [6.7e3, 1e8]  5.81e-16*(Tgas/5.62e4)**(-0.6657*log10(Tgas/5.62e4))
```

`<->` declares the forward reaction with the given rate and adds its
reverse, right after it, from detailed balance:

  k_rev(T) = k_fwd(T) / K_c(T),  K_c = exp(-dG°/RT) (P° / kB T)^dn,

with P° = 1 bar, and dG° and dn = n_products - n_reactants those of the
forward reaction. The Gibbs energies come from the NASA 7-coefficient
polynomials of `data/thermo_burcat.toml` (Burcat's database, Goos,
Burcat & Ruscic, for the species of the shipped networks, by name), or
else the NIST-JANAF Shomate fits of `data/thermoNIST.dat`. Outside a species' tabulated range its G is
taken at the nearest limit. A species without data is a generation error
(e.g. `e-`, which has no entry).

- The reverse copies the forward's temperature range and `where()`
  conditions (each piece of a piecewise rate is reversed piecewise). K_c
  is evaluated at the actual gas temperature, also where the forward rate
  is held at its range limit.
- Allowed for expression and `auto` rates between kinetic species. It is
  rejected for photo/surface/`J_LW` rates, curly reactants, `qeff`
  reactions, and when the reverse is also written in the network.
- The converter's `--reverse` reproduces legacy `-reverse` (every
  eligible reaction written as `<->`, without duplicating reverses
  already present). `codegen rates.check_reverse` writes `reverse_check.txt`
  (legacy `-checkReverse`): each reaction's reverse in the network or its
  absence, and dH(300 K). A `*` marks a missing reverse with |dH| <= 3000 K,
  >= 2 reactants and products, and no electrons.

## 4. Rate expressions

A small, backend-neutral expression language. Grammar:

```
expr    := term (("+"|"-") term)*
term    := factor (("*"|"/") factor)*
factor  := ["+"|"-"] base ["**" factor]
base    := number | symbol | call | "(" expr ")"
call    := name "(" expr ("," expr)* ")"
number  := decimal with optional e-notation exponent   # 1.3e-32, 6.0e-10, 300
```

- **Only `e` exponents.** Fortran `d` notation and `.LE.`-style tokens are not
  part of the language (the legacy converter translates them).
- No assignments, no conditionals, no loops. Anything that needs branching is
  either a piecewise reaction (§3.1) or an extension function (§4.4).

### 4.1 Built-in functions

`exp, log, log10, sqrt, abs, min, max, tanh` (natural log is `log`), and the
functional conditional `where(cond, a, b)` — `a` where `cond` holds, else
`b`. `cond` is a comparison (`<, <=, >, >=`) between expressions, valid only
as the first argument of `where`; **only the selected branch is evaluated**
(a ternary/`merge` in the backends), so it is a safe division-by-zero guard.
It exists for the rare cases the `[Tmin, Tmax]` field cannot express — a gate
on a non-`Tgas` quantity, or a guard on a computed value:

```
He+ -> He++ + e-    []    where(Trad > 1e3, 5.0e1*Trad**1.63*exp(-5.9e5/Trad), 0.0)
H + 14N + {H} -> 12C + 4He    [1e0, 2.5e9]    where(R11_1 + R12_1 > 0.0, \
                                    R10*R11_1/(R11_1 + R12_1), 0.0)
```

Piecewise-in-`Tgas` rates should use the range field, not `where`. The legacy
`if(cond): expr` rate guard converts to `where(cond, expr, 0.0)`.

### 4.2 Predefined state symbols

| symbol    | meaning                                   |
|-----------|-------------------------------------------|
| `Tgas`    | gas temperature [K]                       |
| `Tgas_eV` | gas temperature [eV] (= `Tgas * kB/eV`)   |
| `Tdust`   | dust temperature [K]                      |
| `ntot`    | total number density [cm^-3]              |
| `Hnuclei` | total H-nuclei number density [cm^-3]     |
| `J_LW`    | Lyman-Werner band (11.2-13.6 eV) mean intensity of the photon bins [eV s^-1 cm^-2 Hz^-1 sr^-1] (§5.4) |
| `n(X)`    | number density of species `X` [cm^-3] (an isotope with its mass number first: `n(60Fe)`) |

`ntot` and `Hnuclei` are stoichiometry-derived sums generated per network
(with exact symbolic Jacobian contributions); in the gas+total abundance
representation (`dust.gas_total`, see the manual) their per-species weights
adapt automatically to avoid double counting.

Deliberately minimal: legacy shortcuts (`Te`, `lnTe`, `invT`, `T32`,
`sqrTgas`, ...) are **not predefined** — they are ordinary `var` definitions,
which the converter emits automatically when a legacy network uses them:

```
var Te   = Tgas_eV
var lnTe = log(Tgas_eV)
var invT = 1.0/Tgas
var T32  = Tgas/300.0
```

`n(X)` inside a rate or `var` makes the rate density-dependent; the symbolic
Jacobian picks up the extra terms automatically (this removes the whole
`@noTab*` class of hacks).

### 4.3 Directives: `param`, `var`, `table`

```
param crate = 1.3e-17      # cosmic-ray ionization rate [s^-1]
param Av    = 0.0          # visual extinction
```

`param` declares an **environment parameter**: constant during one solver
call, settable per cell/batch through the runtime API. Replaces
`@common:user_*` (the `user_` prefix is gone; the generated API namespaces
them). The value in the file is the default.

Some parameters are predefined: they can be used in expressions without
declaration, and the runtime physics modules read them: `G0` (local
radiation field [Habing], default 0), `d2g` (dust-to-gas mass ratio,
default 0.01), `Av` (visual extinction [mag], default 0), `Tcmb`,
`zredshift`, and those of specific thermal processes. `Tdust` is used as a
state symbol (§4.2). A generated network has **only the parameters it
uses**: its own `param`s, the predefined ones its expressions or its
configured physics read, plus `gamma_ad` and `Tcmb` (the temperature the
rates see is capped to [Tcmb, 1e9 K]; the thermal terms see it capped to
[2.73 K, 1e9 K]). A parameter of a configured process
for which no default exists (e.g. `tff` for compression heating, `hubble0`
and `omega0` for expansion cooling or the ST93 recombination) is REQUIRED:
the solver rejects a call that leaves it unset.

```
var kl11 = 10**(-27.029 + 3.801*log10(Tgas) - 29487.0/Tgas)
```

`var` declares a named sub-expression usable in any rate. File-scope,
order-independent, cycle-free. Vars are pure expressions over state symbols,
params, tables, and other vars; codegen decides evaluation order and CSE.

```
table shield_H2 = "tables/shield_WG11.dat"
```

`table` binds a name to a tabulated function loaded at runtime init;
call it like a function: `shield_H2(NH2, Tgas)`. Dimensionality (1D/2D) comes
from the file header. Replaces `@tabvar`. The grid may be non-uniform
(rectilinear); the function is linear (1D) or bilinear (2D) in its
arguments and clamped at the edges. With `log` after the file name

```
table fshield_CO = "data/shield_CO_unigrid.dat" log
```

the interpolation is in log10 of the axes: the file keeps its values of
the axis quantities (all > 0), and the function is called with their
log10, e.g. `fshield_CO(log10(NH2), log10(NCO))`; the tabulated values
themselves are interpolated as they are.

### 4.4 Extension functions: `extern`

```
extern H2_solomonLW          # provided by the 'photo' module
extern H_recombination_on_dust
```

`extern` declares a function implemented by a physics module of the runtime
library (dust recombination, H2 self-shielding, sticking/evaporation rates,
...). The preprocessor errors out if the build configuration does not enable a
module providing it. State arguments are implicit (the legacy `(n, Tgas)`
argument lists disappear); genuine arguments (e.g. a flux, a scale factor) are
kept, and a **bare species name is a valid argument** — it resolves to the
species index (legacy `idx_CO` becomes just `CO`):

```
CO -> CO_ice    []    krate_stickSi(CO)
CO_ice -> CO    []    krate_evaporation(CO)
```

Inside `n(...)` and extern arguments, a token that exactly matches a species
name is a species literal, never an expression — so `n(e-)` and `n(H+)` are
unambiguous (legacy Fortran-safe index names decode accordingly:
`idx_Hj` → `H+`, `idx_Dk` → `D-`, `idx_E` → `e-`). An extern argument may
also be a **species list** in brackets, for thermochemical helpers that take
whole stoichiometries (legacy nucleation networks):

```
Si2O2 -> SiO + SiO   []   cluster_growth_rate(SiO, 1) * revKc_with_GFE([Si2O2], [SiO, SiO])
```

Kinetic-nucleation helpers (Gail & Sedlmayr 2013, ch. 13; generated):

- `cluster_growth_rate(M, N[, stick])`: monomer M added to an N-cluster,
  v_th pi (r1 + r1 N^(1/3))^2 stick, with the legacy monomer radii (TiO2,
  Al2O3, MgO, SiO; another monomer is a generation error).
  `general_cluster_growth_rate(M, N1, N2[, stick])` is the same for two
  clusters.
- `revKc_with_GFE([R...], [P...])`: K_c of R -> P, exp(-dG/RT) (P0/kB
  T)^dn, from the Gibbs energies of `data/thermochemistry/<NAME>.gfe`
  (linear in T, clamped to 100-3040 K), at the actual gas temperature.
- `cluster_destruction_rate(M, N[, stick])` in N-mer -> (N-1)-mer + M:
  detailed balance, `cluster_growth_rate(M, N-1[, stick])` times the
  reaction's own `revKc_with_GFE`.
- `troe_falloff(k0, kinf, Fc, M)`: simplified Troe falloff, 0 when k0 or
  kinf underflow.

### 4.5 Physical constants

A small set of predefined named constants is available in expressions:
`Rgas_kJ` (gas constant, kJ mol^-1 K^-1), `boltzmann_eV`, `boltzmann_erg`,
`N_avogadro`, `pi`, `planck_eV` (eV s), `eV_to_erg`. Their values are
CODATA 2018 (kromex/constants.py). Backends inline their values; anything
else belongs in a `param` or `var`.

## 5. Photochemistry

A photo-reaction is an ordinary reaction whose rate is the `photo(...)`
constructor. The runtime computes `k = ∫ σ(E) J(E)/E dE` over the configured
photon bins; a shielding factor and the photons absorbed per event are
declared with the reaction (§5.4).

```
H -> H+ + e-          []    photo(xsec=auto)
H- -> H + e-          []    photo(xsec="swri", Emin=0.755)
H2+ -> H+ + H         []    photo(xsec="leiden", Emin=2.65)
H2+ -> H+ + H+ + e-   []    photo(xsec = 10**(-16.926 - 4.528e-2*E_eV \
                                  + 2.238e-4*E_eV**2 + 4.245e-7*E_eV**3), \
                                  Emin=30.0, Emax=90.0)
H2 -> H + H           []    photo(xsec = H2_sigmaLW(E_eV), Emin=14.159, Emax=17.7)
```

- `xsec=` sources: `auto` (cross-section database lookup by
  reactants/products), a quoted database/file reference (`"swri"`,
  `"leiden"`, `"xsecs/myfile.dat"`, `"verner96"`), or an inline expression of
  the reserved symbol `E_eV` (photon energy, eV) returning cm^2.
- CO: use `xsec="leiden"` (`data/database/leiden_xsecs/CO__C_O.dat`,
  Heays et al. 2017 from the line data of Visser et al. 2009, resampled
  to 0.1 nm keeping the integrated cross section): 2.4e-10 s^-1 in the
  unshielded Draine field, the standard rate. The SWRI table misses most
  of CO's predissociating lines (7.5e-11 s^-1). The lines saturate at
  columns where a binned cross section does not: give the self-shielding
  as `shield=`, e.g. the Visser et al. (2009) functions of
  `data/shield_CO_unigrid.dat` (a network `table`).
- `auto` photoionizations (`data/database/photoionization.dat`): the
  total cross section of Verner's `phfit2.f` (version 2, 1996): the
  Verner et al. (1996) outer-shell fit up to the next inner-shell edge,
  above it the sum over the shells of the Verner & Yakovlev (1995) fits,
  each from its own threshold (`sigma_vy95`), up to 5e4 eV and
  extrapolated above.
- Auger channels: `X -> X^(q+k) + k e-` with `photo(xsec=auto)` (k > 1,
  next to `X -> X^(q+1) + e-`, which the network must have) takes the
  multiple ionizations of the inner shells: a vacancy in shell s ends
  with j electrons out in all (the photoelectron included) with the
  probability P_s(j) of Kaastra & Mewe (1993,
  `data/database/auger_yields.dat`), and each j goes to the network's
  largest k not above j. So the channels share the total cross section
  (the opacity is the same whichever are present), a network with
  `X -> X+ + e-` only keeps it all (as before), and the outer-shell fit
  below the edge is single ionization. Shells without a yield (outer
  shells, inner ones without Auger decay) have P(1) = 1; an ion with no
  inner shell (H- and He-like) cannot have an Auger channel (an error).
  A k no shell reaches has a rate of 0 (a warning).
- The energy an `auto` absorption deposits (photoheating, secondary
  ionization) is, below the inner-shell edge, the photoelectron's
  E - Emin; above it, averaged over the shells by their share of the
  cross section, each shell's photoelectron E - Eth_s plus the mean
  energy of its Auger electrons (Kaastra & Mewe; the fluorescence photons
  escape), whichever channels the network has.
- `Emin`/`Emax` (eV) bound the integration window for **any** source,
  `auto` included; `Emin` defaults to the database threshold when available.
- The `[]` range field on photo lines refers to `Tgas` as usual (photo rates
  are usually T-independent; energy limits live inside `photo()`).
- This replaces `@photo_start/stop` blocks, `@xsecFile=`, and the bare
  `krome_kph_auto` rate tokens.
- `photo()` declares only the cross section. Whether the radiation field is
  held constant during a chemistry step (legacy behavior) or **co-evolved
  inside the solver** (photon bins as log-flux state variables) is a
  build/solver option, not a network property: `radiation.evolve` (see the
  manual).

### 5.1 Piecewise cross sections

Repeated photo lines for the same reaction define **σ(E) piecewise**: the
energy windows must not overlap (a parse error — an overlap would
double-count the integral), each segment may use a different source (database
segments and inline expressions can be mixed), and the rate is the sum of the
segment integrals. This is the energy-space analogue of piecewise T-rates.

```
H2 -> H2+ + e-    []    photo(xsec = 9.56e-17*(E_eV/15.4) - 9.4e-17, Emin=15.4, Emax=16.5)
H2 -> H2+ + e-    []    photo(xsec = 2.16e-17*(E_eV/15.4) - 1.48e-17, Emin=16.5, Emax=17.7)
H2 -> H2+ + e-    []    photo(xsec = 1.51e-17*(E_eV/15.4)**(-2.71),   Emin=17.7, Emax=30)
```

### 5.2 X-ray chemistry

There is **no separate X-ray construct**: X-ray chemistry is ordinary
photochemistry with photon bins extending to X-ray energies, using the same
`photo()` cross sections (which is why the `auto` resolver must carry the
full energy coverage, inner shells included, rather than the truncated
legacy ranges). Legacy `@XRAY_start/stop` blocks and their tabulated rates
(`ratexH.dat`-style primary+secondary tables) are retired; the converter
maps X-ray `auto` entries onto `photo(xsec=auto)` and drops duplicates of
reactions already present as photo. Secondary ionization by fast
photoelectrons becomes a photo-module feature (energy-deposition model),
not a rate-table property.

### 5.3 Grain-surface chemistry

Surface species carry an `_ice` / `_dust` label (any case: `X_DUST` too);
`H_c_dust` is chemisorbed H. Reactions between them use the extension
functions (§4.4) or the `surface(...)` constructor, resolved by the surface
model on the dust population:

```
H2O -> H2O_ice               []  surface(auto)                 # adsorption
H2O_ice -> H2O               []  surface(auto)                 # desorption
H_dust + O_dust -> OH_dust   []  surface(auto)                 # two-body
O_dust + CO_dust -> CO2_dust []  surface(Ea=1000.0, branch=0.5)
H_dust + CO_dust -> HCO_dust []  surface(prob=5.03e-5)
H_dust -> H_c_dust           []  surface(auto)                 # chemisorption
ebind HO2 = 4000.0, 4000.0   # bare / ice E_bind [K] (illustrative value)
CO -> CO_ice                 []  krate_stickSi(CO)             # legacy functions
CO -> CO_ice                 []  surface(stick=1.0)            # given sticking
surface_grains single(n = xdust, a2 = gsize2)  # one grain size (below)
```

- **`surface_grains single(n = <expr>, a2 = <expr>)`** (or `a = <expr>`)
  replaces the default population with one grain size: n_d grains cm^-3 of
  radius^2 a2 [cm^2], expressions of parameters and constants (a per-cell
  property). This is legacy's `-useSemenov` population (`user_xdust`,
  `user_gsize2`): sites S = 4 pi n_d a^2 / (3 A)^2 and adsorption n_d pi a^2
  v_th S. It ignores a runtime DustBins population (Tdust is the cell's),
  and `surface.bins` is rejected with it. `surface_grains mrn` states the
  default. The converter writes it for legacy Semenov-framework networks.
- **Population:** the effective MRN grains holding rho d2g, or,
  when a DustBins population is given at runtime, its bins (rates
  area-weighted, each bin at its own Tdust). With codegen
  `surface.bins = NB` every surface species and reaction is instead cloned
  per bin (`X_b<j>_dust`, legacy @surface), and the runtime binning must
  have NB bins. Sites: S = 4 pi sum n_d a^2 / (3 A)^2. Water-ice coverage:
  f_ice = min(n(H2O ice) / S, 1). Every rate mixes the bare and ice surfaces
  by f_ice.
- **Adsorption** (gas X -> X_ice): v_th sum n_d pi a^2 S(Tgas, Tdust), with
  the HM79 sticking, or with the reaction's own sticking coefficient:
  `surface(stick=<expr>)` (legacy `dust_adsorption_rate(..., stick, ...)`).
- **Desorption** (X_ice -> X, or an isomer of X): nu exp(-E_bind/Tdust).
  With codegen `surface.cr_desorption = "PARAM"` it adds cosmic-ray spot heating
  (Hasegawa & Herbst 1993), f70 nu exp(-E_bind / 70 K), where f70 =
  3.16e-19 PARAM / 1.3e-17.
- **Two-body** (Langmuir-Hinshelwood):
  - k = branch kappa (h_A + h_B) / S, with h = nu exp(-E_diff/Tdust).
    E_diff/E_bind is 2/3 with the KROME binding energies and 0.77 with
    Semenov's (`surface.binding_energies`).
  - A factor 1/2 applies for identical reactants.
  - Reaction-diffusion competition (Chang+07, Garrod & Pauly 11): kappa =
    P_r / (P_r + h_A + h_B + d_A + d_B), with P_r = max(nu) prob, and prob =
    max(exp(-Ea/Tdust), tunnelling through Ea over 1 A) for `Ea=`, the
    given `prob=`, or the database probability for `surface(auto)`.
    `surface.competition = false` gives the legacy kappa = prob.
    `surface.tunnelling = false` drops the tunnelling term, prob =
    exp(-Ea/Tdust) (legacy Semenov framework).
- **nu:** Hasegawa et al. 1992, sqrt(2 n_s E_bind / (pi^2 m)); or
  `surface.nu = "LEGACY"`, 1e12 s^-1.
- **Reactive desorption:**
  - As written (default): a gas-product line (`H_dust + O_dust -> OH`)
    takes the database fraction f, and its surface-product partner 1 - f.
  - `surface.reactive_desorption = "FIXED"`: the gas branch is generated from the
    database fractions.
  - `surface.reactive_desorption = "MINISSALE"` (Minissale+16): CD = exp(-E_bind N /
    (eps dH_R)), divided by 10 on ice. Only the surface-product reactions
    are written then.
- **Chemisorption** (H_dust <-> H_c_dust, H + H_c, H_c + H_c): the legacy
  PC / CP / CC rate tables vs Tdust.
- **`ebind X = E_bare, E_ice`** completes or overrides the binding-energy
  set. The data files are matched by exact name first, then
  case-insensitively (Semenov's file writes HE, SI, MGH2). A process that needs a missing binding energy is a generation error
  (legacy used 0 or 1e99 silently).

### 5.4 Shielding, photons per event, and `J_LW` rates

The photochemistry removes the photons it uses: the `radiation.evolve` sink and
the `get_absorption` / `get_opacity` queries count, for every channel,
(rate) / `pdiss` absorbed photons per reaction event.

```
var NH2 = 1.87e21*(n(H2)*1e-3)**(2.0/3.0)
# Solomon process (line pumping, 11.2-13.6 eV, self-shielded), either as a
# J_LW rate ...
H2 -> H + H   []  H2_solomonLW(J_LW) * fselfH2_WG11(NH2)   ; pdiss=0.15
# ... or (not both) as a photo() rate with an effective line cross section
# (3.748e-18 cm^2 reproduces H2_solomonLW for a flat J)
H2 -> H + H   []  photo(xsec=3.748e-18, Emin=11.2, Emax=13.6, \
                        shield=fselfH2_WG11(NH2))              ; pdiss=0.15
# direct photodissociation (Lyman/Werner continua, 14.16-17.7 eV): a
# continuum absorption, attenuated as exp(-sigma N) by the opacity itself,
# with NO self-shielding factor
H2 -> H + H   []  photo(xsec=H2_sigmaLW(E_eV), Emin=14.159, Emax=17.7)
```

- `photo(..., shield=<expr>)`: a shielding factor s(T, n) multiplying the
  photo rate, k = s sum_b sigma_b Phi_b. It scales the reaction's
  photoheating, secondary ionization and photon removal equally. It is
  meant for line processes whose saturation the flat-in-bin continuum
  cannot represent (self-shielding functions such as fselfH2*); a
  continuum absorber is attenuated by the opacity alone. Dust shielding
  factors (`shield_dust(gamma)`, legacy networks) are the approximation
  for a field that nothing else attenuates: when the photon bins are
  attenuated by the caller (radiative transfer, or `get_opacity` /
  `get_absorption` with dust), the dust must not also be in that
  opacity, or it is counted twice.
- `J_LW`: the band-mean intensity of the photon bins over 11.2-13.6 eV,
  J_LW = sum_b J_b ln(overlap_b) / ln(13.6/11.2), with each bin constant in
  J. A flat J gives J_LW = J, and without photon bins J_LW = 0. A rate using
  it acts on the bins as a photo reaction: it must be **linear** in J_LW
  (J_LW times a J_LW-free factor, or through a linear extern such as
  `H2_solomonLW`), with one kinetic reactant (the absorber), no temperature
  range, and a `pdiss`. Its shielding is part of the rate expression. Its
  photons are taken only from the Lyman-Werner bins, in proportion to the
  photons each has in the band. `J_LW` may appear only in reaction rates
  (not in `var`s or cooling/heating terms).
- Bins: put photon-bin edges at 11.2 and 13.6 eV. In a bin straddling a
  band edge, `J_LW` counts the in-band part with J constant across the
  bin, while photo() cross sections take the photons as uniform in E, so
  the two Solomon forms agree exactly only with aligned edges (the
  difference is confined to the two edge bins).
- Heating: the Solomon process (a `J_LW` or `H2_solomonLW` rate, or a
  photo() H2 -> H + H with its window inside 11.2-13.6 eV) deposits no
  photoelectron energy. Its heating is the Solomon term of
  `thermal.heating = ["CHEM"]` ((18.7 f + 0.4) eV per dissociation, f the thermalized
  fraction) or of `PHOTOAV` ((16.85 f + 0.4) eV), not both. The direct
  photodissociation continuum keeps the photo() deposit E - Emin and gets
  no pumping heat.
- `; pdiss=`: reaction events per absorbed photon. It is 1 by default for
  `photo()`, and required for `J_LW` rates. For the Solomon process it is
  the dissociation probability per pump (~0.15): every pump absorbs a
  photon.
- Column form (`get_opacity`): when a shielding factor (or a `J_LW` rate)
  contains exactly one `fselfH2`, `fselfH2_WG11` or `fselfH2_R14` call, the
  optical depth is the curve of growth sigma_b/pdiss int_0^N s dN', with
  the call's column argument swept up to the num2col column of the absorber.
  Otherwise it is s(N) sigma_b N / pdiss.
- A rate that reads a flux parameter or constant instead of `J_LW` (e.g.
  `H2_solomonLW(myflux_eV)`) is not connected to the bins and removes no
  photons.


**`h2pd()`: H2 photodissociation from line data.** `H2 -> H + H [] h2pd()
; pdiss=<p>` is legacy's `@next_H2photodissociation` (`kpd_H2`). The rate
is the sum over the H2 X-B lines of `data/H2pdB.dat`:

    k = sum_lines pre_l J_b(dE_l) z_X(l)(Tgas)

- `pre` holds pi e^2 f eps / (m_e c dE), with f the oscillator strength
  and eps the dissociation probability.
- `J_b` is the bin's intensity, constant across the bin as in legacy:
  J_b = Phi_b h / (4 pi ln(Ehi/Elo)).
- `z` are the normalized Boltzmann populations of the ground vibrational
  levels. They are computed relative to the lowest level, so they do not
  underflow at low T as legacy's did.
- A line outside every bin contributes nothing (legacy stopped).

How it is built:
- The reaction becomes one photo channel per X level. Each channel has
  its lines as bin cross sections and z_X(Tgas) as its shielding factor.
  So the photoheating rules, the `radiation.evolve` sink, the opacity queries
  and both Jacobians treat it like any photo reaction.
- It is the Solomon process: no photoelectron deposit. With CHEM it gets
  the Solomon heating; legacy's H2pd got none.
- `; pdiss=` is required. The data fold each line's eps into its
  prefactor, so the photons absorbed per dissociation need an effective
  value (0.15 is the usual Solomon one). Legacy removed no photons.

## 6. Cosmic-ray reactions

No special block: a CR-induced rate is a plain expression proportional to the
`crate` param.

```
param crate = 1.3e-17
H  -> H+ + e-   []    0.46 * crate
H2 -> H + H     []    0.10 * crate
```

## 6.1 Custom cooling and heating terms

A network may define named volumetric cooling/heating terms
[erg cm^-3 s^-1] over the same expression language (they may use `n(X)`,
`var`s, tables). They are *defined* here because they reference network
species, but *enabled* in the build configuration like any other thermal
process. Replaces `@cooling_start/stop` + `@cooling:` blocks.

```
var a       = 1.35 - 0.3*log10(Tgas*1e-3)
var lambda3 = 1.32e-23*(Tgas*1e-3)**a
cooling H2O_rot = (n(H2) + 1.39*n(H)) * n(H2O) * lambda3
```

## 6.2 Network-level modifiers: `rate_scale` and `freeze`

Two declarative replacements for the legacy raw-Fortran blocks
`@reactionModifier` / `@ODEmodifier`. Being part of the network, they
are visible to the preprocessor: the analytic Jacobian stays exact, the
selfcheck verifies them, and the same code serves CPU and GPU.

```
rate_scale = 1 / N_avogadro          # NACRE rates are per mole
freeze e-, H, H2, H+ where zredshift > 2.243e3
freeze He                            # always
```

- `rate_scale = <expr>` (at most one): every rate coefficient of the
  network is multiplied by the expression (legacy `k(:) = k(:)/C`). It
  may use anything a rate may use; the preprocessor folds it into each
  rate expression (T-pieces included) after `auto` and `surface`
  resolution. Not allowed together with `photo(...)` reactions: a photo
  rate is sigma x flux, so a scale belongs in the cross section.
- `freeze <species, ...> [where <condition>]` (any number): the listed
  species keep their abundance, i.e. their dn/dt is held at zero while
  the condition holds (legacy `if (cond) dn(idx_X) = 0`). The reactions
  themselves still run, so the other species see them. The solver zeroes
  the species' rows of the full right-hand side (chemistry, secondary
  ionization, radiation coupling) and of the analytic Jacobian. The
  condition is a single comparison of **parameters and constants only**,
  so it is fixed for the duration of a solver call (a switch on the
  evolving state would make the right-hand side discontinuous); it is
  evaluated per cell, from that cell's parameters. Without `where` the
  species are always frozen.

Per-reaction switches need no directive: `where(cond, 0.0, rate)` in the
rate itself (legacy `if (cond) k(a:b) = 0`).

## 7. Database rates: `auto`

The bare rate token `auto` asks the preprocessor to resolve the rate from the
bundled rate database (photoionization, collisional ionization, radiative
recombination, KIDA, ...) by matching reactants and products, exactly as
legacy `auto`. The resolved rate is dumped into the build report so it is
never a black box. A reaction the database lacks, or whose entry needs
machinery KROME-X does not have, stops the generation with an error (as
in legacy), never a rate of 0. On a surface reaction (a species with an
`_ice` / `_dust` label) a bare `auto` is `surface(auto)` (§5.3).

The database is `data/database/*.dat`, read as `hydroHe.dat` (the
primordial H/He set), then `recombination_badnell.dat`, then the other
files in sorted order; the first entry that matches wins. Electron
recombinations X^q+ + e- -> X^(q-1)+ are Badnell's total radiative
(Badnell 2006) plus dielectronic (Badnell et al. 2003 and the sequence
papers) fits for the ground level, wherever Badnell has the ion, before
the astrochemistry databases (KIDA's X+ + e- fits have no dielectronic
part); Verner's radiative fits otherwise (`radiative_recombination.dat`;
both files generated from the verbatim sources in
`data/database/sources/`). Collisional ionization is Voronov
(1997), as distributed in Verner's `cfit.f`.

```
He -> He+ + e-    []    photo(xsec=auto)
O+ + e- -> O      []    auto
```

## 8. Validation (parser guarantees)

Errors: unknown symbol/function, undeclared param, mass or charge imbalance,
duplicate reaction with overlapping or gapped T ranges, cyclic `var`
definitions, `extern` not provided by an enabled module, malformed numbers
(including `d` exponents), a `freeze` of an unknown species or with a
state-dependent condition, a repeated `rate_scale` or one combined with
photo rates. Warnings: unused `var`/`param`/`table`, rate pieces
whose union is narrower than a configurable "expected" range.

## 9. What is deliberately *not* in the network file

| legacy                                | new home                                   |
|---------------------------------------|--------------------------------------------|
| `@noTab*`, `@storeonce*`, `-useTabs`  | build config (tabulation is a solver choice)|
| `@format`                             | gone (fixed syntax)                        |
| `@CR_start/stop`                      | plain expressions with `crate` (§6)        |
| `@photo_*`, `@xsecFile`               | `photo(...)` constructor (§5)              |
| `@XRAY_start/stop`, `ratexH` tables   | `photo(...)` with X-ray energy bins (§5.2) |
| `@ice:`, `@surface_*`                 | explicit reactions with `_ice`/`_dust` species (§5.3) |
| `@cooling:`/`@heating:` blocks        | `cooling`/`heating` directives (§6.1)      |
| `@ghost:`                             | `species` directive (§2)                   |
| `@ref:` / `#@ref:` comments           | `ref` directive (§3)                       |
| `@tabvar:`                            | `table` + `var` (§4.3)                     |
| `@common:user_*`                      | `param` (§4.3); `user_Tdust` → builtin `Tdust` |
| `qeff` column                         | `; qeff=` annotation (§3)                  |
| `-cooling=`, `-heating=`, `-dust=`    | build config                               |
| `@reactionModifier`                   | `rate_scale` (§6.2) or `where()` in the rates |
| `@ODEmodifier`                        | `freeze ... where` (§6.2); raw code has no equivalent |
| index column / `-useFileIdx`          | auto-assigned + exported name map          |

## 10. Reserved syntax (not in v1)

- `id NNN:` explicit reaction labels.
- Three-phase (surface / mantle) grain chemistry.
- Further per-reaction annotations (e.g. `outside=zero`, `nocheck=charge`)
  beyond `qeff`, if experience shows we need them.

## 11. Legacy conversion

`kromex/convert_legacy.py` translates old `react_*` files mechanically:
`d`→`e` literals, operator T-limits → intervals, `@format` columns → arrow
syntax, `@var`/`@common`/blocks → `var`/`param`/constructors, shortcut
symbols → emitted `var` preamble, photon products dropped, `E` → `e-`,
species names recased (`HE` → `He`), `idx_X`/state arguments cleaned from
extension-function calls, and **commented-out legacy reaction blocks are
converted too** (kept as comments, but in the new syntax, so they can be
re-enabled by just uncommenting). Modifier blocks are translated when
they match the corpus patterns: `k(:) = k(:)/C` → `rate_scale = 1 / C`,
`if (cond) k(a:b) = 0` → `where(cond, 0.0, rate)` on reactions a..b
(legacy sequential numbering), `if (cond) then; dn(idx_X) = 0; ...` →
`freeze X, ... where cond`, rate references `k(N)` → the referenced rate
inlined; any other block stays a marked comment. With `--check` it numerically compares
every plain analytic rate before/after conversion over the declared T range.
