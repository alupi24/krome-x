// KROME-X: the public C++ interface.
//
//   #include "krome.hpp"      // the only header a C++ caller includes
//
// Everything a caller uses is declared and documented in this file. The
// implementation (header-only, so that the network's kernels inline and
// compile for devices) lives in core/, and its names are internal: they
// may change between versions. The C interface is krome.h (libkrome.a),
// with the optional Fortran (krome.f90) and Python (krome.py) bindings on
// top of it. docs/manual.md is the user manual.
//
// Build: compile your code with -I<project> and link net_<name>_data.o
// (no external library: the integrator, a port of CVODE's BDF, is part of
// the headers); the project Makefile's `example` target shows the line.
// Optional: -DKROME_USE_LAPACK (+ LAPACK, Makefile LAPACK=),
// -DKROME_USE_KOKKOS (+ Kokkos, Makefile KOKKOS=), -fopenmp (cells in
// parallel).
//
// Units: cgs throughout. n [cm^-3], Tgas [K], times [s], rho [g cm^-3],
// photon flux per bin Phi [photons cm^-2 s^-1], energies of photon bins
// [eV], thermal rates [erg cm^-3 s^-1].
//
// ===========================================================================
// The network (generated: net_<name>_core.hpp, see NETWORK.md in the
// project for this network's values)
// ===========================================================================
//
//   NSPEC, NREA, NPAR, NPHOTO, NELEM   numbers of species, reactions,
//                                      parameters, photo reactions, elements
//   IDX_<species>                      index into n[NSPEC] (species sorted
//                                      alphabetically by code point; H+ ->
//                                      IDX_Hj, H- -> IDX_Hk, e- -> IDX_E)
//   NETWORK_NAME                       the network name
//   species_names[NSPEC]               "H", "H+", ...
//   species_mass_g[NSPEC]              mass [g]
//   species_mass_amu[NSPEC]            mass [amu]
//   species_charge[NSPEC]              charge [e]
//   element_names[NELEM], species_elements[NSPEC * NELEM]
//                                      composition (atoms of element e in
//                                      species i at [i * NELEM + e])
//   reaction_names[NREA]               "H + e- -> H+ + e- + e-", ...
//
//   struct Params                      the runtime parameters of one cell,
//                                      one field per parameter (p.tff,
//                                      p.crate, ...). Only the parameters
//                                      this network uses exist. A field
//                                      without a default is NaN: the
//                                      physics chosen at generation needs
//                                      it, and the solver rejects the cell
//                                      (FLAG_PARAM_UNSET) until it is set.
//   param_names[NPAR]                  the field names, in the order of the
//                                      C interface's par[] rows
//   param_required_by[NPAR]            what needs the parameter ("" = it
//                                      has a default)
//   double* param_ptr(Params&, int i)  the field of parameter i
//   int params_missing(const Params&)  index of the first unset required
//                                      parameter, or -1
//
// The physics (coolers, heaters, shielding, the adiabatic index, the
// evolved blocks) is fixed when the project is generated (codegen flags);
// at runtime only the numerics below and SolverOpts::evolve_T change.

#pragma once
#include <cstdio>
#include <string>
#include <vector>

// the network: the one already included (a net_<name>.hpp), else this
// project's (krome_network.hpp, generated)
#ifndef KROME_NETWORK_INCLUDED
#include "krome_network.hpp"
#endif
#include "core/special_rates.hpp"
#include "core/h2_shielding.hpp"
#include "core/dust_rates.hpp"
#include "core/cooling.hpp"
#include "core/photo.hpp"
#include "core/dust_bins.hpp"

namespace krome {

// physical constants: krome::constants, defined by the generated network
// header (kromex/constants.py: CODATA 2018, IAU), e.g. constants::p_mass

// ===========================================================================
// Return flags
// ===========================================================================
// Every solver function returns 0 on success. A negative integrator flag
// (FLAG_TOO_MUCH_WORK = -1, FLAG_TOO_MUCH_ACC = -2, FLAG_ERR_FAILURE = -3,
// FLAG_CONV_FAILURE = -4, ...) means the integration failed after the
// restarts; the input state is then left untouched and, by default, a
// KROME_ERROR_REPORT is written (SolverOpts::error_report). The
// integrator's failures (-1 .. -27, CVODE's values: KROME's integrator is
// a port of CVODE's BDF, core/cvode/bdf_cell.hpp):
constexpr int FLAG_TOO_MUCH_WORK = -1;            // max_steps exhausted before t_end
constexpr int FLAG_TOO_MUCH_ACC = -2;             // tolerances too small for machine precision
constexpr int FLAG_ERR_FAILURE = -3;              // repeated error-test failures
constexpr int FLAG_CONV_FAILURE = -4;             // repeated Newton convergence failures
constexpr int FLAG_LSETUP_FAIL = -6;              // Newton matrix setup / factorization failed
constexpr int FLAG_LSOLVE_FAIL = -7;              // linear solve failed
constexpr int FLAG_RHSFUNC_FAIL = -8;             // right-hand side failed unrecoverably
constexpr int FLAG_FIRST_RHSFUNC_ERR = -9;        // right-hand side failed at the initial state
constexpr int FLAG_REPTD_RHSFUNC_ERR = -10;       // repeated recoverable right-hand side errors
constexpr int FLAG_UNREC_RHSFUNC_ERR = -11;       // right-hand side failed unrecoverably in a restart
constexpr int FLAG_NLS_FAIL = -16;                // nonlinear solver failure
constexpr int FLAG_ILL_INPUT = -22;               // illegal input (e.g. error weights not positive)
constexpr int FLAG_TOO_CLOSE = -27;               // t_end too close to the start
// The KROME flags below (<= -100) mean the inputs were rejected before
// integrating, except FLAG_NONPOSITIVE_T: the integration reached Tgas <=
// 0 (trial states there are refused, so the integrator shortens its step;
// if it cannot avoid them, or the final Tgas is <= 0, the solve fails with
// this flag).
constexpr int FLAG_OK = 0;
constexpr int FLAG_NO_DUST = -102;         // cooling DUSTBINS, Tgas evolved,
                                           // and no SolverOpts::dust
constexpr int FLAG_NO_FIELD = -103;        // photo reactions without a
                                           // field (SolverOpts::photo/Phi),
                                           // or radiation.evolve without a
                                           // PhotoModule
constexpr int FLAG_UNSUPPORTED = -104;     // the requested path (device)
                                           // does not support this
                                           // configuration or is not built
constexpr int FLAG_PARAM_UNSET = -105;     // a required parameter is NaN
constexpr int FLAG_CONSERVE_GAS_TOTAL = -106;  // SolverOpts::conserve with
                                           // the dust.gas_total representation
constexpr int FLAG_SURFACE_BINS = -107;    // surface.bins = NB without NB
                                           // dust bins (SolverOpts::dust)
constexpr int FLAG_NO_LAPACK = -108;       // LINSOL_LAPACK in a build
                                           // without KROME_USE_LAPACK
constexpr int FLAG_NO_KOKKOS = -109;       // device path before
                                           // Kokkos::initialize()
constexpr int FLAG_NONPOSITIVE_T = -110;   // the integration reached a
                                           // non-positive (or non-finite)
                                           // Tgas: inputs untouched
constexpr int FLAG_DEVICE_MEMORY = -111;   // PATH_DEVICE: not even one cell
                                           // fits in the device memory

// a one-line description of a flag, never nullptr
inline const char* flag_message(int flag);

// ===========================================================================
// Solver options
// ===========================================================================

// Jacobian of the Newton iteration
enum JacMode {
  JAC_AUTO,      // analytic at fixed Tgas, difference quotients when Tgas
                 // is evolved (the robust choice there: JAC_ANALYTIC with
                 // Tgas evolved is faster for large networks, but not yet
                 // robust at tiny atol or with evolving dust)
  JAC_ANALYTIC,  // the generated analytic chemistry Jacobian (difference
                 // quotients only for the thermal and auxiliary rows)
  JAC_FD         // difference-quotient Jacobian (CVODE's; legacy MF=222)
};

// linear solver of the Newton systems
enum LinSol {
  LINSOL_AUTO,    // LAPACK from LAPACK_MIN_NEQ equations when built with
                  // it, else the integrator's own LU
  LINSOL_DENSE,   // the integrator's own dense LU (SUNDIALS', ported)
  LINSOL_LAPACK,  // LAPACK (vendor BLAS); FLAG_NO_LAPACK without it
  LINSOL_SPARSE   // a sparse LU with static pivots in the generated
                  // fill-reducing order (bdf_cell.hpp Opts::sparse), the
                  // Tgas / aux equations as a dense border
};
constexpr int LAPACK_MIN_NEQ = 64;
#ifdef KROME_USE_LAPACK
constexpr bool HAS_LAPACK = true;
#else
constexpr bool HAS_LAPACK = false;
#endif

// what to evolve and how (all numerics: the physics is fixed at
// generation). The defaults are the recommended settings.
struct SolverOpts {
  // evolve Tgas with the chemistry (the thermal processes chosen at
  // generation); false: fixed-Tgas chemistry at the input temperature
  bool evolve_T = false;
  JacMode jac = JAC_AUTO;
  LinSol linsol = LINSOL_AUTO;
  // tolerances: relative, and absolute per state block. rtol and atol_n
  // start from the project's defaults (codegen solver.rtol / solver.atol, 1e-8 and
  // 1e-20 unless given)
  double rtol = DEFAULT_RTOL;
  // (atol_n also bounds the key-element atoms locked in each dust bin,
  // the dust.evolve state)
  double atol_n = DEFAULT_ATOL_N;  // species [cm^-3]
  double atol_T = 1e-6;      // Tgas and the dust temperatures [K]
  double atol_u = 1e-10;     // ln Phi (radiation.evolve)
  // per-species tolerances (legacy -customATOL / -customRTOL): arrays of
  // NSPEC values replacing atol_n / rtol for the species (the other state
  // slots keep the scalars); nullptr, or a negative entry, takes the
  // scalar. The arrays are the caller's (not copied) and start from the
  // project's defaults (codegen solver.atol_file / solver.rtol_file: the listed
  // species, -1 for the others; none unless given). A per-species rtol
  // runs on every path
  const double* atol_species = DEFAULT_ATOL_SPECIES_PTR;
  const double* rtol_species = DEFAULT_RTOL_SPECIES_PTR;
  // integrator limits: the step budget (exhausting it is final:
  // FLAG_TOO_MUCH_WORK), and the failures tolerated within one step
  long max_steps = 100000;
  int max_err_fails = 7;
  int max_conv_fails = 10;
  // after an error-test, convergence or repeated rhs failure, restart
  // the integrator at the reached state (fresh history and Jacobian, first
  // step a tenth of the one it had reduced to), up to this many times
  // (legacy re-entered DLSODES the same way)
  int max_restarts = 5;
  // n = max(n, 0) after a successful solve (legacy)
  bool clip_negative = true;
  // mass-fraction solves (solve_x, solve_cells_x): rescale the returned
  // fractions to their input sum (legacy -useX forced mass conservation
  // this way); off: x = m n / rho as integrated
  bool renormalize_x = false;
  // post-solve correction of numerical drift: 0 off, 1 the electrons
  // restore the input total charge, 2 the minimal relative correction
  // restoring the input element totals and charge, 3 the electrons
  // neutralize the ions (even when frozen; legacy -computeElectrons)
  int conserve = 0;
  // on a failed integration, write a readable report (legacy
  // KROME_ERROR_REPORT: flag, state, dn/dt, sorted rates and fluxes, the
  // input state as C++/Fortran initializers, integrator statistics) to
  // error_file, replacing the previous one. In an MPI run give each rank
  // its own file name.
  bool error_report = true;
  const char* error_file = "KROME_ERROR_REPORT";
  // photon bins (networks with photo reactions): the PhotoModule built
  // from the bins (see Radiation / PhotoModule below) and the photon flux
  // per bin Phi[photo->B], held fixed during the solve, or co-evolved and
  // written back when the network was generated with radiation.evolve
  const PhotoModule* photo = nullptr;
  double* Phi = nullptr;
  // binned dust (cooling DUSTBINS, surface.bins, dust.evolve_Td): the size
  // distribution and per-bin temperatures (see DustBins below); with
  // dust.evolve_Td the Td are the initial state and are written back
  DustBins* dust = nullptr;
};

// what a solve did (optional output of solve() / solve_x())
struct SolveStats {
  long steps = 0;          // integration steps
  long rhs_evals = 0;      // right-hand-side evaluations
  long jac_evals = 0;      // Jacobian evaluations
  long rhs_evals_jac = 0;  // rhs evaluations of the DQ Jacobian (JAC_FD)
  long lin_setups = 0;     // Newton-matrix factorizations
  long newton_iters = 0;   // Newton iterations
  long newton_fails = 0;   // Newton convergence failures
  long err_fails = 0;      // step attempts failing the error test
  int restarts = 0;        // restarts after failures (max_restarts)
  int flag = 0;            // the returned flag
};

// many cells in one call (solve_cells)
enum BatchPath {
  PATH_AUTO,     // PER-CELL, or DEVICE from min_device cells when a device
                 // backend is compiled in and Kokkos is initialized
  PATH_PERCELL,  // an (OpenMP-parallel) loop of independent solve() calls
  PATH_DEVICE    // the Kokkos backend (-DKROME_USE_KOKKOS): all cells in one
                 // call to the device (GPU or OpenMP), one thread per cell,
                 // each with its own integration to its own dt; cells that
                 // fail there are re-solved on the host
};

struct CellsOpts {
  BatchPath path = PATH_AUTO;
  // keep the device buffers between calls (sized for the largest call so
  // far; memory only: every call writes its own inputs). Off: allocated
  // and freed in every call. kk::release_device_cache() frees them;
  // Kokkos::finalize does
  bool device_cache = true;
  // the device memory a call may use [bytes]: calls of more cells are
  // integrated in chunks that fit (on the same buffers). 0: 80 % of the
  // device's free memory (the cached buffers the call reuses included),
  // on CUDA and HIP (elsewhere unbounded); < 0: unbounded. An allocation
  // that fails anyway halves the chunk; when not even one cell fits,
  // PATH_AUTO runs the cells on the host and PATH_DEVICE returns
  // FLAG_DEVICE_MEMORY
  long long device_memory = 0;
  // the device path's memory layout: interleaved across cells (element i
  // of cell c at c + i ncell: coalesced on a GPU) for calls (or chunks)
  // of at least this many cells, contiguous per cell below; < 0:
  // contiguous, on CUDA interleaved above one wave of the contiguous
  // kernel (the cells resident at once: 27 648 on an A100). Numerics only: the layouts agree to roundoff (bitwise
  // without FMA contraction)
  int min_interleave = -1;
  // AUTO: calls from this many cells go to the device backend, when
  // compiled in (smaller ones run per cell on the CPU, faster below ~1e4
  // cells: an A100 matched ~20 CPU cores per cell only from ~1e4 cells)
  int min_device = 10000;
  // per-cell photon fields Phi[ncell][photo->B] (instead of
  // SolverOpts::Phi): read, or updated in place with radiation.evolve.
  // Without them every cell starts from SolverOpts::Phi; with
  // radiation.evolve a single-cell call then updates it, and in a call of
  // several cells the evolved fields are not returned
  double* Phi = nullptr;
  // per-cell dust temperatures Td[ncell][dust->nb] [K] (with
  // SolverOpts::dust): each cell's bins start from its row (with
  // dust.evolve_Td the initial state) and the row receives the cell's final
  // Td (the balance at the final state, or the evolved state); a failed
  // cell's row is untouched. nullptr: every cell starts from dust->Td and
  // the results are discarded
  double* Td = nullptr;
  // per-cell grain radii adust[ncell][dust->nb] [cm] (dust.evolve): each
  // cell starts from its row and the row receives the evolved radii (0:
  // an evaporated bin); a failed cell's row is untouched. nullptr: every
  // cell starts from dust->a and the results are discarded. Hydro codes
  // advect them with the gas (the grain counts are d2g x the population's
  // grains per gas mass)
  double* adust = nullptr;
  // optional per-cell return flags flags[ncell]
  int* flags = nullptr;
  // optional per-cell solver statistics stats[ncell] (a cell re-solved on
  // the host gets the host solve's)
  SolveStats* stats = nullptr;
};

// the steady-state solver (solve_equilibrium)
struct EquilOpts {
  double t0 = 100.0 * constants::seconds_per_year;   // first output [s]
  double t_max = 5e8 * constants::seconds_per_year;  // give up after [s]
  double factor = 3.0;             // output times t0 * factor^k
  double ertol = 1e-5;             // converged when the relative change
  double eatol = 1e-12;            // between outputs, with n floored at
                                   // eatol * ntot, is below ertol
};

// ===========================================================================
// Integration
// ===========================================================================

// Integrate one cell from t = 0 to t_end [s].
//   n[NSPEC]  number densities [cm^-3], updated in place
//   Tgas      [K]: updated with o.evolve_T, else the fixed temperature
//   params    the cell's parameters
//   o         options; o.Phi is updated with radiation.evolve, o.dust->Td with
//             dust.evolve_Td
//   stats     optional counters
// Returns 0, a KROME flag (inputs rejected) or a negative integrator flag
// (failed after the restarts). On failure n, Tgas, Phi and Td are left
// untouched.
inline int solve(double* n, double& Tgas, const Params& params,
                 double t_end, const SolverOpts& o = SolverOpts(),
                 SolveStats* stats = nullptr);

// solve() with mass fractions x[NSPEC] (legacy -useX) and the gas density
// rho [g cm^-3]: n = rho x / m on the way in, x = m n / rho with the same
// rho on the way out (no renormalization). The solver integrates number
// densities. On failure x is untouched.
inline int solve_x(double* x, double rho, double& Tgas, const Params& params,
                   double t_end, const SolverOpts& o = SolverOpts(),
                   SolveStats* stats = nullptr);

// Integrate ncell independent cells, each over its own dt[c] [s]:
//   n[ncell * NSPEC], Tgas[ncell], dt[ncell]
//   p: one Params for every cell, or params[ncell] (per-cell overload)
//   co: path, per-cell fields and per-cell flags (co.flags[ncell])
// Returns 0, or the flag of the lowest-index failing cell. A failed cell
// is left untouched; the others are updated. Cells that fail on the
// device are re-solved on the host, so only cells that really fail
// report.
inline int solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                       const Params& p, const SolverOpts& o,
                       const CellsOpts& co = CellsOpts());
inline int solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                       const Params* params, const SolverOpts& o,
                       const CellsOpts& co = CellsOpts());

// solve_cells() with mass fractions x[ncell * NSPEC] and the gas density
// rho[ncell] (see solve_x); a failed cell keeps its input x
inline int solve_cells_x(int ncell, double* x, const double* rho,
                         double* Tgas, const double* dt, const Params& p,
                         const SolverOpts& o,
                         const CellsOpts& co = CellsOpts());
inline int solve_cells_x(int ncell, double* x, const double* rho,
                         double* Tgas, const double* dt,
                         const Params* params, const SolverOpts& o,
                         const CellsOpts& co = CellsOpts());

// Integrate to steady state (legacy krome_equilibrium): outputs at
// t0 * factor^k, each interval a solve() from the current state;
// converged when the relative change between outputs is below e.ertol
// (and |dT|/T when Tgas evolves). Isothermal unless o.evolve_T.
// Returns 0 (converged), +1 (t_max reached first: the state at t_max is
// returned) or a negative flag (inputs untouched). t_reached (optional):
// the time integrated [s].
inline int solve_equilibrium(double* n, double& Tgas, const Params& params,
                             const SolverOpts& o = SolverOpts(),
                             const EquilOpts& e = EquilOpts(),
                             double* t_reached = nullptr);

#ifdef KROME_USE_KOKKOS
namespace kk {
// the device backend directly (solve_cells with PATH_DEVICE is the usual
// entry, and falls back to the host when this returns FLAG_UNSUPPORTED or
// a failure). Needs Kokkos::initialize(). pstride: 0 = params[0] for
// every cell, 1 = params[ncell]. Returns 0, FLAG_UNSUPPORTED,
// FLAG_PARAM_UNSET, FLAG_NO_KOKKOS, or the flag of the first cell that
// failed on the device and again on the host (the others written).
inline int solve_cells_kokkos(int ncell, double* n, double* Tgas,
                              const double* dt, const Params* params,
                              int pstride, const SolverOpts& o,
                              const CellsOpts& co = CellsOpts());
}  // namespace kk
#endif

// ===========================================================================
// Queries: the cell as the solver sees it (temperature capped to
// [Tcmb, 1e9 K] for the rates, to [2.73 K, 1e9 K] for the thermal terms
// and gamma; o.photo / o.Phi / o.dust as in a solve). Functions
// returning int return 0, or FLAG_NO_FIELD / FLAG_PARAM_UNSET /
// FLAG_NO_DUST, or -1 when the network has no photon bins (get_opacity,
// get_absorption).
// ===========================================================================

// rate coefficients k[NREA] (including the photo rates)
inline int get_rates(const double* n, double Tgas, const Params& p,
                     const SolverOpts& o, double* k);

// reaction fluxes k * prod(n) [cm^-3 s^-1], fl[NREA]
inline int get_fluxes(const double* n, double Tgas, const Params& p,
                      const SolverOpts& o, double* fl);

// the dust temperature the rates see [K] (legacy krome_get_table_Tdust):
// the dust.tables value, else the Tdust parameter (the binned-dust mean
// with o.dust); -1 when the network has no dust temperature (the binned
// dust's own per-bin Td: DustBins::balance, krome_opts_get_dust_Td)
inline int get_Tdust(const double* n, double Tgas, const Params& p,
                     const SolverOpts& o, double* Td);

// thermal balance [erg cm^-3 s^-1]: total heating and cooling, the
// configured terms one by one (heats[] in heating_name order, cools[] in
// cooling_name order) and dTgas/dt [K s^-1]; every output may be nullptr
inline int get_thermal(const double* n, double Tgas, const Params& p,
                       const SolverOpts& o, double* heat, double* cool,
                       double* heats = nullptr, double* cools = nullptr,
                       double* dTdt = nullptr);

// names of the configured heating / cooling terms (index i of heats[] /
// cools[]), nullptr beyond the last
inline const char* heating_name(int i);
inline const char* cooling_name(int i);

// adiabatic index at (n, Tgas): gamma_ad, or the thermal.gamma model
inline double get_gamma(const double* n, double Tgas, const Params& p);
// mean molecular weight rho / (ntot m_p)
inline double get_mu(const double* n);
// mass density [g cm^-3], total number density [cm^-3], H nuclei [cm^-3]
inline double get_rho(const double* n);
inline double get_ntot(const double* n);
inline double get_Hnuclei(const double* n);
// net charge density sum q n [e cm^-3]; the electron density that would
// make the cell neutral, max(sum over ions q n, 0) (legacy get_electrons)
inline double get_charge(const double* n);
inline double get_electrons_neutral(const double* n);
// element totals tot[NELEM] [atoms cm^-3], in element_names order
inline void get_element_totals(const double* n, double* tot);
// free-fall time sqrt(3 pi / (32 G rho)) [s], from n or from the mass
// density rho [g cm^-3]; Jeans length [cm]
inline double get_free_fall_time(const double* n);
inline double get_free_fall_time_rho(double rho);
inline double get_jeans_length(const double* n, double Tgas);
// the column density [cm^-2] of a number density num [cm^-3] in the cell,
// with the network's estimator (radiation.num2col; legacy krome_num2col)
inline double num2col(double num, const double* n, double Tgas,
                      const Params& p);

// the line coolants: the coolZ.dat metals of -cooling=Z, then the
// thermal.coolant_files coolants (NCOOLANT; names coolant_name(i), legacy CII,
// FeII, O2; coolant_species[i]; levels coolant_nlev[i], energies [K] and
// degeneracies coolant_level_E / _g[coolant_lev0[i] + l]).
// get_coolant_cooling: each coolant's line cooling cool[NCOOLANT]
// [erg cm^-3 s^-1] (legacy krome_cooling<metal>; the Z and coolant-file
// terms are their sums). get_coolant_populations: the level populations
// of coolant i [cm^-3] in the file's numbering (legacy
// krome_popcool_dump)
inline const char* coolant_name(int i);
inline int get_coolant_cooling(const double* n, double Tgas, const Params& p,
                               const SolverOpts& o, double* cool);
inline int get_coolant_populations(const double* n, double Tgas,
                                   const Params& p, const SolverOpts& o,
                                   int i, double* pop);

// photoheating per photo reaction heat[NPHOTO] [erg cm^-3 s^-1], summing to
// the PHOTO heating term (legacy krome_get_photoBin_heats gave it per
// absorber); the field over the Draine ISRF (legacy
// krome_get_ratioFluxDraine). -1 without photon bins
inline int get_photo_heating(const double* n, double Tgas, const Params& p,
                             const SolverOpts& o, double* heat);
inline int get_ratio_flux_draine(const SolverOpts& o, double* ratio);

// ---- legacy conveniences (krome_user) ------------------------------------
// the species names separated by spaces, for an output header (legacy
// krome_get_names_header)
inline std::string names_header();
// the nbest largest reaction fluxes of the cell, one line each: the
// 1-based reaction index, its equation and the flux [cm^-3 s^-1] (legacy
// krome_print_best_flux); returns 0 or a query flag
inline int print_best_flux(const double* n, double Tgas, const Params& p,
                           int nbest, const SolverOpts& o = SolverOpts(),
                           std::FILE* f = stdout);
// every reaction's flux, absolute and relative to the largest and to the
// sum, one line per reaction then a blank line, labelled by xvar (e.g. the
// density; legacy krome_explore_flux)
inline int explore_flux(std::FILE* f, const double* n, double Tgas,
                        const Params& p, double xvar,
                        const SolverOpts& o = SolverOpts());
// the species table: index, name, mass [g], charge, abundance, their sum,
// and Tgas (legacy krome_get_info)
inline void print_info(const double* n, double Tgas, std::FILE* f = stdout);
// (scale_Z(n, Z), legacy krome_scale_Z, is generated with the network:
// its neutral atomic metals at the Asplund+09 abundances times 10^Z)
// every coolant's level populations, one row per level: name, level
// (1-based), Tgas, population [cm^-3], the coolant's total (legacy
// krome_popcool_dump)
inline int print_coolant_populations(std::FILE* f, const double* n,
                                     double Tgas, const Params& p,
                                     const SolverOpts& o = SolverOpts());
// one row: Tgas, the total cooling, every cooling term (legacy
// krome_dump_cooling)
inline int dump_cooling(std::FILE* f, const double* n, double Tgas,
                        const Params& p, const SolverOpts& o = SolverOpts());
// the rate coefficients on npts log temperatures Tmin..Tmax at n = 1e-40:
// rows Tgas, reaction (1-based), k (legacy krome_dump_rates)
inline int dump_rates(std::FILE* f, double Tmin, double Tmax, int npts,
                      const Params& p, const SolverOpts& o = SolverOpts());
// every parameter, 'name: value' (legacy krome_print_phys_variables)
inline void print_params(const Params& p, std::FILE* f = stdout);

// optical depth per photon bin tau[B] of the cell's columns (num2col
// estimate, legacy krome_get_opacity): photo reactions and J_LW channels,
// plus dust when d2g > 0 (the bins of o.dust, else kdust[B] [cm^2 per g of
// dust], may be nullptr)
inline int get_opacity(const double* n, double Tgas, const Params& p,
                       const SolverOpts& o, double d2g, const double* kdust,
                       double* tau);

// photon removal coefficient per bin kabs[B] [cm^-1] (the radiation.evolve
// sink, plus the binned dust of o.dust, plus rho * kdust[B] when kdust is
// given as the dust opacity per GAS mass [cm^2 g^-1]); times a length it
// is legacy krome_get_opacity_size
inline int get_absorption(const double* n, double Tgas, const Params& p,
                          const SolverOpts& o, const double* kdust,
                          double* kabs);

// energy released per reaction event q[NREA] [MeV] (the qeff= annotations,
// 0 elsewhere); times get_fluxes: legacy krome_stars_energy
inline int get_qeff(const double* n, double Tgas, const Params& p,
                    double* q);

// number densities <-> mass fractions: x = m n / rho(n); n = rho x / m
inline void n2x(const double* n, double* x);
inline void x2n(const double* x, double rho, double* n);

// ===========================================================================
// Photon bins and binned dust (setup objects, defined in core/)
// ===========================================================================
//
// struct Radiation                       the photon bins of a run
//   std::vector<double> Elo, Ehi, Phi    bin edges [eV], flux per bin
//   static Radiation log_bins(B, emin, emax)
//
// struct PhotoModule                     bin-averaged cross sections of the
//                                        network's photo reactions
//   void init(const Radiation& rad)      build them for these bins (once;
//                                        then o.photo = &module)
//   int B                                number of bins
//   double habing(const double* Phi)     Habing flux G0 of a field
//
// struct DustBins                        MRN size distribution per grain
//                                        type, per-bin Qabs and Tdust
//   int nd = 10; bool useSi, useC;       configuration (bins per type,
//   double alow, aup, phi;               grain types, MRN range [cm] and
//                                        slope), then:
//   void init()                          build the bins
//   void interp_qabs(const Radiation&)   Qabs at the photon bins
//   std::vector<double> Td               per-bin dust temperatures [K]
//
// Build the fields for Phi with the spectra:: functions below.

namespace spectra {
// bin edges Elo[B], Ehi[B] [eV]: linear / logarithmic in [lower, upper]
inline void bins_lin(int B, double lower, double upper, double* Elo,
                     double* Ehi);
inline void bins_log(int B, double lower, double upper, double* Elo,
                     double* Ehi);
// a specific intensity J [eV s^-1 cm^-2 Hz^-1 sr^-1] constant across a
// bin <-> the bin's photon flux Phi = 4 pi J ln(Ehi/Elo) / h
inline double J_to_phi(double J, double Elo, double Ehi);
inline double phi_to_J(double phi, double Elo, double Ehi);
// blackbody B_nu(Tbb) [eV s^-1 cm^-2 Hz^-1 sr^-1] at E [eV]
inline double planck_J(double E, double Tbb);
// the fields phi[B], each integrated exactly over its bin:
//   blackbody at Tbb (no dilution)
inline void blackbody(int B, const double* Elo, const double* Ehi,
                      double Tbb, double* phi);
//   Draine 1978 ISRF (5-13.6 eV)
inline void draine(int B, const double* Elo, const double* Ehi,
                   double* phi);
//   J = J21 1e-21 (13.6 eV / E)^alpha erg s^-1 cm^-2 Hz^-1 sr^-1
inline void j21(int B, const double* Elo, const double* Ehi, double J21,
                double alpha, double* phi);
//   a table J(E) at nE ascending energies [eV], linear or log-log between
//   the nodes, zero outside
inline void table(int B, const double* Elo, const double* Ehi, int nE,
                  const double* Et, const double* Jt, bool loglog,
                  double* phi);
//   Haardt & Madau 2012 UV background at redshift z; -1 outside
//   z in [0, 15.66]
inline int hm2012(int B, const double* Elo, const double* Ehi, double z,
                  double* phi);
// the energy range where B_nu(Tbb) > eps (legacy BBlog_auto); -1 if none
inline int blackbody_range(double Tbb, double* lower, double* upper,
                           double eps = 1e-6);
// bin-mean dust opacity kbin[B] [cm^2 g^-1]: from a table kappa(E) at nE
// ascending energies [eV] (-1 when a bin is outside the table), or from
// the embedded Draine R_V = 3.1 curve
inline int dust_kappa_bins(int nE, const double* Et, const double* kt,
                           int B, const double* Elo, const double* Ehi,
                           double* kbin);
inline int dust_kappa_draine(int B, const double* Elo, const double* Ehi,
                             double* kbin);
// G0 and Av that best explain a field as an attenuated Draine field
// (legacy krome_find_G0_Av); xH: H mass fraction. -1 with fewer than 2
// usable bins
inline int find_G0_Av(int B, const double* Elo, const double* Ehi,
                      const double* phi, const double* kbin, double d2g,
                      double xH, double* G0, double* Av);
}  // namespace spectra

}  // namespace krome

// ---- implementation -------------------------------------------------------
#include "core/solve_cell.hpp"
#include "core/solve_cells.hpp"
#include "core/query.hpp"
#include "core/spectra.hpp"
#ifdef KROME_USE_KOKKOS
#include "core/solve_kokkos.hpp"
#endif
