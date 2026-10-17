// KROME-X solver: one BDF integrator (core/cvode/bdf_cell.hpp, CVODE's BDF
// ported) for every combination of evolved ingredients (dense LU, or
// LAPACK's; the per-cell device backend runs the same integrator).
//
//   state y = [ n_0..n_{S-1} | Tgas     | u_0..u_{B-1} | Td_0..Td_{D-1} ]
//               always         evolve_T   EvolveConfig   EvolveConfig
//                              (runtime)  ::rad          ::Td (with Tgas)
//   (the photon-field and dust-temperature blocks are fixed at generation,
//   codegen radiation.evolve / dust.evolve_Td; Tgas is the one runtime switch)
//   u_b = ln Phi_b: the co-evolved photon flux per bin cannot overshoot
//   below zero (docs/solver-design.md section 2)
//
//   dn/dt   = chemistry(n, k(T, n, Phi)) [+ secondary-ionization fluxes]
//   dT/dt   = (gamma_ad - 1)(heat - cool) / (kB ntot)       (thermal_aux)
//   du_b/dt = -c kabs_b(n)                                 (pure absorption)
//   dTd/dt  = pre dT/dt / (pre + intBB'(Td) besc)  (legacy -usedTdust)
//
// Blocks absent from the state are inputs: a fixed Tgas, a field Phi held
// at the caller's values (photochemistry, photoheating and grains still
// use it), and per-bin Tdust from the radiative-balance bisection.
//
// Jacobian (SolverOpts::jac):
//   JAC_FD        the integrator's difference quotients (CVODE's dense
//                 DQ Jacobian; the analog of legacy DLSODES MF=222)
//   JAC_ANALYTIC  generated analytic chemistry block + analytic radiation
//                 coupling blocks; difference quotients for the Tgas/Tdust
//                 columns and for the thermal, secondary-ionization and
//                 dust-temperature ("aux") rows
//   JAC_AUTO      (default) ANALYTIC at fixed Tgas, FD when Tgas is evolved.
//                 With Tgas evolved, ANALYTIC is the hybrid (difference
//                 quotients for the Tgas column and the thermal rows): 1.6x
//                 faster than FD on DarkCloud with cooling, 12-18 % slower on
//                 10-28 equations, but not robust yet (2026-10-06): its aux
//                 rows (thermal, Td) vs trace species use fd_step_n_, whose
//                 1e-9 ntot floor is far above a trace species' value
//                 (evolving dust: ~100x the Jacobians between 160 and 2700
//                 K), while a step relative to n_j alone falls below the
//                 row's roundoff (legacy collapseCO at atol 1e-40: d(Tgas)/
//                 d(H2) = 0, Newton fails every step). The fix: CVODE's
//                 increments (error weights, step size, rhs norm) there
//
// Implementation header: include krome.hpp (the public API), which
// declares solve(), solve_x(), solve_equilibrium() and the option structs.
#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <mutex>
#include <string>
#include <limits>
#include <utility>
#include <vector>
// the integrator: KROME's port of CVODE's BDF (LAPACK's dense LU with
// KROME_USE_LAPACK, the project Makefile's LAPACK= option)
#include "cvode/bdf_cell.hpp"

namespace krome {

constexpr double KB_ERG = constants::boltzmann_erg;  // legacy KROME value

inline const char* flag_message(int flag) {
  switch (flag) {
    case FLAG_OK: return "success";
    case FLAG_NO_DUST: return "cooling DUSTBINS with Tgas evolved needs SolverOpts::dust";
    case FLAG_NO_FIELD: return "photo reactions need a PhotoModule and a field Phi (or --evolve-rad without a PhotoModule)";
    case FLAG_UNSUPPORTED: return "configuration not supported by the requested path, or backend not built";
    case FLAG_PARAM_UNSET: return "a required parameter is unset (NaN)";
    case FLAG_CONSERVE_GAS_TOTAL: return "conservation is not available with --gas-total";
    case FLAG_SURFACE_BINS: return "--surface-bins needs SolverOpts::dust with the generated number of bins";
    case FLAG_NO_LAPACK: return "LINSOL_LAPACK in a build without LAPACK (KROME_USE_LAPACK)";
    case FLAG_NO_KOKKOS: return "device path used before Kokkos::initialize()";
    case FLAG_NONPOSITIVE_T: return "the integration ended at a non-positive or non-finite Tgas (input untouched)";
    case FLAG_DEVICE_MEMORY: return "device path: not even one cell fits in the device memory (CellsOpts::device_memory)";
    case FLAG_TOO_MUCH_WORK: return "max_steps exhausted before t_end (FLAG_TOO_MUCH_WORK)";
    case FLAG_TOO_MUCH_ACC: return "tolerances too small for machine precision (FLAG_TOO_MUCH_ACC)";
    case FLAG_ERR_FAILURE: return "repeated error-test failures (FLAG_ERR_FAILURE)";
    case FLAG_CONV_FAILURE: return "repeated Newton convergence failures (FLAG_CONV_FAILURE)";
    case FLAG_LSETUP_FAIL: return "Newton matrix setup / factorization failed (FLAG_LSETUP_FAIL)";
    case FLAG_LSOLVE_FAIL: return "linear solve failed (FLAG_LSOLVE_FAIL)";
    case FLAG_RHSFUNC_FAIL: return "right-hand side failed (FLAG_RHSFUNC_FAIL)";
    case FLAG_FIRST_RHSFUNC_ERR: return "right-hand side failed at the initial state (FLAG_FIRST_RHSFUNC_ERR)";
    case FLAG_REPTD_RHSFUNC_ERR: return "repeated recoverable right-hand side errors (FLAG_REPTD_RHSFUNC_ERR)";
    case FLAG_UNREC_RHSFUNC_ERR: return "right-hand side failed unrecoverably (FLAG_UNREC_RHSFUNC_ERR)";
    case FLAG_ILL_INPUT: return "illegal input (FLAG_ILL_INPUT)";
    case FLAG_NLS_FAIL: return "nonlinear solver failure (FLAG_NLS_FAIL)";
    case FLAG_TOO_CLOSE: return "t_end too close to the start (FLAG_TOO_CLOSE)";
    default: break;
  }
  if (flag > 0)
    return "not converged: t_max reached first (solve_equilibrium)";
  return flag <= -100 ? "unknown KROME flag" : "integrator failure";
}

// the name of a flag (the error report)
inline const char* flag_name(int flag) {
  switch (flag) {
    case FLAG_OK: return "FLAG_OK";
    case FLAG_NO_DUST: return "FLAG_NO_DUST";
    case FLAG_NO_FIELD: return "FLAG_NO_FIELD";
    case FLAG_UNSUPPORTED: return "FLAG_UNSUPPORTED";
    case FLAG_PARAM_UNSET: return "FLAG_PARAM_UNSET";
    case FLAG_CONSERVE_GAS_TOTAL: return "FLAG_CONSERVE_GAS_TOTAL";
    case FLAG_SURFACE_BINS: return "FLAG_SURFACE_BINS";
    case FLAG_NO_LAPACK: return "FLAG_NO_LAPACK";
    case FLAG_NO_KOKKOS: return "FLAG_NO_KOKKOS";
    case FLAG_NONPOSITIVE_T: return "FLAG_NONPOSITIVE_T";
    case FLAG_DEVICE_MEMORY: return "FLAG_DEVICE_MEMORY";
    case FLAG_TOO_MUCH_WORK: return "FLAG_TOO_MUCH_WORK";
    case FLAG_TOO_MUCH_ACC: return "FLAG_TOO_MUCH_ACC";
    case FLAG_ERR_FAILURE: return "FLAG_ERR_FAILURE";
    case FLAG_CONV_FAILURE: return "FLAG_CONV_FAILURE";
    case FLAG_LSETUP_FAIL: return "FLAG_LSETUP_FAIL";
    case FLAG_LSOLVE_FAIL: return "FLAG_LSOLVE_FAIL";
    case FLAG_RHSFUNC_FAIL: return "FLAG_RHSFUNC_FAIL";
    case FLAG_FIRST_RHSFUNC_ERR: return "FLAG_FIRST_RHSFUNC_ERR";
    case FLAG_REPTD_RHSFUNC_ERR: return "FLAG_REPTD_RHSFUNC_ERR";
    case FLAG_UNREC_RHSFUNC_ERR: return "FLAG_UNREC_RHSFUNC_ERR";
    case FLAG_NLS_FAIL: return "FLAG_NLS_FAIL";
    case FLAG_ILL_INPUT: return "FLAG_ILL_INPUT";
    case FLAG_TOO_CLOSE: return "FLAG_TOO_CLOSE";
    default: return "?";
  }
}

// the thermal processes are FIXED at generation time: the generated
// ThermalConfig (codegen thermal.cooling / thermal.heating / thermal.floor) holds
// them as constexpr flags, so disabled processes are compiled out and
// the caller cannot switch individual processes. At runtime only thermal
// evolution as a whole is optional (SolverOpts::evolve_T).
using TC = ThermalConfig;

// what the thermal terms see: parameters and the optional photo field /
// binned dust
struct ThermCtx {
  Params p;
#ifdef KROME_PHOTO_HPP_INCLUDED
  // the photo module: host-side object (setup, Jacobian blocks, queries)
  // and the raw view of its tables the kernels read (host or device data)
  const PhotoModule* photo = nullptr;
  bool photo_on = false;
  PhotoView ph;
  const double* Phi = nullptr;
  // the bin sums of a field held fixed during the solve (solve() computes
  // them once; nullptr: evaluate from Phi)
  const PhotoSumsView* psum = nullptr;
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  // the binned dust (dust_on): the population as a view (during a solve
  // its radii are the state's, dust.acur, and its counts those of the
  // solve start, dust.xds_fix), the population's own Td (the input; fixed
  // without dust.evolve_Td state or balance), and the per-bin balance Td at
  // the current state (computed with the rates; nullptr with dust.evolve_Td,
  // whose Td are state)
  bool dust_on = false;
  DustView dust;
  const double* dust_Td = nullptr;
  const double* dust_tdb = nullptr;
#endif
  double k[NREA + 1];  // +1 keeps zero-reaction (cooling-only) networks legal
  // radiative channels (spec 5.4), evaluated with the rates: photo
  // shielding factors and the J_LW unit rates
  double shield[NPHOTO + 1];
  double khat[NLW + 1];
};

// runtime parameters a configured process reads: templates, so a network
// without that process (hence without that Params member) never
// instantiates them (the calls below sit in discarded if-constexpr
// branches)
template <class P>
KROME_HD inline double heat_compress_(double Tgas, const double* n, const P& p) {
  return ntot_(n) * KB_ERG * Tgas / p.tff;  // legacy heat_compress
}
template <class P>
KROME_HD inline double heat_visc_(const double* n, const P& p) {
  return 2.25 * p.visc_nu * rho_(n) * p.visc_omega * p.visc_omega;
}
template <class P>
KROME_HD inline double habing_param_(const P& p) {
  return p.G0 * std::exp(-2.5 * p.Av);  // legacy exp(-2.5 Av)
}
template <class P>
KROME_HD inline double cooling_floor_T_(const P& p) {
  return p.Tfloor;
}

// the individual thermal terms [erg cm^-3 s^-1]. The lists are fixed; a
// term exists in a build when its process is configured (heat_active_ /
// cool_active_), and the query API reports the active ones in this order.
enum HeatTerm {
  H_CUSTOM, H_CR, H_CHEM, H_DH, H_COMPRESS, H_VISC, H_PHOTODUST,
  H_NETPHOTODUST, H_PHOTOAV, H_PHOTO, H_ZCIE, H_NUCLEAR, NHEAT_ALL
};
// the coolant-file coolants (thermal.coolant_files) take the slots after the fixed
// terms, one each: C_FILE0 + i, named by the generated network
#ifndef KROME_COOLFILE_NAMES
#define KROME_COOLFILE_NAMES
#define KROME_COOLFILE_ACTIVE
constexpr int NCOOL_FILE = 0;
template <class P>
KROME_HD inline void cooling_files_(double, const double*, const P&,
                                    double*) {}
#endif
enum CoolTerm {
  C_CUSTOM, C_ATOMIC, C_H2, C_HD, C_MOL, C_Z, C_ZCIE, C_ZEXT,
  C_FLOOR, C_CHEM, C_COMPTON, C_CIE, C_EXPANSION, C_CONT, C_FF, C_DH,
  C_DUST, C_DUSTTAB, C_DUSTBINS, C_FILE0,
  NCOOL_ALL = C_FILE0 + NCOOL_FILE
};
inline const char* heat_term_names[NHEAT_ALL] = {
    "CUSTOM", "CR", "CHEM", "DH", "COMPRESS", "VISCOUS", "PHOTODUST",
    "PHOTODUSTNET", "PHOTOAV", "PHOTO", "ZCIE", "NUCLEAR"};
inline const char* cool_term_names[NCOOL_ALL] = {
    "CUSTOM", "ATOMIC", "H2", "HD", "MOL", "Z", "Z_CIE",
    "Z_EXTENDED", "FLOOR", "CHEM", "COMPTON", "CIE", "EXPANSION", "CONT",
    "FF", "DH", "DUST", "DUSTTAB", "DUSTBINS" KROME_COOLFILE_NAMES};
// PHOTO: photoheating, or the photoelectron heat of radiation.secondary_ionization
#ifdef KROME_HAS_SECION
constexpr bool heat_photo_active_ = true;
#else
constexpr bool heat_photo_active_ = TC::heat_photo;
#endif
constexpr bool heat_active_[NHEAT_ALL] = {
    HAS_CUSTOM_HEATING, TC::heat_CR, TC::heat_chem, TC::heat_dH,
    TC::heat_compress, TC::heat_visc, TC::heat_photoDust,
    TC::heat_netPhotoDust, TC::heat_photoAv, heat_photo_active_,
    TC::heat_ZCIE, TC::heat_nuclear};
constexpr bool cool_active_[NCOOL_ALL] = {
    HAS_CUSTOM_COOLING, TC::cool_atomic, TC::cool_H2, TC::cool_HD,
    TC::cool_mol, TC::cool_Z, TC::cool_ZCIE, TC::cool_Zext,
    TC::cooling_floor, TC::cool_chem, TC::cool_compton, TC::cool_CIE,
    TC::cool_expansion,
    TC::cool_cont, TC::cool_ff, TC::cool_dH, TC::cool_dust,
    TC::cool_dustTab, TC::cool_dust_bins KROME_COOLFILE_ACTIVE};
constexpr int count_active_(const bool* a, int n) {
  int m = 0;
  for (int i = 0; i < n; ++i) m += a[i] ? 1 : 0;
  return m;
}
constexpr int NHEAT = count_active_(heat_active_, NHEAT_ALL);
constexpr int NCOOL = count_active_(cool_active_, NCOOL_ALL);

// the cooling floor (thermal.floor): the floored terms evaluated at
// Tfloor, subtracted from the cooling at Tgas, i.e. each floored cooler
// becomes cool(Tgas) - cool(Tfloor) (legacy individual floors). For
// Z_EXTENDED legacy floors only the metal-line part of the tanh merge,
// with the merge weight at Tgas: f2(Tgas) cooling_Z(Tfloor).
KROME_HD inline double floor_cool_(double Tgas, double Tf, const double* n,
                          ThermCtx* c) {
  double f = 0.0;
#ifdef KROME_COOLING_HPP_INCLUDED
  if constexpr (TC::floor_atomic)
    f += cooling_Atomic(Tf, n) + cooling_recomb_(Tf, n, c->p);
  if constexpr (TC::floor_H2) f += cooling_H2(Tf, n) * h2_opacity_(Tf, n, c->p);
  if constexpr (TC::floor_HD) f += cooling_HD(Tf, n, c->p);
#endif
  if constexpr (TC::floor_mol) f += cooling_mol_floor(Tf, n, c->p);
  if constexpr (TC::floor_Z) f += cooling_Z(Tf, n, c->p);
  if constexpr (TC::floor_file) {
    double fc[NCOOL_FILE > 0 ? NCOOL_FILE : 1];
    cooling_files_(Tf, n, c->p, fc);
    for (int i = 0; i < NCOOL_FILE; ++i) f += fc[i];
  }
  if constexpr (TC::floor_ZCIE) f += cooling_ZCIE(Tf, n, c->p);
  if constexpr (TC::floor_Zext) {
    const double f2 = (std::tanh(1e-3 * (-Tgas + 1e4)) + 1.0) * 0.5;
    f += f2 * cooling_Z(Tf, n, c->p);
  }
  (void)Tgas; (void)Tf; (void)n; (void)c;
  return f;
}

// every thermal term into H[NHEAT_ALL] / C[NCOOL_ALL] (inactive ones 0);
// the secondary-ionization fluxes into dnx[NSPEC] (zeroed here). dustTd:
// per-bin dust temperatures when they are ODE state (legacy -usedTdust
// variant); nullptr = solve them by bisection inside DustBins::cooling.
// FLOOR is the (negative) cooling-floor correction; CHEM heating is the
// H2 formation heating, CHEM cooling the collisional dissociation.
KROME_HD inline void thermal_terms_(double Tgas, const double* n, const double* k,
                           ThermCtx* c, double* dnx, const double* dustTd,
                           double* H, double* C) {
  (void)dustTd;
  (void)k;
  for (int i = 0; i < NSPEC; ++i) dnx[i] = 0.0;
  for (int i = 0; i < NHEAT_ALL; ++i) H[i] = 0.0;
  for (int i = 0; i < NCOOL_ALL; ++i) C[i] = 0.0;
  H[H_CUSTOM] = heating_custom(Tgas, n, c->p);
  if constexpr (TC::heat_CR) H[H_CR] = heating_CR(n, k, c->p);
  if constexpr (TC::heat_chem) H[H_CHEM] = heating_chem(Tgas, n, k, c->p);
  if constexpr (TC::heat_dH) H[H_DH] = heating_dH(Tgas, n, k, c->p);
  if constexpr (TC::heat_compress)
    H[H_COMPRESS] = heat_compress_(Tgas, n, c->p);
  if constexpr (TC::heat_visc) H[H_VISC] = heat_visc_(n, c->p);
  if constexpr (TC::heat_nuclear)
    H[H_NUCLEAR] = heating_nuclear(Tgas, n, k, c->p);
  if constexpr (TC::heat_photoDust || TC::heat_netPhotoDust) {
    // Habing flux: the photobin field's 6-13.6 eV integral when radiation
    // is present, else the G0 param attenuated by dust (legacy exp(-2.5 Av))
    double Ghab = habing_param_(c->p);
#ifdef KROME_PHOTO_HPP_INCLUDED
    if (c->photo_on)
      Ghab = photo_habing(c->ph, c->Phi, c->psum)
          * habing_shield_(Tgas, n, c->p);
#endif
    if constexpr (TC::heat_photoDust)
      H[H_PHOTODUST] = heating_photoDust(Tgas, n, c->p, Ghab);
    if constexpr (TC::heat_netPhotoDust)
      H[H_NETPHOTODUST] = heating_netPhotoDust(Tgas, n, c->p, Ghab);
  }
  if constexpr (TC::heat_photoAv)
    H[H_PHOTOAV] = heating_photoAv(Tgas, n, k, c->p);
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (c->photo_on) {
#ifdef KROME_HAS_SECION
    // radiation.secondary_ionization: the photoelectron deposit is split into
    // secondary ionizations (species fluxes) and heat (Quadri+25)
    H[H_PHOTO] = photo_secondary(c->ph, c->Phi, n, dnx, c->shield, c->psum);
#else
    if constexpr (TC::heat_photo)
      H[H_PHOTO] = photo_heating(c->ph, c->Phi, n, c->shield, c->psum);
#endif
  }
#endif
  // Shen+13 ZCIE heating, always with the ZCIE cooling; under Z_EXTENDED
  // it gets the ZCIE cooling's merge weight (the network tracks the
  // metals below 1e4 K: unweighted it would double count their heating)
  if constexpr (TC::heat_ZCIE) {
    double hz = heating_ZCIE(Tgas, n, c->p);
    if constexpr (TC::cool_Zext)
      hz *= (std::tanh(1e-3 * (Tgas - 1e4)) + 1.0) * 0.5;
    H[H_ZCIE] = hz;
  }

  C[C_CUSTOM] = cooling_custom(Tgas, n, c->p);
  // radiative line coolers; with a cooling floor the floored ones are
  // cut at Tfloor (C_FLOOR), so the gas cannot radiate itself below it
#ifdef KROME_COOLING_HPP_INCLUDED
  if constexpr (TC::cool_atomic)
    C[C_ATOMIC] = cooling_Atomic(Tgas, n) + cooling_recomb_(Tgas, n, c->p);
  if constexpr (TC::cool_H2)
    C[C_H2] = cooling_H2(Tgas, n) * h2_opacity_(Tgas, n, c->p);
  if constexpr (TC::cool_HD) C[C_HD] = cooling_HD(Tgas, n, c->p);
#endif
  if constexpr (TC::cool_mol) C[C_MOL] = cooling_mol(Tgas, n, c->p);
  if constexpr (TC::cool_Z) C[C_Z] = cooling_Z(Tgas, n, c->p);
  if constexpr (NCOOL_FILE > 0) cooling_files_(Tgas, n, c->p, C + C_FILE0);
  if constexpr (TC::cool_ZCIE) C[C_ZCIE] = cooling_ZCIE(Tgas, n, c->p);
  if constexpr (TC::cool_Zext) {  // legacy cooling_ZExtended tanh merge
    const double f1 = (std::tanh(1e-3 * (Tgas - 1e4)) + 1.0) * 0.5;
    const double f2 = (std::tanh(1e-3 * (-Tgas + 1e4)) + 1.0) * 0.5;
    C[C_ZEXT] = f1 * cooling_ZCIE(Tgas, n, c->p)
        + f2 * cooling_Z(Tgas, n, c->p);
  }
  if constexpr (TC::cooling_floor)
    C[C_FLOOR] = -floor_cool_(Tgas, cooling_floor_T_(c->p), n, c);
#ifdef KROME_COOLING_HPP_INCLUDED
  if constexpr (TC::cool_compton) C[C_COMPTON] = cooling_Compton(Tgas, n, c->p);
  if constexpr (TC::cool_CIE) C[C_CIE] = cooling_CIE(Tgas, n, c->p);
  if constexpr (TC::cool_expansion)
    C[C_EXPANSION] = cooling_expansion(Tgas, n, c->p);
  if constexpr (TC::cool_cont) C[C_CONT] = cooling_Continuum(Tgas, n);
#endif
  if constexpr (TC::cool_chem) C[C_CHEM] = cooling_chem(Tgas, n, k, c->p);
  if constexpr (TC::cool_ff) C[C_FF] = cooling_ff(Tgas, n);
  if constexpr (TC::cool_dH) C[C_DH] = cooling_dH(Tgas, n, k, c->p);
#ifdef KROME_DUST_RATES_HPP_INCLUDED
  if constexpr (TC::cool_dust)  // the Tdust parameter without tables
    C[C_DUST] = cooling_dustExchange(Tgas, table_Tdust(Tgas, n, c->p), n,
                                     c->p);
#endif
  if constexpr (TC::cool_dustTab) C[C_DUSTTAB] = cooling_dustTab(Tgas, n, c->p);
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if constexpr (TC::cool_dust_bins) {
    const double* phi = nullptr;
#ifdef KROME_PHOTO_HPP_INCLUDED
    phi = c->Phi;
#endif
    // dust.evolve_Td: the gas-grain exchange at the state Td (energy
    // conserving); else the balance Td of this evaluation (radiated =
    // exchanged at the balance)
    C[C_DUSTBINS] = dustTd
        ? dust_cooling_exchange(c->dust, Tgas, n, c->p, dustTd)
        : dust_cooling_state(c->dust, Tgas, n, c->p,
                             c->dust_tdb ? c->dust_tdb : c->dust_Td, phi);
  }
#endif
}

// the adiabatic index the thermal ODE uses
KROME_HD inline double gamma_ad_(double Tgas, const double* n, const Params& p) {
  return gamma_index(Tgas, n, p);
}

// everything beyond pure chemistry: dTgas/dt from the thermal terms, and
// the secondary-ionization fluxes in dnx[NSPEC] (see thermal_terms_)
KROME_HD inline double thermal_aux(double Tgas, const double* n, const double* k,
                          ThermCtx* c, double* dnx,
                          const double* dustTd = nullptr) {
  double H[NHEAT_ALL], C[NCOOL_ALL];
  thermal_terms_(Tgas, n, k, c, dnx, dustTd, H, C);
  double heat = 0.0, cool = 0.0;
  for (int i = 0; i < NHEAT_ALL; ++i) heat += H[i];
  for (int i = 0; i < NCOOL_ALL; ++i) cool += C[i];
  return (gamma_ad_(Tgas, n, c->p) - 1.0) * (heat - cool)
      / (KB_ERG * ntot_(n));
}


// ---- one cell's rhs and Jacobian: the host's and the device's -----------
//
// RhsCtx: the state layout of one integration and the work arrays of its
// rhs and Jacobian, as plain pointers (host: SolveCtx's vectors; device:
// a cell's scratch rows or thread-local arrays), so that the functions
// below (KROME_HD) are the one implementation of both backends.
struct RhsCtx {
  ThermCtx* th = nullptr;  // parameters, rates, shields, photo, dust
  double Tfix = 0.0;  // Tgas when it is not evolved
  int NEQ = NSPEC;
  int iT = -1;        // Tgas slot
  int iu = -1, B = 0; // first ln Phi slot, number of photo bins evolved
  int id = -1, D = 0; // first Tdust slot, number of dust bins evolved
  int ia = -1, A = 0; // first grain-atoms slot, number of bins (dust.evolve)
  double xds = 0.0;   // the grain count scale of the solve
  bool hit_T0 = false;      // the rhs refused a trial state at Tgas <= 0
  const double* kfix = nullptr;  // fixed Tgas: rates_fixed of the solve
  // evolved Tgas: rates_fixed at the rate temperature kT_T, reused while
  // it stays the same (the species columns of a difference-quotient
  // Jacobian); nullptr: off (binned dust, whose Td-dependent rates change
  // with the state at the same Tgas, or the device)
  double* kT = nullptr;
  double kT_T = 0.0;
  bool kT_ok = false;
  double* Phi = nullptr;    // [B] the evolved field of the state
  double* kabs = nullptr;   // [B] photon removal coefficients
  double* arad = nullptr;   // [A] the radii of the state (from N_b)
  double* dnd = nullptr;    // [NSPEC] the key elements' exchange with grains
  double* tdb = nullptr;    // [nb] balance Td of the binned dust (without
                            // dust.evolve_Td; nullptr: none)
  double* awcur = nullptr;  // [nb] area weights of the state's radii (the
                            // dust hook's)
  double* U = nullptr;      // [NSPEC * NAGG] the split Jacobian's coupling
  double *k0 = nullptr, *kp = nullptr, *kunit = nullptr;  // [NREA + 1]
  double *dnf = nullptr, *dnx = nullptr, *dnx0 = nullptr,
         *dnx1 = nullptr;                                 // [NSPEC]
  double *td0 = nullptr, *td1 = nullptr;                  // [D + 1]
  double *da0 = nullptr, *da1 = nullptr;                  // [A]
  double *dn0 = nullptr, *dn1 = nullptr;                  // [NSPEC]
  double* fm = nullptr;                                   // [NEQ]
};

// the doubles of the RhsCtx work arrays (rhs_bind_ lays them out); nb: the
// dust bins (0: none)
KROME_HD constexpr long rhs_work_doubles(int NEQ, int B, int nb, int D,
                                         int A) {
  return 2L * B + 3L * A + 2L * nb + static_cast<long>(NSPEC) * NAGG
         + 3L * (NREA + 1) + 7L * NSPEC + 2L * (D + 1) + NEQ;
}
// point the work arrays of c into w[rhs_work_doubles(...)] (nb: the dust
// bins; tdb only when the balance Td is used); returns the end of them
KROME_HD inline double* rhs_bind_(RhsCtx* c, double* w, int nb,
                                  bool balance) {
  c->Phi = w;
  w += c->B;
  c->kabs = w;
  w += c->B;
  c->arad = w;
  w += c->A;
  c->da0 = w;
  w += c->A;
  c->da1 = w;
  w += c->A;
  c->awcur = w;
  w += nb;
  c->tdb = balance ? w : nullptr;
  w += nb;
  c->U = w;
  w += static_cast<long>(NSPEC) * NAGG;
  double** k3[3] = {&c->k0, &c->kp, &c->kunit};
  for (double** q : k3) {
    *q = w;
    w += NREA + 1;
  }
  double** n7[7] = {&c->dnd, &c->dnf, &c->dnx, &c->dnx0, &c->dnx1, &c->dn0,
                    &c->dn1};
  for (double** q : n7) {
    *q = w;
    w += NSPEC;
  }
  c->td0 = w;
  w += c->D + 1;
  c->td1 = w;
  w += c->D + 1;
  c->fm = w;
  return w + c->NEQ;
}

// the temperature the rates see: the state (or the fixed input) capped
// to [Tcmb, 1e9 K], as legacy coe did
KROME_HD inline double state_T_(const double* y, const RhsCtx* c) {
  const double T = c->iT >= 0 ? y[c->iT] : c->Tfix;
  return std::fmin(std::fmax(T, c->th->p.Tcmb), 1e9);
}

// the temperature the thermal terms (heating, cooling, gamma, the dust
// energy exchange) see: capped to [TTHERMAL_MIN, 1e9 K] with the fixed
// legacy floor (fex), NOT Tcmb. A cap at Tcmb would freeze every cooling
// at cool(Tcmb) below Tcmb, so dT/dt stays negative and Tgas runs through
// zero (collapseZ at Tcmb = 43.7 K); this way the cooling keeps falling
// with Tgas down to 2.73 K
constexpr double TTHERMAL_MIN = 2.73;
KROME_HD inline double thermal_T_(double T) {
  return std::fmin(std::fmax(T, TTHERMAL_MIN), 1e9);
}
KROME_HD inline double thermal_T_(const double* y, const RhsCtx* c) {
  return thermal_T_(c->iT >= 0 ? y[c->iT] : c->Tfix);
}

KROME_HD inline void unpack_rad_(const double* y, RhsCtx* c) {
  for (int b = 0; b < c->B; ++b) c->Phi[b] = std::exp(y[c->iu + b]);
}

// rate coefficients at state y: generated rates + photo rates of the
// field, with their shielding factors (kept in th->shield for the
// photoheating, secondary ionization and photon removal at this state)
KROME_HD inline void eval_rates_(double T, const double* y, RhsCtx* c,
                                 double* k) {
  ThermCtx& tc = *c->th;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  // dust.evolve: the grain radii of this state, before anything that
  // uses the cross sections (the balance, the rates' area weights);
  // tc.dust.acur points at them
  if (c->ia >= 0) {
    for (int b = 0; b < c->A; ++b)
      c->arad[b] = dust_radius_of_atoms(tc.dust, b, y[c->ia + b], c->xds);
    dust_area_weights(tc.dust, c->awcur);
  }
  // binned dust without dust.evolve_Td: the per-bin Td are the radiative-
  // collisional balance at this state, computed before the rates (which
  // see them through the hook) and used by the DUSTBINS cooling
#ifdef KROME_HAS_PAR_d2g
  if (c->tdb != nullptr) {
    const double* phi = nullptr;
#ifdef KROME_PHOTO_HPP_INCLUDED
    phi = tc.Phi;
#endif
    dust_balance(tc.dust, T, y, tc.p, phi, c->tdb);
    tc.p.hk.dust.Td = c->tdb;
    tc.dust_tdb = c->tdb;
  }
#endif
#endif
  if (c->kfix != nullptr) {  // = rates(): zeros + rates_fixed, then rates_var
    const double* kf = c->kfix;
    for (int i = 0; i < NREA; ++i) k[i] = kf[i];
    rates_var(T, y, tc.p, k);
  } else if (c->kT != nullptr) {  // the same, rates_fixed cached per Tgas
    if (!c->kT_ok || T != c->kT_T) {
      for (int i = 0; i < NREA; ++i) c->kT[i] = 0.0;
      rates_fixed(T, y, tc.p, c->kT);
      c->kT_T = T;
      c->kT_ok = true;
    }
    for (int i = 0; i < NREA; ++i) k[i] = c->kT[i];
    rates_var(T, y, tc.p, k);
  } else {
    rates(T, y, tc.p, k);
  }
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (tc.photo_on) {
    if constexpr (PHOTO_SHIELDED) photo_shields(T, y, tc.p, tc.shield);
    photo_rates(tc.ph, tc.Phi, k, tc.shield, tc.psum);
  }
#endif
}

#ifdef KROME_PHOTO_HPP_INCLUDED
// photon removal coefficient kabs[B] of the photochemistry at state y
// (th->shield must hold the shields at y: eval_rates_ at y)
KROME_HD inline void photo_absorption_(double T, const double* y, RhsCtx* c,
                                       double* kabs) {
  ThermCtx& tc = *c->th;
  if constexpr (NLW > 0) lw_unit_rates(T, y, tc.p, tc.khat);
  photo_absorption(tc.ph, y, tc.shield, NLW > 0 ? tc.khat : nullptr, kabs);
}
#endif

// the thermal rows at state y (rates k at the rate temperature): returns
// dTgas/dt, fills the secondary-ionization fluxes dnx[NSPEC] and, with
// dust.evolve_Td, dTdust/dt into dtd[D]; every thermal term at thermal_T_
KROME_HD inline double aux_(const double* y, const double* k, RhsCtx* c,
                            double* dnx, double* dtd) {
  const double T = thermal_T_(y, c);
  const double* Tds = c->id >= 0 ? y + c->id : nullptr;
  const double Tdot = thermal_aux(T, y, k, c->th, dnx, Tds);
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if constexpr (EvolveConfig::Td)
    if (c->id >= 0) {
      const double* phi = nullptr;
#ifdef KROME_PHOTO_HPP_INCLUDED
      phi = c->th->Phi;
#endif
      dust_dTdust(c->th->dust, T, y, c->th->p, Tds, phi, dtd);
    }
#else
  (void)dtd;
#endif
  return Tdot;
}

#ifdef KROME_DUST_BINS_HPP_INCLUDED
// dust.evolve: dN/dt of the grain slots into dad and the key-element
// exchange into dnd[NSPEC] (zeroed here), at state y (the rates, hence the
// radii and the balance Td, evaluated at y by eval_rates_)
KROME_HD inline void dust_rhs_(const double* y, RhsCtx* c, double* dad,
                               double* dnd) {
  for (int i = 0; i < NSPEC; ++i) dnd[i] = 0.0;
  const DustView& d = c->th->dust;
#ifdef KROME_HAS_PAR_d2g  // every dust.evolve network has d2g
  const double* Tds = c->id >= 0 ? y + c->id
      : (c->tdb != nullptr ? c->tdb : c->th->dust_Td);
  dust_dadt(d, state_T_(y, c), y, c->th->p, Tds, dad, dnd);
#else
  (void)y; (void)dad;
#endif
  // the state is N_b: dN/dt = 3 N/a da/dt = n_d 4 pi a^2 rho0/m da/dt,
  // exactly what dadt took from / gave to the key element
  for (int b = 0; b < c->A; ++b) {
    const double ab = c->arad[b];
    dad[b] *= 3.0 * dust_atoms_per_a3(d, b, c->xds) * ab * ab;
  }
}
#endif

// network freeze directives: zero the species rows of f (length NSPEC
// prefix of the state) while their condition holds
KROME_HD inline void freeze_rows_(const RhsCtx* c, double* f) {
  if constexpr (NFREEZE_MAX > 0) {
    int idx[NFREEZE_MAX > 0 ? NFREEZE_MAX : 1];
    const int m = frozen_species(c->th->p, idx);
    for (int q = 0; q < m; ++q) f[idx[q]] = 0.0;
  } else {
    (void)c; (void)f;
  }
}

// one-sided difference step for density column j of the aux rows (thermal,
// dust, secondary ionization): relative to n_j, but at least 1e-9 of the
// total density. The aux terms depend on n_j also through sums (ntot,
// Hnuclei, x_e); a step relative to a trace species alone would move those
// sums below roundoff and yield an exactly zero derivative. Terms linear
// in n_j are exact for any step.
KROME_HD inline double fd_step_n_(const double* y, int j) {
  double s = 0.0;
  for (int i = 0; i < NSPEC; ++i) s += std::fabs(y[i]);
  return std::fmax(std::fmax(1e-7 * std::fabs(y[j]), 1e-9 * s), 1e-30);
}

#if defined(KROME_PHOTO_HPP_INCLUDED) && defined(KROME_HAS_SECION)
// secondary ionization by fast photoelectrons is chemistry, active at any
// Tgas: with Tgas evolved its species fluxes come from thermal_aux (with
// the photoelectron heat); at fixed Tgas they are computed here alone
KROME_HD inline void secion_fixedT_(const double* n, RhsCtx* c, double* dnx) {
  ThermCtx& tc = *c->th;
  for (int i = 0; i < NSPEC; ++i) dnx[i] = 0.0;
  if (tc.photo_on) {
    if constexpr (PHOTO_SHIELDED)
      photo_shields(state_T_(n, c), n, tc.p, tc.shield);
    photo_secondary(tc.ph, tc.Phi, n, dnx, tc.shield, tc.psum);
  }
}
#endif

KROME_HD inline void full_rhs_(const double* y, double* f, RhsCtx* c) {
  ThermCtx& tc = *c->th;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  // dust-chemistry hook: with Tdust as ODE state, the rates see the
  // current state slice (else eval_rates_ installs the balance Td of the
  // state; with dust.evolve_Td and Tgas fixed, the caller's fixed Td)
  if (c->id >= 0) tc.p.hk.dust.Td = y + c->id;
#endif
  if (c->iu >= 0) unpack_rad_(y, c);
  const double T = state_T_(y, c);
  eval_rates_(T, y, c, tc.k);
  rhs(y, tc.k, f);
  if (c->iT >= 0) {
    double* dnx = c->dnx;
    f[c->iT] = aux_(y, tc.k, c, dnx, c->id >= 0 ? f + c->id : nullptr);
    for (int i = 0; i < NSPEC; ++i) f[i] += dnx[i];
  }
#if defined(KROME_PHOTO_HPP_INCLUDED) && defined(KROME_HAS_SECION)
  else {
    double* dnx = c->dnx;
    secion_fixedT_(y, c, dnx);
    for (int i = 0; i < NSPEC; ++i) f[i] += dnx[i];
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (c->ia >= 0) {  // grain atoms and the key-element exchange
    dust_rhs_(y, c, f + c->ia, c->dnd);
    for (int i = 0; i < NSPEC; ++i) f[i] += c->dnd[i];
  }
#endif
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (c->iu >= 0) {
    if constexpr (PHOTO_SHIELDED)  // secion_fixedT_ may have moved them
      photo_shields(T, y, tc.p, tc.shield);
    photo_absorption_(T, y, c, c->kabs);
    for (int b = 0; b < c->B; ++b) f[c->iu + b] = -CLIGHT * c->kabs[b];
  }
#endif
  freeze_rows_(c, f);
}

// a trial state the rhs can evaluate: a non-positive gas (or dust)
// temperature is not physical (the rhs reports it as recoverable, so the
// integrator retries with a smaller step)
KROME_HD inline bool rhs_trial_ok_(const double* y, const RhsCtx* c) {
  if (c->iT >= 0 && !(y[c->iT] > 0.0)) return false;
  for (int b = 0; b < c->D; ++b)
    if (!(y[c->id + b] > 0.0)) return false;
  return true;
}

// the integrator's view of the rhs on the host: a non-finite rhs is
// reported as a recoverable failure, so the integrator retries with a
// smaller step instead of accepting it
inline int rhs_host_(const double* y, double* f, RhsCtx* c) {
  if (!rhs_trial_ok_(y, c)) {
    c->hit_T0 = true;
    return 1;
  }
  full_rhs_(y, f, c);
  for (int i = 0; i < c->NEQ; ++i)
    if (!std::isfinite(f[i])) return 1;
  return 0;
}

// a dense NEQ x NEQ matrix, element (i, j) at d[i rs + j cs]: the host
// integrator's Newton matrix is column-major (rs 1, cs n); the device's
// has the cell's stride (interleaved: rs S, cs n S)
struct DenseMat_ {
  double* d;
  long rs, cs;
  KROME_HD DenseMat_(double* d_, int n) : d(d_), rs(1), cs(n) {}
  KROME_HD DenseMat_(double* d_, long rs_, long cs_)
      : d(d_), rs(rs_), cs(cs_) {}
};
#define KROME_JEL(J, i, j) \
  ((J).d[static_cast<long>(i) * (J).rs + static_cast<long>(j) * (J).cs])

#ifdef KROME_DUST_BINS_HPP_INCLUDED
// dust.evolve Jacobian entries: the grain columns by central differences
// of the full rhs, and the dependence of the grain rows and of the
// key-element exchange on n_j (difference quotients, Td held at the
// state's); yp, fp: NEQ scratch
KROME_HD inline void jac_dust_(const double* y, DenseMat_ J, RhsCtx* c,
                               double* yp, double* fp) {
  const int NEQ = c->NEQ;
  double* fm = c->fm;
  for (int i = 0; i < NEQ; ++i) yp[i] = y[i];
  for (int b = 0; b < c->A; ++b) {
    const int col = c->ia + b;
    const double h = std::fmax(1e-6 * std::fabs(y[col]), 1e-20);
    yp[col] = y[col] + h;
    full_rhs_(yp, fp, c);
    yp[col] = y[col] - h;
    full_rhs_(yp, fm, c);
    yp[col] = y[col];
    for (int i = 0; i < NEQ; ++i)
      KROME_JEL(J, i, col) = (fp[i] - fm[i]) / (2.0 * h);
  }
  eval_rates_(state_T_(y, c), y, c, c->kp);  // radii, Td at y
  double *da0 = c->da0, *da1 = c->da1, *dn0 = c->dn0, *dn1 = c->dn1;
  dust_rhs_(y, c, da0, dn0);
  for (int j = 0; j < NSPEC; ++j) {
    const double h = fd_step_n_(y, j);
    yp[j] = y[j] + h;
    dust_rhs_(yp, c, da1, dn1);
    yp[j] = y[j];
    for (int b = 0; b < c->A; ++b)
      KROME_JEL(J, c->ia + b, j) += (da1[b] - da0[b]) / h;
    for (int i = 0; i < NSPEC; ++i)
      if (dn0[i] != 0.0 || dn1[i] != 0.0)
        KROME_JEL(J, i, j) += (dn1[i] - dn0[i]) / h;
  }
}
#endif

// a species' tolerance: its SolverOpts::atol_species / rtol_species
// entry, or the scalar when there is no array or the entry is negative
inline double species_tol_(const double* v, int i, double scalar) {
  return v != nullptr && !(v[i] < 0.0) ? v[i] : scalar;
}

// a per-species rtol that differs from the scalar somewhere: only then do
// the integrator use per-slot error weights (otherwise CVODE's SV form)
inline bool rtol_varies_(const SolverOpts& o) {
  if (o.rtol_species == nullptr) return false;
  for (int i = 0; i < NSPEC; ++i)
    if (species_tol_(o.rtol_species, i, o.rtol) != o.rtol) return true;
  return false;
}

// JAC_ANALYTIC assembly (see the header comment) into J, zeroed by the
// caller; yp, fp: NEQ scratch
KROME_HD inline int jac_assemble_(const double* y, DenseMat_ J, RhsCtx* c,
                                  double* yp, double* fp) {
  ThermCtx& tc = *c->th;
  const int NEQ = c->NEQ;
  if (c->iu >= 0) unpack_rad_(y, c);
  const double T = state_T_(y, c);
  // reference rates at y, kept apart from tc.k (the difference quotients
  // below re-evaluate the rhs, which overwrites tc.k)
  double* k0 = c->k0;
  eval_rates_(T, y, c, k0);

  // chemistry block: analytic (the split Jacobian J_direct + U W^T, as
  // the generated jacobian())
  double* U = c->U;
  jacobian_split(T, y, tc.p, k0, J.d, U, J.rs, J.cs);
  for (int i = 0; i < NSPEC; ++i)
    for (int a = 0; a < NAGG; ++a) {
      const double u = U[i * NAGG + a];
      if (u == 0.0) continue;
      for (int j = 0; j < NSPEC; ++j) KROME_JEL(J, i, j) += u * agg_w_(a, j);
    }

#ifdef KROME_PHOTO_HPP_INCLUDED
  if (tc.photo_on && (PHOTO_SHIELDED || (NLW > 0 && c->iu >= 0))) {
    // n-dependence of the shielding factors and of the J_LW unit rates
    // (difference quotients): d k_r/d n_j = d s_r/d n_j sum_b sigma_rb Phi_b
    // enters the chemistry rows; with a co-evolved field both also enter
    // the removal rows, d(du_b)/d n_j = -c sigma_rb n_r ds_r/dn_j / pdiss
    // (photo) and -c lww_b n_q dkhat_q/dn_j / pdiss (J_LW; the J_LW rates'
    // own dk/dn is in the generated chemistry block)
    const PhotoView& ph = tc.ph;
    double s0[NPHOTO + 1], s1[NPHOTO + 1], kh0[NLW + 1], kh1[NLW + 1];
    double u0[NPHOTO + 1];
    for (int r = 0; r < NPHOTO; ++r) {
      s0[r] = tc.shield[r];
      u0[r] = 0.0;
      for (int b = 0; b < ph.B; ++b)
        u0[r] += ph.sigma[static_cast<size_t>(r) * ph.B + b] * tc.Phi[b];
    }
    if constexpr (NLW > 0) lw_unit_rates(T, y, tc.p, kh0);
    double* kunit = c->kunit;
    double* dnf = c->dnf;
    double ys[NSPEC];
    for (int i = 0; i < NSPEC; ++i) ys[i] = y[i];
    for (int j = 0; j < NSPEC; ++j) {
      const double h = fd_step_n_(y, j);
      ys[j] = y[j] + h;
      if constexpr (PHOTO_SHIELDED) photo_shields(T, ys, tc.p, s1);
      if constexpr (NLW > 0)
        if (c->iu >= 0) lw_unit_rates(T, ys, tc.p, kh1);
      ys[j] = y[j];
      for (int r = 0; r < NPHOTO; ++r) {
        if (!photo_shielded[r]) continue;
        const double ds = (s1[r] - s0[r]) / h;
        if (ds == 0.0) continue;
        for (int i = 0; i < NREA; ++i) kunit[i] = 0.0;
        kunit[photo_reaction[r]] = 1.0;
        rhs(y, kunit, dnf);
        for (int i = 0; i < NSPEC; ++i)
          if (dnf[i] != 0.0) KROME_JEL(J, i, j) += ds * u0[r] * dnf[i];
        if (c->iu >= 0) {
          const double f = -CLIGHT * y[photo_reactant[r]] * ds
              / photo_pdiss[r];
          for (int b = 0; b < c->B; ++b)
            KROME_JEL(J, c->iu + b, j) +=
                f * ph.sigma[static_cast<size_t>(r) * c->B + b];
        }
      }
      if constexpr (NLW > 0)
        if (c->iu >= 0)
          for (int q = 0; q < NLW; ++q) {
            const double dk = (kh1[q] - kh0[q]) / h;
            if (dk == 0.0) continue;
            const double f = -CLIGHT * y[lw_reactant[q]] * dk / lw_pdiss[q];
            for (int b = 0; b < c->B; ++b)
              KROME_JEL(J, c->iu + b, j) += f * ph.lww[b];
          }
    }
  }
  if (c->iu >= 0) {
    const PhotoView& ph = tc.ph;
    // d(du_b)/d(n_j) = -c s_r sigma_rb / pdiss_r for each photo reaction
    // with reactant j, -c khat_q lww_b / pdiss_q for each J_LW rate
    // (constant in Phi; d(du_b)/d(u) = 0 for pure absorption)
    for (int r = 0; r < NPHOTO; ++r)
      for (int b = 0; b < c->B; ++b) {
        const double srb = ph.sigma[static_cast<size_t>(r) * c->B + b];
        if (srb != 0.0)
          KROME_JEL(J, c->iu + b, photo_reactant[r]) +=
              -CLIGHT * srb * tc.shield[r] / photo_pdiss[r];
      }
    double kh[NLW + 1];
    if constexpr (NLW > 0) {
      lw_unit_rates(T, y, tc.p, kh);
      for (int q = 0; q < NLW; ++q)
        for (int b = 0; b < c->B; ++b)
          if (ph.lww[b] != 0.0)
            KROME_JEL(J, c->iu + b, lw_reactant[q]) +=
                -CLIGHT * kh[q] * ph.lww[b] / lw_pdiss[q];
    }
    // d(dn_i)/d(u_b) = S_ir s_r sigma_rb Phi_b n_reactant (photo) and
    // S_iq khat_q lww_b Phi_b n_reactant (J_LW): the stoichiometric
    // column S_ir n_reactant comes from rhs() with a unit k on reaction r
    // (rhs is linear in k), so no product lists are needed
    double* kunit = c->kunit;
    double* dnf = c->dnf;
    for (int r = 0; r < NPHOTO + NLW; ++r) {
      const bool isph = r < NPHOTO;
      const int q = r - NPHOTO;
      for (int i = 0; i < NREA; ++i) kunit[i] = 0.0;
      kunit[isph ? photo_reaction[r] : lw_reaction[q]] = 1.0;
      rhs(y, kunit, dnf);
      for (int b = 0; b < c->B; ++b) {
        const double srb = isph
            ? ph.sigma[static_cast<size_t>(r) * c->B + b] * tc.shield[r]
            : kh[q] * ph.lww[b];
        if (srb == 0.0) continue;
        const double fb = srb * c->Phi[b];
        for (int i = 0; i < NSPEC; ++i)
          if (dnf[i] != 0.0) KROME_JEL(J, i, c->iu + b) += fb * dnf[i];
      }
    }
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (c->ia >= 0) jac_dust_(y, J, c, yp, fp);
#endif
  if (c->iT < 0) {  // fixed Tgas: no thermal rows
#if defined(KROME_PHOTO_HPP_INCLUDED) && defined(KROME_HAS_SECION)
    // the secondary-ionization rows: difference quotients vs n_j and,
    // with a co-evolved field, vs u_b (as the evolved-Tgas aux rows)
    if (tc.photo_on) {
      double* dnx0 = c->dnx0;
      double* dnx1 = c->dnx1;
      for (int i = 0; i < NEQ; ++i) yp[i] = y[i];
      secion_fixedT_(y, c, dnx0);
      for (int j = 0; j < NSPEC; ++j) {
        const double h = fd_step_n_(y, j);
        yp[j] = y[j] + h;
        secion_fixedT_(yp, c, dnx1);
        yp[j] = y[j];
        for (int i = 0; i < NSPEC; ++i)
          if (dnx0[i] != 0.0 || dnx1[i] != 0.0)
            KROME_JEL(J, i, j) += (dnx1[i] - dnx0[i]) / h;
      }
      for (int b = 0; b < c->B; ++b) {
        const double keep = c->Phi[b];
        const double du = 1e-7;  // d u_b of a relative Phi change of 1e-7
        c->Phi[b] = keep * (1.0 + du);
        secion_fixedT_(y, c, dnx1);
        c->Phi[b] = keep;
        for (int i = 0; i < NSPEC; ++i)
          if (dnx0[i] != 0.0 || dnx1[i] != 0.0)
            KROME_JEL(J, i, c->iu + b) += (dnx1[i] - dnx0[i]) / du;
      }
    }
#endif
    (void)NEQ; (void)yp; (void)fp;
    return 0;
  }

  // Tgas and Tdust columns: central differences of the full rhs
  double* fm = c->fm;
  for (int i = 0; i < NEQ; ++i) yp[i] = y[i];
  for (int q = 0; q < 1 + c->D; ++q) {
    const int col = q == 0 ? c->iT : c->id + q - 1;
    const double h = std::fmax(1e-6 * std::fabs(y[col]), 1e-8);
    yp[col] = y[col] + h;
    full_rhs_(yp, fp, c);
    yp[col] = y[col] - h;
    full_rhs_(yp, fm, c);
    yp[col] = y[col];
    for (int i = 0; i < NEQ; ++i)
      KROME_JEL(J, i, col) = (fp[i] - fm[i]) / (2.0 * h);
  }

  // aux rows (Tgas, Tdust, secondary-ionization fluxes) vs n_j: difference
  // quotient with the rates re-evaluated at the perturbed densities, so
  // k-bearing heating terms get the dk/dn of n-dependent rates
  if (c->iu >= 0) unpack_rad_(y, c);
  // the Tgas columns above left the shields at T -/+ h: back to y
  if constexpr (PHOTO_SHIELDED) {
#ifdef KROME_PHOTO_HPP_INCLUDED
    if (tc.photo_on) photo_shields(T, y, tc.p, tc.shield);
#endif
  }
  double* kp = c->kp;
  double* dnx0 = c->dnx0;
  double* dnx1 = c->dnx1;
  double* td0 = c->td0;
  double* td1 = c->td1;
  // the base state's rates (and dust balance / hooks: the calls above may
  // have left those of another state)
  eval_rates_(T, y, c, k0);
  const double f0 = aux_(y, k0, c, dnx0, td0);
  for (int j = 0; j < NSPEC; ++j) {
    const double h = fd_step_n_(y, j);
    yp[j] = y[j] + h;
    eval_rates_(T, yp, c, kp);
    const double f1 = aux_(yp, kp, c, dnx1, td1);
    yp[j] = y[j];
    KROME_JEL(J, c->iT, j) = (f1 - f0) / h;
    for (int b = 0; b < c->D; ++b)
      KROME_JEL(J, c->id + b, j) = (td1[b] - td0[b]) / h;
    for (int i = 0; i < NSPEC; ++i)
      if (dnx0[i] != 0.0 || dnx1[i] != 0.0)
        KROME_JEL(J, i, j) += (dnx1[i] - dnx0[i]) / h;
  }

  // aux rows vs u_b: photoheating, photoelectric Habing flux, secondary
  // fluxes and grain radiation absorption all respond to Phi
  for (int b = 0; b < c->B; ++b) {
    const double keep = c->Phi[b];
    const double du = 1e-7;  // d u_b of a relative Phi change of 1e-7
    c->Phi[b] = keep * (1.0 + du);
    eval_rates_(T, y, c, kp);
    const double f1 = aux_(y, kp, c, dnx1, td1);
    c->Phi[b] = keep;
    KROME_JEL(J, c->iT, c->iu + b) += (f1 - f0) / du;
    for (int q = 0; q < c->D; ++q)
      KROME_JEL(J, c->id + q, c->iu + b) += (td1[q] - td0[q]) / du;
    for (int i = 0; i < NSPEC; ++i)
      if (dnx0[i] != 0.0 || dnx1[i] != 0.0)
        KROME_JEL(J, i, c->iu + b) += (dnx1[i] - dnx0[i]) / du;
  }
  return 0;
}

// the analytic Jacobian the integrator uses (J zeroed by the caller):
// jac_assemble_, then the frozen species' rows zeroed
KROME_HD inline int jac_cell_(const double* y, DenseMat_ J, RhsCtx* c,
                              double* yp, double* fp) {
  const int flag = jac_assemble_(y, J, c, yp, fp);
  if constexpr (NFREEZE_MAX > 0) {  // frozen species: zero Jacobian rows
    int idx[NFREEZE_MAX > 0 ? NFREEZE_MAX : 1];
    const int m = frozen_species(c->th->p, idx);
    for (int q = 0; q < m; ++q)
      for (int j = 0; j < c->NEQ; ++j) KROME_JEL(J, idx[q], j) = 0.0;
  }
  return flag;
}

// the same into any J (zeroed here): for callers outside the integrator
inline int jac_host_(const double* y, DenseMat_ J, RhsCtx* c, double* t1,
                     double* t2) {
  for (int j = 0; j < c->NEQ; ++j)
    for (int i = 0; i < c->NEQ; ++i) KROME_JEL(J, i, j) = 0.0;
  return jac_cell_(y, J, c, t1, t2);
}

// ---- the host solver's context: RhsCtx over its own vectors -------------

struct SolveCtx : RhsCtx {
  ThermCtx tc;
  std::vector<double> ewt_r;     // per-slot rtol (rtol_species)
  std::vector<double> work_;     // the RhsCtx work arrays (rhs_bind_)
  std::vector<double> kfix_;     // rates_fixed (solve_caches_)
  std::vector<double> kT_;       // rates_fixed per Tgas (solve_caches_)
#ifdef KROME_PHOTO_HPP_INCLUDED
  PhotoSums psum;
  PhotoSumsView psumv;
#endif
  SolveCtx() { th = &tc; }
  SolveCtx(const SolveCtx&) = delete;
  SolveCtx& operator=(const SolveCtx&) = delete;
};

// build the state layout and work buffers for the given options; returns
// 0, or -102 / -103 / -105 for inconsistent inputs (see solve)
inline int init_ctx_(SolveCtx& c, double Tgas, const Params& params,
                     const SolverOpts& o) {
  if (params_missing(params) >= 0) return -105;
  c.tc.p = params;
  c.Tfix = Tgas;
  c.NEQ = NSPEC;
  if (o.evolve_T) c.iT = c.NEQ++;
#ifdef KROME_PHOTO_HPP_INCLUDED
  c.tc.photo = o.photo;
  c.tc.photo_on = o.photo != nullptr;
  if (o.photo) c.tc.ph = o.photo->view();
  if (o.photo) {
    if (o.Phi == nullptr) return -103;
    if constexpr (EvolveConfig::rad) {
      c.B = o.photo->B;
      c.iu = c.NEQ;
      c.NEQ += c.B;
    } else {
      c.tc.Phi = o.Phi;
    }
  } else if (EvolveConfig::rad) {
    return -103;  // generated to co-evolve a field, none supplied
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  // surface.bins: the per-bin surface species need that many bins
  if constexpr (SURF_NB > 0)
    if (o.dust == nullptr || o.dust->nb != SURF_NB) return -107;
  c.tc.dust_on = o.dust != nullptr;
  if (o.dust != nullptr) {
    // the population as the physics sees it: its own radii and counts
    // (solve() points the radii at the state's with dust.evolve)
    c.tc.dust = o.dust->view();
    c.tc.dust.acur = nullptr;
    c.tc.dust.xds_fix = -1.0;
    c.tc.dust_Td = o.dust->Td.data();
  }
  // configured binned-dust cooling needs its population
  if (TC::cool_dust_bins && o.evolve_T && o.dust == nullptr) return -102;
  // Tdust as state (dust.evolve_Td) only while Tgas evolves: with thermal
  // evolution off the dust temperatures stay fixed as well
  if (EvolveConfig::Td && o.evolve_T) {
    c.D = o.dust->nb;
    c.id = c.NEQ;
    c.NEQ += c.D;
  }
  // dust.evolve: the grain radii as state (thermal evolution or not);
  // the element conservation (conserve 2) would undo the exchange with
  // the grains
  if constexpr (EvolveConfig::dust) {
    if (o.dust == nullptr) return -102;
    if (o.conserve == 2) return -106;
    c.A = o.dust->nb;
    c.ia = c.NEQ;
    c.NEQ += c.A;
  }
#endif
  // the work arrays; without dust.evolve_Td the bins' Td are always the
  // balance of the state (a network whose physics uses the bins has the
  // d2g parameter)
  int nb = 0;
  bool balance = false;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  nb = o.dust != nullptr ? o.dust->nb : 0;
#ifdef KROME_HAS_PAR_d2g
  balance = !EvolveConfig::Td && o.dust != nullptr;
#endif
#endif
  c.work_.assign(rhs_work_doubles(c.NEQ, c.B, nb, c.D, c.A), 0.0);
  if (rhs_bind_(&c, c.work_.data(), nb, balance)
      != c.work_.data() + c.work_.size())
    std::abort();  // rhs_work_doubles and rhs_bind_ disagree
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (c.B > 0) {  // the evolved field starts at the caller's
    for (int b = 0; b < c.B; ++b) c.Phi[b] = o.Phi[b];
    c.tc.Phi = c.Phi;
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (o.dust != nullptr) dust_area_weights(c.tc.dust, c.awcur);
  for (int b = 0; b < c.A; ++b) c.arad[b] = o.dust->a[b];
#endif
  for (int r = 0; r <= NPHOTO; ++r) c.tc.shield[r] = 1.0;
  for (int q = 0; q <= NLW; ++q) c.tc.khat[q] = 0.0;
  return 0;
}

// pack n, Tgas, the evolved field and dust temperatures into y[NEQ]
inline void pack_state_(const SolveCtx& c, const double* n, double Tgas,
                        const SolverOpts& o, double* y) {
  for (int i = 0; i < NSPEC; ++i) y[i] = n[i];
  if (c.iT >= 0) y[c.iT] = Tgas;
  for (int b = 0; b < c.B; ++b)
    y[c.iu + b] = std::log(std::fmax(c.Phi[b], 1e-300));
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  for (int b = 0; b < c.D; ++b) y[c.id + b] = o.dust->Td[b];
  for (int b = 0; b < c.A; ++b)
    y[c.ia + b] = o.dust->atoms_of_radius(b, o.dust->a[b], c.xds);
#else
  (void)o;
#endif
}

// the per-cell state the helpers inside the rates read (Params::hk), set
// in the context's own copy of the parameters: n_global = the abundances
// at the start of the call (the input n is only written back after the
// integration), the photon field of the photobin-based helpers, and the
// binned-dust hook of a supplied population (surface rates, table_Tdust)
inline void install_hooks_(SolveCtx& c, const double* n0,
                           const SolverOpts& o) {
  StateHooks& h = c.tc.p.hk;
  h = StateHooks();
  h.n_global = n0;
#ifdef KROME_PHOTO_HPP_INCLUDED
  h.phi = c.tc.Phi;
  h.nbins = o.photo ? o.photo->B : 0;
  h.elo = o.photo ? o.photo->Elo.data() : nullptr;
  h.ehi = o.photo ? o.photo->Ehi.data() : nullptr;
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (o.dust != nullptr) {
    o.dust->hook_install(c.tc.p);
    c.tc.p.hk.dust.aw = c.awcur;  // the area weights of the state's radii
  }
#endif
  (void)o;
}

// per-solve caches of what stays constant during one solve() (exact: the
// same numbers, computed once): at fixed Tgas the rates of Tgas and the
// parameters (rates_fixed at the capped temperature), and with a static
// photon field (no radiation.evolve) its per-reaction bin sums. Not used by
// the queries (the device path keeps its own per-cell sums).
inline void solve_caches_(SolveCtx& c, const double* n) {
  if constexpr (NREA_FIXED > 0) {
    if (c.iT < 0) {
      c.kfix_.assign(NREA + 1, 0.0);
      rates_fixed(state_T_(n, &c), n, c.tc.p, c.kfix_.data());
      c.kfix = c.kfix_.data();
    }
#ifndef KROME_NO_RATE_REUSE
    // evolved Tgas: rates_fixed reused while Tgas is the same (exact; not
    // with binned dust, whose Td-dependent rates follow the state)
    else if (!c.tc.dust_on && c.tdb == nullptr && c.id < 0 && c.ia < 0
             && c.tc.dust.nb == 0) {
      c.kT_.assign(NREA + 1, 0.0);
      c.kT = c.kT_.data();
      c.kT_ok = false;
    }
#endif
  }
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (c.tc.photo_on && c.iu < 0) {
    c.tc.photo->sums(c.tc.Phi, c.psum);
    c.psumv = c.psum.view();
    c.tc.psum = &c.psumv;
  }
#endif
  (void)c; (void)n;
}

// post-solve conservation of cell n against its input n0 (SolverOpts::
// conserve). Mode 1: electrons restore the input net charge. Mode 2: the
// minimal relative correction delta_s = w_s sum_k A_ks lambda_k, w_s = n_s
// (a species at zero stays zero), that restores the input element totals
// and charge (rows A: species_elements, species_charge); frozen and opaque
// species do not move. The small system (A W A^T) lambda = residual is
// solved by Gaussian elimination with partial pivoting (dependent rows are
// skipped); projection and clipping are iterated a few times. Mode 3:
// the electrons neutralize the ions, n(e-) = max(sum q n over the other
// species, 0), even while a freeze directive holds e- fixed (legacy
// -computeElectrons, applied after the solve instead of in the rhs).
inline void conserve_cell_(const double* n0, double* n, const Params& p,
                           int mode) {
  if (mode <= 0) return;
  if (mode == 3) {
#ifdef KROME_HAS_E
    double qi = 0.0;
    for (int i = 0; i < NSPEC; ++i)
      if (i != IDX_E) qi += species_charge[i] * n[i];
    n[IDX_E] = std::fmax(qi, 0.0);
#endif
    (void)n0; (void)p;
    return;
  }
  bool fixed[NSPEC];
  for (int i = 0; i < NSPEC; ++i) fixed[i] = species_opaque[i];
  if constexpr (NFREEZE_MAX > 0) {
    int idx[NFREEZE_MAX > 0 ? NFREEZE_MAX : 1];
    const int m = frozen_species(p, idx);
    for (int q = 0; q < m; ++q) fixed[idx[q]] = true;
  }
  (void)p;
  double q0 = 0.0, q = 0.0;
  for (int i = 0; i < NSPEC; ++i) {
    q0 += species_charge[i] * n0[i];
    q += species_charge[i] * n[i];
  }
  (void)q0; (void)q;  // unused without electrons
  if (mode == 1) {
#ifdef KROME_HAS_E
    if (!fixed[IDX_E]) n[IDX_E] = std::fmax(n[IDX_E] + (q - q0), 0.0);
#endif
    return;
  }
  constexpr int NR = NELEM + 1;  // elements + charge
  auto A = [](int k, int i) -> double {
    return k < NELEM ? species_elements[i * NELEM + k] : species_charge[i];
  };
  for (int it = 0; it < 3; ++it) {
    double r[NR], M[NR][NR];
    for (int k = 0; k < NR; ++k) {
      double t0 = 0.0, t1 = 0.0;
      for (int i = 0; i < NSPEC; ++i) {
        t0 += A(k, i) * n0[i];
        t1 += A(k, i) * n[i];
      }
      r[k] = t0 - t1;
      for (int l = 0; l < NR; ++l) {
        double s = 0.0;
        for (int i = 0; i < NSPEC; ++i)
          if (!fixed[i]) s += A(k, i) * A(l, i) * n[i];
        M[k][l] = s;
      }
    }
    // Gaussian elimination with partial pivoting; rows without freedom
    // (or dependent on others) get lambda = 0
    double lam[NR];
    for (int k = 0; k < NR; ++k) lam[k] = 0.0;
    double scale = 0.0;
    for (int k = 0; k < NR; ++k) scale = std::fmax(scale, std::fabs(M[k][k]));
    bool used[NR] = {};
    for (int col = 0; col < NR; ++col) {
      int piv = -1;
      double best = 1e-14 * scale;
      for (int row = col; row < NR; ++row)
        if (std::fabs(M[row][col]) > best) best = std::fabs(M[row][col]),
                                           piv = row;
      if (piv < 0) continue;
      if (piv != col) {
        for (int l = 0; l < NR; ++l) std::swap(M[piv][l], M[col][l]);
        std::swap(r[piv], r[col]);
      }
      used[col] = true;
      for (int row = col + 1; row < NR; ++row) {
        const double f = M[row][col] / M[col][col];
        for (int l = col; l < NR; ++l) M[row][l] -= f * M[col][l];
        r[row] -= f * r[col];
      }
    }
    for (int col = NR - 1; col >= 0; --col) {
      if (!used[col]) continue;
      double s = r[col];
      for (int l = col + 1; l < NR; ++l) s -= M[col][l] * lam[l];
      lam[col] = s / M[col][col];
    }
    bool clipped = false;
    for (int i = 0; i < NSPEC; ++i) {
      if (fixed[i]) continue;
      double d = 0.0;
      for (int k = 0; k < NR; ++k) d += A(k, i) * lam[k];
      n[i] += n[i] * d;
      if (n[i] < 0.0) n[i] = 0.0, clipped = true;
    }
    if (!clipped) break;
  }
}


// the host integrator's problem: the solver's rhs and Jacobian
// (rhs_host_, jac_cell_; the integrator zeroes J) on its arrays
struct HostProb_ {
  SolveCtx* c;
  int rhs(double, const double* y, double* f) { return rhs_host_(y, f, c); }
  int jac(double, const double* y, const double*, double* J, double* t1,
          double* t2, double*) {
    return jac_cell_(y, DenseMat_(J, c->NEQ), c, t1, t2);
  }
};

// the sparse LU's tables (generated: Codegen.sparse_lu_tables) as the
// integrator's descriptor; built where it is used (on a device, from the
// tables' device accessors)
KROME_HD inline bdf::Sparse sparse_tables_() {
  bdf::Sparse s;
  s.ns = NSPEC;
  s.nnz = SPARSE_NNZ;
  s.perm = sp_perm;
  s.rp = sp_rp;
  s.ci = sp_ci;
  s.di = sp_di;
  return s;
}

// whether the sparse pattern holds every entry the analytic Jacobian
// assembles: the difference-quotient blocks of photo shielding, J_LW
// rates, secondary ionization and dust exchange add species entries
// outside it, which the sparse LU drops (an approximate Newton matrix:
// the solution is unchanged, the Newton convergence may be slower)
constexpr bool sparse_exact_() {
#ifdef KROME_HAS_SECION
  return false;
#else
  return !EvolveConfig::dust && !PHOTO_SHIELDED && NLW == 0;
#endif
}

// LINSOL_AUTO's choice of the sparse LU: an exact pattern, a factorization
// of at most a tenth of the dense one's flops (the fill decides, not the
// size: measured, 55-species ions 3x and 452-species cloud 2x faster than
// LAPACK, the 9-34-species benchmark networks, 45-87% filled, up to 1.5x
// slower), and a border of at most SPARSE_MAX_BORDER equations (Tgas, a
// few aux slots: photon or dust bins make the dense border dominate)
constexpr int SPARSE_MAX_BORDER = 4;
KROME_HD inline bool sparse_auto_(int NEQ) {
  return sparse_exact_() && NEQ - NSPEC <= SPARSE_MAX_BORDER
         && 2.0 * NSPEC * NSPEC * NSPEC / 3.0 >= 10.0 * SPARSE_FLOPS;
}

// the linear solver of a host solve (SolverOpts::linsol resolved): AUTO is
// the sparse LU where sparse_auto_, else LAPACK from LAPACK_MIN_NEQ
// equations in a build with it, else the dense LU
inline LinSol linsol_of_(LinSol l, int NEQ) {
  if (l != LINSOL_AUTO) return l;
  if (sparse_auto_(NEQ)) return LINSOL_SPARSE;
  return HAS_LAPACK && NEQ >= LAPACK_MIN_NEQ ? LINSOL_LAPACK : LINSOL_DENSE;
}

// the charge as the invariant the Newton iteration keeps with the
// difference-quotient Jacobian (bdf::Opts::inv_row, the electrons'
// component of each correction): without it the quotients of trace
// species leak charge at roundoff, which a collapse amplifies until it
// exceeds the electrons (0.2-0.4 % of the cells of the collapse benchmark
// lost them, 2026-10-07). Not with dust.gas_total (the slots hold totals) or
// freezable species (a frozen ion breaks the invariant);
// -DKROME_NO_CHARGE_NEWTON: off (CVODE's own iteration, the fixtures of
// tests/test_bdf_cell.py)
KROME_HD inline void charge_newton_(bdf::Opts& o) {
#if defined(KROME_HAS_E) && !defined(KROME_GAS_TOTAL) \
    && !defined(KROME_NO_CHARGE_NEWTON)
  if constexpr (NFREEZE_MAX == 0) {
    o.inv_row = IDX_E;
    o.inv_n = NSPEC;
    o.inv_c = species_charge;
  }
#else
  (void)o;
#endif
}

// the integrator's workspace, kept per thread and reused across calls
// (memory only: every integration starts from a fresh state)
struct HostWork_ {
  std::vector<double> work;
  std::vector<int> iwork;
  bool busy = false;
};
inline thread_local HostWork_ host_work_;

// integrate y (NEQ values, in place) over [0, dt] with the integrator
// (bdf_cell.hpp: CVODE's BDF, ported); error-test, convergence and
// repeated-rhs failures restart at the reached state (the integrator
// returns the last accepted state on failure), from a tenth of the step
// size it had reduced to, up to SolverOpts::max_restarts times: a restart
// from the initial-step estimate would repeat the failed attempts (a
// thermal collapse that needs a far smaller first step, legacy
// collapseCO). atol: the absolute tolerance per slot. Returns 0 or an
// integrator flag; diag (optional): the final state, for the error report.
inline int integrate_(SolveCtx& c, double* y, const double* atol, double dt,
                      const SolverOpts& o, bool dq_jac, LinSol ls,
                      SolveStats* st, bdf::Diag* diag) {
  const int N = c.NEQ;
  HostWork_ tmp;  // used only when the thread's workspace is busy
  HostWork_& W = host_work_.busy ? tmp : host_work_;
  W.busy = true;
  const bool sparse = ls == LINSOL_SPARSE;
  const long wd = bdf::work_doubles(N, sparse, NSPEC, SPARSE_NNZ);
  if (static_cast<long>(W.work.size()) < wd) W.work.resize(wd);
  if (static_cast<long>(W.iwork.size()) < bdf::work_ints(N))
    W.iwork.resize(bdf::work_ints(N));
  HostProb_ P{&c};
  bdf::Opts bo;
  bo.rtol = o.rtol;
  bo.atol = atol;
  bo.rtol_v = rtol_varies_(o) ? c.ewt_r.data() : nullptr;
  bo.max_steps = o.max_steps;
  bo.max_err_fails = o.max_err_fails;
  bo.max_conv_fails = o.max_conv_fails;
  bo.dq_jac = dq_jac;
  charge_newton_(bo);
  bo.lapack = ls == LINSOL_LAPACK;
  bo.sparse = sparse;
  if (sparse) bo.sp = sparse_tables_();
  bdf::Stats bs;
  int nr = 0;
  const int flag = bdf::integrate(P, N, y, dt, bo, o.max_restarts,
                                  W.work.data(), W.iwork.data(), &bs, &nr, 1,
                                  diag);
  W.busy = false;
  if (st != nullptr) {
    st->steps += bs.nst;
    st->rhs_evals += bs.nfe;
    st->jac_evals += bs.nje;
    st->rhs_evals_jac += bs.nfeDQ;
    st->lin_setups += bs.nsetups;
    st->newton_iters += bs.nni;
    st->newton_fails += bs.nnf;
    st->err_fails += bs.netf;
    st->restarts += nr;
  }
  return flag;
}

// ---- the error report (legacy krome_dump -> KROME_ERROR_REPORT) ---------

inline std::mutex error_report_mutex_;  // one writer at a time (OpenMP)
// the cell being integrated by solve_cells (for the report), -1 outside
inline thread_local int report_cell_ = -1;
inline thread_local int report_ncell_ = 0;

// write the report of a failed solve: the flag, the state the integrator
// returned (the last successful one), dn/dt there and the input state,
// copy-paste initial conditions, the rate coefficients and the reaction
// fluxes sorted by size, and the solver statistics. Returns 0, or -1 when
// the file cannot be opened.
inline int write_error_report_(const bdf::Diag& d, const double* y,
                               SolveCtx& c, const double* n_in, double T_in,
                               double t_end, const SolverOpts& o, int flag,
                               const SolveStats& st) {
  if (o.error_file == nullptr || o.error_file[0] == '\0') return -1;
  const int NEQ = c.NEQ;
  std::vector<double> f(NEQ, 0.0), k(NREA + 1, 0.0), fl(NREA + 1, 0.0);
  full_rhs_(y, f.data(), &c);
  const double T = state_T_(y, &c);
  eval_rates_(T, y, &c, k.data());
  fluxes(y, k.data(), fl.data());
  // the integrator's view of the failure
  const double tcur = d.tcur, hlast = d.hlast, hcur = d.hcur,
               tolsf = d.tolsf;
  const int qlast = d.qlast, qcur = d.qcur;
  const long netf = d.last.netf, ncfn = d.last.nnf;
  int imax = -1;  // component with the largest weighted local error
  if (d.ewt != nullptr && d.ele != nullptr) {
    double best = -1.0;
    for (int i = 0; i < NEQ; ++i) {
      const double e = std::fabs(d.ele[i] * d.ewt[i]);
      if (e > best) best = e, imax = i;
    }
  }
  auto slot_name = [&](int i) -> const char* {
    if (i < 0) return "?";
    if (i < NSPEC) return species_names[i];
    if (i == c.iT) return "Tgas";
    if (c.iu >= 0 && i >= c.iu && i < c.iu + c.B) return "ln Phi (bin)";
    if (c.id >= 0 && i >= c.id && i < c.id + c.D) return "Tdust (bin)";
    if (c.ia >= 0 && i >= c.ia && i < c.ia + c.A) return "grain radius (bin)";
    return "?";
  };
  std::vector<int> ord(NREA);
  for (int r = 0; r < NREA; ++r) ord[r] = r;

  // numbers as C++ / Fortran literals (a replay must compile)
  auto cnum = [](double v) -> std::string {
    if (std::isnan(v)) return "std::numeric_limits<double>::quiet_NaN()";
    if (std::isinf(v))
      return v > 0 ? "std::numeric_limits<double>::infinity()"
                   : "-std::numeric_limits<double>::infinity()";
    char b[40];
    std::snprintf(b, sizeof b, "%.17e", v);
    return b;
  };
  auto fnum = [](double v) -> std::string {
    if (std::isnan(v)) return "ieee_value(0d0, ieee_quiet_nan)";
    if (std::isinf(v))
      return v > 0 ? "ieee_value(0d0, ieee_positive_inf)"
                   : "ieee_value(0d0, ieee_negative_inf)";
    char b[40];
    std::snprintf(b, sizeof b, "%.17e", v);
    std::string t = b;
    t[t.find('e')] = 'd';
    return t;
  };
  const char* jac_names[] = {"JAC_AUTO", "JAC_ANALYTIC", "JAC_FD"};
  const char* ls_names[] = {"LINSOL_AUTO", "LINSOL_DENSE", "LINSOL_LAPACK",
                            "LINSOL_SPARSE"};
  const PhotoModule* ph = o.photo != nullptr && o.Phi != nullptr ? o.photo
                                                                  : nullptr;
  const DustBins* du = o.dust;
  Params pp = c.tc.p;
  const double rho_in = rho_(n_in), ntot_in = ntot_(n_in);
  double q_in = 0.0;
  for (int i = 0; i < NSPEC; ++i) q_in += species_charge[i] * n_in[i];

  std::lock_guard<std::mutex> lock(error_report_mutex_);
  FILE* fp = std::fopen(o.error_file, "w");
  if (fp == nullptr) return -1;
  std::fprintf(fp, "KROME ERROR REPORT (network %s)\n\n", NETWORK_NAME);
  std::fprintf(fp, "solve() failed with flag %d (%s) after %d restart(s)\n",
               flag, flag_name(flag), st.restarts);
  std::fprintf(fp, "  %s\n", flag_message(flag));
  if (report_cell_ >= 0)
    std::fprintf(fp, "cell %d of %d (0-based) of a solve_cells / "
                 "krome_solve_cells call\n", report_cell_, report_ncell_);
  std::fprintf(fp, "integration time reached %.6e s of %.6e s\n", tcur,
               t_end);

  // ---- what the caller passed
  std::fprintf(fp, "\nCell properties (as passed to KROME)\n");
  std::fprintf(fp, "  %-24s %.6e s\n", "t_end (dt)", t_end);
  std::fprintf(fp, "  %-24s %.6e K (%s)\n", "Tgas", T_in,
               c.iT >= 0 ? "evolved" : "fixed");
  std::fprintf(fp, "  %-24s %.6e g cm^-3\n", "rho", rho_in);
  std::fprintf(fp, "  %-24s %.6e cm^-3\n", "ntot", ntot_in);
  std::fprintf(fp, "  %-24s %.6e cm^-3\n", "H nuclei", Hnuclei_(n_in));
  std::fprintf(fp, "  %-24s %.6e\n", "mu",
               ntot_in > 0.0 ? rho_in / (ntot_in * constants::p_mass) : 0.0);
  std::fprintf(fp, "  %-24s %.6e e cm^-3\n", "net charge", q_in);
  std::fprintf(fp, "  parameters\n");
  for (int i = 0; i < NPAR; ++i)
    std::fprintf(fp, "    %-22s %.6e%s\n", param_names[i],
                 *param_ptr(pp, i),
                 param_required_by[i][0] != '\0' ? "  (required)" : "");
  if (ph != nullptr) {
    std::fprintf(fp, "  photon field: %d bins, Habing flux G0 = %.6e\n",
                 ph->B, ph->habing(o.Phi));
    std::fprintf(fp, "    %5s %14s %14s %14s\n", "bin", "Elo [eV]",
                 "Ehi [eV]", "Phi [cm-2 s-1]");
    for (int b2 = 0; b2 < ph->B; ++b2)
      std::fprintf(fp, "    %5d %14.6e %14.6e %14.6e\n", b2, ph->Elo[b2],
                   ph->Ehi[b2], o.Phi[b2]);
  } else {
    std::fprintf(fp, "  photon field: none\n");
  }
  if (du != nullptr) {
    std::fprintf(fp, "  dust bins: %d (%d per type%s%s), a = %.3e-%.3e cm, "
                 "MRN slope %.3f\n", du->nb, du->nd, du->useSi ? ", Si" : "",
                 du->useC ? ", C" : "", du->alow, du->aup, du->phi);
    for (int b2 = 0; b2 < du->nb && b2 < static_cast<int>(du->Td.size());
         ++b2)
      std::fprintf(fp, "    bin %3d  a = %.4e cm  Td = %.6e K\n", b2,
                   du->asize[b2], du->Td[b2]);
  } else {
    std::fprintf(fp, "  dust bins: none\n");
  }
  std::fprintf(fp, "  solver options\n");
  std::fprintf(fp, "    evolve_T %d, jac %s, linsol %s\n", o.evolve_T ? 1 : 0,
               jac_names[o.jac], ls_names[o.linsol]);
  std::fprintf(fp, "    rtol %.3e, atol_n %.3e, atol_T %.3e, atol_u %.3e\n",
               o.rtol, o.atol_n, o.atol_T, o.atol_u);
  if (o.atol_species != nullptr || o.rtol_species != nullptr)
    std::fprintf(fp, "    per-species tolerances:%s%s (see the replay)\n",
                 o.atol_species != nullptr ? " atol" : "",
                 o.rtol_species != nullptr ? " rtol" : "");
  std::fprintf(fp, "    max_steps %ld, max_err_fails %d, max_conv_fails %d, "
               "max_restarts %d\n", o.max_steps, o.max_err_fails,
               o.max_conv_fails, o.max_restarts);
  std::fprintf(fp, "    clip_negative %d, conserve %d\n",
               o.clip_negative ? 1 : 0, o.conserve);

  // ---- the state at the failure
  std::fprintf(fp, "\nSpecies (state returned by the integrator at the "
               "failure)\n");
  std::fprintf(fp, "%5s %-20s %14s %14s %14s\n", "#", "name", "n [cm^-3]",
               "dn/dt", "n input");
  for (int i = 0; i < NSPEC; ++i)
    std::fprintf(fp, "%5d %-20s %14.6e %14.6e %14.6e\n", i, species_names[i],
                 y[i], f[i], n_in[i]);
  if (c.iT >= 0)
    std::fprintf(fp, "%5d %-20s %14.6e %14.6e %14.6e\n", c.iT, "Tgas",
                 y[c.iT], f[c.iT], T_in);

  // ---- the failing call, ready to replay
  std::fprintf(fp, "\nReplay, C++ (a driver that includes krome.hpp; "
               "0-based indices)\n");
  std::fprintf(fp, "  using namespace krome;\n");
  std::fprintf(fp, "  double n[NSPEC] = {0.0};\n");
  for (int i = 0; i < NSPEC; ++i)
    if (n_in[i] != 0.0)
      std::fprintf(fp, "  n[%d] = %s;  // %s\n", i, cnum(n_in[i]).c_str(),
                   species_names[i]);
  std::fprintf(fp, "  double Tgas = %s;\n", cnum(T_in).c_str());
  std::fprintf(fp, "  const double t_end = %s;\n", cnum(t_end).c_str());
  std::fprintf(fp, "  Params p;\n");
  for (int i = 0; i < NPAR; ++i)
    std::fprintf(fp, "  p.%s = %s;\n", param_names[i],
                 cnum(*param_ptr(pp, i)).c_str());
  std::fprintf(fp, "  SolverOpts o;\n");
  std::fprintf(fp, "  o.evolve_T = %s;\n", o.evolve_T ? "true" : "false");
  std::fprintf(fp, "  o.jac = %s;\n  o.linsol = %s;\n", jac_names[o.jac],
               ls_names[o.linsol]);
  std::fprintf(fp, "  o.rtol = %s;\n  o.atol_n = %s;\n  o.atol_T = %s;\n"
               "  o.atol_u = %s;\n", cnum(o.rtol).c_str(),
               cnum(o.atol_n).c_str(), cnum(o.atol_T).c_str(),
               cnum(o.atol_u).c_str());
  for (int q = 0; q < 2; ++q) {  // per-species tolerances
    const double* v = q == 0 ? o.atol_species : o.rtol_species;
    const char* nm = q == 0 ? "atol_species" : "rtol_species";
    if (v == nullptr) {
      std::fprintf(fp, "  o.%s = nullptr;\n", nm);
      continue;
    }
    std::fprintf(fp, "  static const double %s_[NSPEC] = {", nm);
    for (int i = 0; i < NSPEC; ++i)
      std::fprintf(fp, "%s%s", i ? ", " : "", cnum(v[i]).c_str());
    std::fprintf(fp, "};\n  o.%s = %s_;\n", nm, nm);
  }
  std::fprintf(fp, "  o.max_steps = %ld;\n  o.max_err_fails = %d;\n"
               "  o.max_conv_fails = %d;\n  o.max_restarts = %d;\n",
               o.max_steps, o.max_err_fails, o.max_conv_fails,
               o.max_restarts);
  std::fprintf(fp, "  o.clip_negative = %s;\n  o.conserve = %d;\n",
               o.clip_negative ? "true" : "false", o.conserve);
  if (ph != nullptr) {
    std::fprintf(fp, "  Radiation rad;\n");
    const char* nm[3] = {"Elo", "Ehi", "Phi"};
    for (int v = 0; v < 3; ++v) {
      std::fprintf(fp, "  rad.%s = {", nm[v]);
      for (int b2 = 0; b2 < ph->B; ++b2) {
        const double x = v == 0 ? ph->Elo[b2] : v == 1 ? ph->Ehi[b2]
                                                       : o.Phi[b2];
        std::fprintf(fp, "%s%s", b2 ? ", " : "", cnum(x).c_str());
      }
      std::fprintf(fp, "};\n");
    }
    std::fprintf(fp, "  PhotoModule photo;\n  photo.init(rad);\n"
                 "  o.photo = &photo;\n  o.Phi = rad.Phi.data();\n");
  }
  if (du != nullptr) {
    std::fprintf(fp, "  DustBins dust;\n  dust.nd = %d;\n  dust.useSi = %s;\n"
                 "  dust.useC = %s;\n  dust.alow = %s;\n  dust.aup = %s;\n"
                 "  dust.phi = %s;\n  dust.init();\n", du->nd,
                 du->useSi ? "true" : "false", du->useC ? "true" : "false",
                 cnum(du->alow).c_str(), cnum(du->aup).c_str(),
                 cnum(du->phi).c_str());
    if (ph != nullptr && !du->qmid.empty())
      std::fprintf(fp, "  dust.interp_qabs(rad);\n");
    std::fprintf(fp, "  dust.Td = {");
    for (size_t b2 = 0; b2 < du->Td.size(); ++b2)
      std::fprintf(fp, "%s%s", b2 ? ", " : "", cnum(du->Td[b2]).c_str());
    std::fprintf(fp, "};\n");
    if (EvolveConfig::dust) {  // the input grain radii
      std::fprintf(fp, "  dust.a = {");
      for (size_t b2 = 0; b2 < du->a.size(); ++b2)
        std::fprintf(fp, "%s%s", b2 ? ", " : "", cnum(du->a[b2]).c_str());
      std::fprintf(fp, "};\n  dust.set_radii(nullptr);\n");
    }
    std::fprintf(fp, "  o.dust = &dust;\n");
  }
  std::fprintf(fp, "  const int flag = solve(n, Tgas, p, t_end, o);  "
               "// %d here\n", flag);
  if (FORTRAN_BINDING) {
    std::fprintf(fp, "\nReplay, Fortran (use krome, krome_net; 1-based "
                 "indices; one cell)\n");
    std::fprintf(fp, "  real(8) :: n(krome_nspec, 1), Tgas(1), dt(1), "
                 "par(krome_npar, 1)\n");
    if (ph != nullptr)
      std::fprintf(fp, "  real(8) :: elo(%d), ehi(%d), phi(%d, 1)\n", ph->B,
                   ph->B, ph->B);
    if (du != nullptr)
      std::fprintf(fp, "  real(8) :: td(%d)%s\n", static_cast<int>(du->Td.size()),
                   EvolveConfig::dust ? ", ad(size(td))" : "");
    if (o.atol_species != nullptr || o.rtol_species != nullptr)
      std::fprintf(fp, "  real(8) :: tol(krome_nspec)\n");
    std::fprintf(fp, "  type(c_ptr) :: h\n  integer :: ierr\n");
    std::fprintf(fp, "  n = 0d0\n");
    for (int i = 0; i < NSPEC; ++i)
      if (n_in[i] != 0.0)
        std::fprintf(fp, "  n(%d, 1) = %s  ! %s\n", i + 1,
                     fnum(n_in[i]).c_str(), species_names[i]);
    std::fprintf(fp, "  Tgas(1) = %s\n  dt(1) = %s\n", fnum(T_in).c_str(),
                 fnum(t_end).c_str());
    std::fprintf(fp, "  call krome_params_default(par)\n");
    for (int i = 0; i < NPAR; ++i)
      std::fprintf(fp, "  par(%d, 1) = %s  ! %s\n", i + 1,
                   fnum(*param_ptr(pp, i)).c_str(), param_names[i]);
    std::fprintf(fp, "  h = krome_opts_new()\n");
    std::fprintf(fp, "  call krome_set(h, \"evolve_T\", %s)\n",
                 o.evolve_T ? "1d0" : "0d0");
    std::fprintf(fp, "  call krome_set(h, \"jac\", %dd0)\n", o.jac);
    std::fprintf(fp, "  call krome_set(h, \"linsol\", %dd0)\n", o.linsol);
    const std::pair<const char*, double> reals[] = {
        {"rtol", o.rtol}, {"atol_n", o.atol_n}, {"atol_T", o.atol_T},
        {"atol_u", o.atol_u}};
    for (const auto& kv : reals)
      std::fprintf(fp, "  call krome_set(h, \"%s\", %s)\n", kv.first,
                   fnum(kv.second).c_str());
    for (int q = 0; q < 2; ++q) {  // per-species tolerances
      const double* v = q == 0 ? o.atol_species : o.rtol_species;
      if (v == nullptr) continue;
      for (int i = 0; i < NSPEC; ++i)
        std::fprintf(fp, "  tol(%d) = %s  ! %s\n", i + 1, fnum(v[i]).c_str(),
                     species_names[i]);
      std::fprintf(fp, "  call krome_set_%s_species(h, tol)\n",
                   q == 0 ? "atol" : "rtol");
    }
    const std::pair<const char*, long> ints[] = {
        {"max_steps", o.max_steps}, {"max_err_fails", o.max_err_fails},
        {"max_conv_fails", o.max_conv_fails},
        {"max_restarts", o.max_restarts},
        {"clip_negative", o.clip_negative ? 1 : 0},
        {"conserve", o.conserve}};
    for (const auto& kv : ints)
      std::fprintf(fp, "  call krome_set(h, \"%s\", %ldd0)\n", kv.first,
                   kv.second);
    if (du != nullptr) {
      std::fprintf(fp, "  call krome_set_dustbins(h, %d, %s, %s, %s, %s, %s)\n",
                   du->nd, du->useSi ? ".true." : ".false.",
                   du->useC ? ".true." : ".false.", fnum(du->alow).c_str(),
                   fnum(du->aup).c_str(), fnum(du->phi).c_str());
      for (size_t b2 = 0; b2 < du->Td.size(); ++b2)
        std::fprintf(fp, "  td(%d) = %s\n", static_cast<int>(b2) + 1,
                     fnum(du->Td[b2]).c_str());
      std::fprintf(fp, "  call krome_set_dust_Td(h, td)\n");
      if (EvolveConfig::dust) {
        for (size_t b2 = 0; b2 < du->a.size(); ++b2)
          std::fprintf(fp, "  ad(%d) = %s\n", static_cast<int>(b2) + 1,
                       fnum(du->a[b2]).c_str());
        std::fprintf(fp, "  call krome_set_dust_a(h, ad)\n");
      }
    }
    if (ph != nullptr) {
      for (int b2 = 0; b2 < ph->B; ++b2)
        std::fprintf(fp, "  elo(%d) = %s; ehi(%d) = %s; phi(%d, 1) = %s\n",
                     b2 + 1, fnum(ph->Elo[b2]).c_str(), b2 + 1,
                     fnum(ph->Ehi[b2]).c_str(), b2 + 1,
                     fnum(o.Phi[b2]).c_str());
      std::fprintf(fp, "  call krome_set_photo(h, elo, ehi)\n");
      std::fprintf(fp, "  ierr = krome_solve(n, Tgas, dt, h, par=par, "
                   "phi=phi)  ! %d here\n", flag);
    } else {
      std::fprintf(fp, "  ierr = krome_solve(n, Tgas, dt, h, par=par)  "
                   "! %d here\n", flag);
    }
  }

  // ---- rate coefficients and fluxes, sorted
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return std::fabs(k[a]) > std::fabs(k[b]); });
  const double kmax = NREA > 0 ? std::fabs(k[ord[0]]) : 0.0;
  std::fprintf(fp, "\nRate coefficients (sorted) at Tgas = %.6e K\n", T);
  std::fprintf(fp, "%5s %14s %10s  %s\n", "#", "k", "k %", "reaction");
  for (int q = 0; q < NREA; ++q) {
    const int r = ord[q];
    std::fprintf(fp, "%5d %14.6e %10.3e  %s\n", r, k[r],
                 kmax > 0.0 ? 100.0 * std::fabs(k[r]) / kmax : 0.0,
                 reaction_names[r]);
  }
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return std::fabs(fl[a]) > std::fabs(fl[b]); });
  const double fmax = NREA > 0 ? std::fabs(fl[ord[0]]) : 0.0;
  std::fprintf(fp, "\nReaction fluxes (sorted) [k n1 n2 ..., cm^-3 s^-1]\n");
  std::fprintf(fp, "%5s %14s %10s  %s\n", "#", "flux", "flux %", "reaction");
  for (int q = 0; q < NREA; ++q) {
    const int r = ord[q];
    std::fprintf(fp, "%5d %14.6e %10.3e  %s\n", r, fl[r],
                 fmax > 0.0 ? 100.0 * std::fabs(fl[r]) / fmax : 0.0,
                 reaction_names[r]);
  }

  // ---- the solver
  std::fprintf(fp, "\nSolver (BDF, CVODE's ported)\n");
  std::fprintf(fp, "  %-34s %.6e\n", "last step size", hlast);
  std::fprintf(fp, "  %-34s %.6e\n", "step size attempted", hcur);
  std::fprintf(fp, "  %-34s %.6e\n", "current time", tcur);
  std::fprintf(fp, "  %-34s %.6e\n", "tolerance scale factor", tolsf);
  std::fprintf(fp, "  %-34s %d / %d\n", "last / current order", qlast, qcur);
  std::fprintf(fp, "  %-34s %ld\n", "steps", st.steps);
  std::fprintf(fp, "  %-34s %ld\n", "rhs evaluations", st.rhs_evals);
  std::fprintf(fp, "  %-34s %ld\n", "Jacobian evaluations", st.jac_evals);
  std::fprintf(fp, "  %-34s %ld\n", "LU factorizations", st.lin_setups);
  std::fprintf(fp, "  %-34s %ld\n", "Newton iterations", st.newton_iters);
  std::fprintf(fp, "  %-34s %ld\n", "error-test failures (last run)", netf);
  std::fprintf(fp, "  %-34s %ld\n", "convergence failures (last run)", ncfn);
  std::fprintf(fp, "  %-34s %d (%s)\n", "largest local error", imax,
               slot_name(imax));
  std::fprintf(fp, "\nEND KROME ERROR REPORT\n");
  std::fclose(fp);
  return 0;
}

// integrate from t = 0 to t_end [s]. n[NSPEC] [cm^-3] is updated in place;
// Tgas is updated when opts.evolve_T (else it is the fixed temperature),
// opts.Phi with radiation.evolve, opts.dust->Td with dust.evolve_Td (and Tgas
// evolved). On failure (after the restarts, see SolverOpts)
// the inputs are left untouched. Returns 0 on success, a negative
// integrator flag on failure (FLAG_TOO_MUCH_WORK, ...), or
//   -102  cooling DUSTBINS configured, Tgas evolved and no DustBins
//   -103  photo module without a field Phi, or radiation.evolve without a
//         photo module
//   -105  a REQUIRED parameter is unset (NaN; see params_missing)
//   -106  conservation requested with the dust.gas_total representation
//   -107  surface.bins: no DustBins, or not SURF_NB bins
//   -108  linsol LINSOL_LAPACK in a build without LAPACK
inline int solve(double* n, double& Tgas, const Params& params,
                 double t_end, const SolverOpts& o, SolveStats* stats) {
  if (stats != nullptr) *stats = SolveStats();
#ifdef KROME_GAS_TOTAL
  if (o.conserve != 0) {  // slots store totals: not supported
    if (stats != nullptr) stats->flag = -106;
    return -106;
  }
#endif
  SolveCtx c;
  int flag = init_ctx_(c, Tgas, params, o);
  if (flag != 0) {
    if (stats != nullptr) stats->flag = flag;
    return flag;
  }
  const int NEQ = c.NEQ;
  if (o.linsol == LINSOL_LAPACK && !HAS_LAPACK) {
    if (stats != nullptr) stats->flag = -108;
    return -108;
  }
  const LinSol ls = linsol_of_(o.linsol, NEQ);
  install_hooks_(c, n, o);
#if defined(KROME_DUST_BINS_HPP_INCLUDED) && defined(KROME_HAS_PAR_d2g)
  // dust.evolve: the grain counts of the solve are those of its start
  // (rho_gas d2g), and the physics reads the radii from the state
  if (c.A > 0) {
    c.tc.dust.xds_fix = c.xds = rho_(n) * c.tc.p.d2g;
    c.tc.dust.acur = c.arad;
  }
#endif

  std::vector<double> yv(NEQ), avec(NEQ);
  double* y = yv.data();
  pack_state_(c, n, Tgas, o, y);
  solve_caches_(c, n);
  for (int i = 0; i < NSPEC; ++i)
    avec[i] = species_tol_(o.atol_species, i, o.atol_n);
  if (c.iT >= 0) avec[c.iT] = o.atol_T;
  for (int b = 0; b < c.B; ++b) avec[c.iu + b] = o.atol_u;
  for (int b = 0; b < c.D; ++b) avec[c.id + b] = o.atol_T;
  for (int b = 0; b < c.A; ++b) avec[c.ia + b] = o.atol_n;
  if (rtol_varies_(o)) {  // per-species rtol: weights with a per-slot rtol
    c.ewt_r.assign(NEQ, o.rtol);
    for (int i = 0; i < NSPEC; ++i)
      c.ewt_r[i] = species_tol_(o.rtol_species, i, o.rtol);
  }

  const JacMode jac = o.jac != JAC_AUTO ? o.jac
                      : (o.evolve_T ? JAC_FD : JAC_ANALYTIC);
  SolveStats st_local;  // counts for the error report without a caller's
  SolveStats* st = stats != nullptr ? stats : &st_local;
  // the final state of the integrator, for the error report
  std::vector<double> ewt_d, ele_d;
  bdf::Diag diag;
  if (o.error_report) {
    ewt_d.assign(NEQ, 0.0);
    ele_d.assign(NEQ, 0.0);
    diag.ewt = ewt_d.data();
    diag.ele = ele_d.data();
  }
  flag = integrate_(c, y, avec.data(), t_end, o, jac == JAC_FD, ls, st,
                    &diag);
  // a non-positive (or non-finite) final Tgas is a failure, whatever
  // the integrator reported: the inputs stay untouched
  if (flag == 0 && c.iT >= 0 && !(y[c.iT] > 0.0 && std::isfinite(y[c.iT])))
    flag = FLAG_NONPOSITIVE_T;
  // a failure after the rhs refused temperatures <= 0: say so
  if (c.hit_T0 && (flag == FLAG_REPTD_RHSFUNC_ERR || flag == FLAG_RHSFUNC_FAIL
                   || flag == FLAG_UNREC_RHSFUNC_ERR
                   || flag == FLAG_FIRST_RHSFUNC_ERR))
    flag = FLAG_NONPOSITIVE_T;
  if (flag != 0) st->flag = flag;
  if (flag != 0 && o.error_report)
    write_error_report_(diag, y, c, n, Tgas, t_end, o, flag, *st);
  if (flag == 0) {
    double n0[NSPEC];
    for (int i = 0; i < NSPEC; ++i) {
      n0[i] = n[i];
      n[i] = y[i];
      if (o.clip_negative && n[i] < 0.0) n[i] = 0.0;
    }
    conserve_cell_(n0, n, params, o.conserve);
    if (c.iT >= 0) Tgas = y[c.iT];
#ifdef KROME_PHOTO_HPP_INCLUDED
    for (int b = 0; b < c.B; ++b) o.Phi[b] = std::exp(y[c.iu + b]);
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
    for (int b = 0; b < c.D; ++b) o.dust->Td[b] = y[c.id + b];
    // the evolved radii (an evaporated bin: 0, never negative)
    for (int b = 0; b < c.A; ++b)
      o.dust->a[b] = o.dust->radius_of_atoms(b, y[c.ia + b], c.xds);
    if (c.A > 0) o.dust->set_radii(nullptr);  // area weights of the new a
    // balance Td: the one of the accepted final state, for the caller
#ifdef KROME_HAS_PAR_d2g
    if (c.tdb != nullptr) {
      const double* yf = y;
      if (c.iu >= 0) unpack_rad_(yf, &c);  // the final state's field
      const double* phi = nullptr;
#ifdef KROME_PHOTO_HPP_INCLUDED
      phi = c.tc.Phi;
#endif
      dust_balance(c.tc.dust, state_T_(yf, &c), yf, c.tc.p, phi,
                   o.dust->Td.data());
    }
#endif
#endif
  }
  if (stats != nullptr) stats->flag = flag;
  return flag;
}

// ---- mass fractions (legacy -useX) ---------------------------------------

// x_i = m_i n_i / rho with the caller's gas density rho [g cm^-3]: a
// species without mass (none in the usual networks) has no mass fraction
// and maps to n = 0 / x = 0. The solver itself works on number densities.
inline void x_to_n_(const double* x, double rho, double* n) {
  for (int i = 0; i < NSPEC; ++i)
    n[i] = species_mass_g[i] > 0.0 ? rho * x[i] / species_mass_g[i] : 0.0;
}
inline void n_to_x_(const double* n, double rho, double* x) {
  for (int i = 0; i < NSPEC; ++i) x[i] = species_mass_g[i] * n[i] / rho;
}
// the sum of the mass fractions of the species with a mass
inline double x_sum_(const double* x) {
  double s = 0.0;
  for (int i = 0; i < NSPEC; ++i)
    if (species_mass_g[i] > 0.0) s += x[i];
  return s;
}
// the solved x, rescaled to the input sum xin with o.renormalize_x
inline void n_to_x_out_(const double* n, double rho, double xin,
                        const SolverOpts& o, double* x) {
  n_to_x_(n, rho, x);
  if (!o.renormalize_x) return;
  const double s = x_sum_(x);
  if (s > 0.0)
    for (int i = 0; i < NSPEC; ++i) x[i] *= xin / s;
}

// solve() with mass fractions x[NSPEC] and the gas density rho: the same
// integration, converted on the way in (n = rho x / m) and out
// (x = m n / rho, the caller's rho; rescaled to the input sum with
// SolverOpts::renormalize_x; element and charge conservation is
// SolverOpts::conserve). On failure x is untouched.
inline int solve_x(double* x, double rho, double& Tgas, const Params& params,
                   double t_end, const SolverOpts& o, SolveStats* stats) {
  std::vector<double> n(NSPEC);
  const double xin = x_sum_(x);
  x_to_n_(x, rho, n.data());
  const int flag = solve(n.data(), Tgas, params, t_end, o, stats);
  if (flag == 0) n_to_x_out_(n.data(), rho, xin, o, x);
  return flag;
}

// ---------------------------------------------------------------------------
// equilibrium solver (legacy krome_equilibrium)


// integrate to steady state: outputs at t0 * factor^k, each interval a
// solve() from the current state (restarts, clipping, conservation as
// there); converged when max_s |dn_s| / max(n_s, eatol ntot) < ertol (and
// |dT|/T when Tgas evolves) between outputs. Isothermal unless
// opts.evolve_T. Returns 0 (converged), +1 (t_max reached first: the state
// at t_max is returned) or a negative flag (failure: inputs untouched).
// t_reached (optional) receives the time integrated.
inline int solve_equilibrium(double* n, double& Tgas, const Params& params,
                             const SolverOpts& o, const EquilOpts& e,
                             double* t_reached) {
  double nin[NSPEC], prev[NSPEC];
  for (int i = 0; i < NSPEC; ++i) nin[i] = n[i];
  const double Tin = Tgas;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  // the grain state the solves evolve (restored on failure)
  const std::vector<double> ain = o.dust != nullptr ? o.dust->a
                                                    : std::vector<double>();
  const std::vector<double> tin = o.dust != nullptr ? o.dust->Td
                                                    : std::vector<double>();
#endif
  double t = 0.0, dt = e.t0;
  for (;;) {
    for (int i = 0; i < NSPEC; ++i) prev[i] = n[i];
    const double Tprev = Tgas;
    const int flag = solve(n, Tgas, params, dt, o);
    if (flag != 0) {
      for (int i = 0; i < NSPEC; ++i) n[i] = nin[i];
      Tgas = Tin;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
      if (o.dust != nullptr) {
        o.dust->a = ain;
        o.dust->Td = tin;
        o.dust->set_radii(nullptr);
      }
#endif
      if (t_reached != nullptr) *t_reached = t;
      return flag;
    }
    t += dt;
    const double floor = e.eatol * ntot_(n);
    double w = o.evolve_T ? std::fabs(Tgas - Tprev) / Tgas : 0.0;
    for (int i = 0; i < NSPEC; ++i)
      w = std::fmax(w, std::fabs(n[i] - prev[i]) / std::fmax(n[i], floor));
    if (t_reached != nullptr) *t_reached = t;
    if (w < e.ertol) return 0;
    if (t >= e.t_max) return 1;
    dt = std::fmin((e.factor - 1.0) * t, e.t_max - t);
  }
}

} // namespace krome
