// KROME-X binned dust module: MRN size distribution per grain type
// (Si and/or C), per-bin dust temperatures from the radiative-collisional
// balance, and the resulting gas-dust cooling. Ported from legacy KROME
// src/krome_dust.f90 (set_dust_distribution, dust_load_Qabs,
// dust_init_intBB, besc, compute_Tdust).
//
//   balance per bin:  (intBB(Td) - intBB(Tcmb) - intJ) * besc
//                       = pre * (Tgas - Td),   pre = 0.5 fact vgas kB ntot
//   cooling = 4 pi Sum_bins (intBB(Td) - intBB(Tcmb) - intJ) besc xdust a^2
//
// Tdust is solved by bisection at every call (legacy default), or, with
// dust.evolve_Td, evolved per bin by the grains' energy equation
// (dust_dTdust below). Include the generated net_<name>.hpp BEFORE this
// header (uses rho_/ntot_); photo coupling (intJ) takes the Radiation bins
// + Phi.
//
// Physical constants from krome::constants (kromex/constants.py); grain
// bulk densities C 2.25 / Si 3.13 g/cm3 (Zhukovska), key-element
// fractions 10^(Asplund solar): C 8.43-12, Si 7.51-12.
#pragma once
#include <cmath>
#include <memory>
#include <vector>
#include "data_qabs.hpp"

#define KROME_DUST_BINS_HPP_INCLUDED 1

namespace krome {

// ---- the numerics of a dust population over plain arrays ---------------
//
// DustView: what the physics below reads, as pointers (host or device
// memory). DustBins (host) builds one over its own vectors (view()); the
// device path builds one per cell over device arrays. The functions are
// device-callable (KROME_HD) and are the only implementation: the
// DustBins methods call them.
struct DustView {
  int nb = 0;                     // bins
  int nT = 0;                     // intBB samples
  double TbbMin = 1.0, TbbMax = 1e4, TbbMult = 0.0;
  const double* xw = nullptr;     // [nb] xdust = xw * rho_gas * d2g
  const int* isC = nullptr;       // [nb] 1 carbon, 0 silicate
  const double* a = nullptr;      // [nb] the population's radii
  const double* acur = nullptr;   // [nb] the radii of a solve's state, or
                                  // nullptr (then a)
  double xds_fix = -1.0;          // the solve start's rho_gas d2g, or < 0
  const double* ibb = nullptr;    // [nb * nT] intBB tables
  const double* ibb_sig = nullptr;
  int B = 0;                      // photon bins of qmid / emid (0: none)
  const double* qmid = nullptr;   // [nb * B] Qabs at the bin midpoints
  const double* emid = nullptr;   // [B] the bin midpoints [eV]
};

// linear interpolation on the uniform Tbb grid (legacy get_dust_intBB).
// T is clamped into [TbbMin, TbbMax]: legacy indexed out of bounds
// there; the clamp keeps the dTdust state variant from running away
// when integration error pushes a Td outside the table.
KROME_HD inline double dust_tabT_(const DustView& v, const double* tab, int b,
                                  double T) {
  if (T < v.TbbMin) T = v.TbbMin;
  if (T > v.TbbMax) T = v.TbbMax;
  int i = static_cast<int>((T - v.TbbMin) * v.TbbMult);
  if (i < 0) i = 0;
  if (i > v.nT - 2) i = v.nT - 2;
  const double T1 = i / v.TbbMult + v.TbbMin;
  const double f = (T - T1) * v.TbbMult;
  const size_t o = static_cast<size_t>(b) * v.nT + i;
  return f * (tab[o + 1] - tab[o]) + tab[o];
}
KROME_HD inline double dust_intBB(const DustView& v, int b, double T) {
  return dust_tabT_(v, v.ibb, b, T);
}
// slope of the linearly interpolated intBB table at T (the exact
// derivative of what dust_intBB returns)
KROME_HD inline double dust_intBB_slope_(const DustView& v, int b, double T) {
  if (T < v.TbbMin) T = v.TbbMin;
  if (T > v.TbbMax) T = v.TbbMax;
  int i = static_cast<int>((T - v.TbbMin) * v.TbbMult);
  if (i < 0) i = 0;
  if (i > v.nT - 2) i = v.nT - 2;
  const size_t o = static_cast<size_t>(b) * v.nT + i;
  return (v.ibb[o + 1] - v.ibb[o]) * v.TbbMult;
}

// the current radius of bin b (the state during a solve, else a[b]; a
// negative trial value counts as 0) and its square
KROME_HD inline double dust_a_(const DustView& v, int b) {
  return std::fmax(v.acur != nullptr ? v.acur[b] : v.a[b], 0.0);
}
KROME_HD inline double dust_a2_(const DustView& v, int b) {
  const double x = dust_a_(v, b);
  return x * x;
}
// the area weights xw a^2 at the current radii into aw[nb]
KROME_HD inline void dust_area_weights(const DustView& v, double* aw) {
  for (int b = 0; b < v.nb; ++b) aw[b] = v.xw[b] * dust_a2_(v, b);
}

// Planck opacity of bin b (legacy kpla_dust) [cm^2 per grain]
KROME_HD inline double dust_kpla(const DustView& v, double Tdust, int b) {
  return dust_tabT_(v, v.ibb_sig, b, Tdust) * dust_a2_(v, b) * constants::pi;
}

// escape probability (legacy besc; gas opacity term disabled as legacy) at
// the per-bin dust temperatures Tds
KROME_HD inline double dust_besc(const DustView& v, double Tgas,
                                 double xd_scale, double lj,
                                 const double* Tds) {
  double tau_d = 0.0;
  for (int b = 0; b < v.nb; ++b) {
    double Tff = Tds[b];
    const double xd = v.xw[b] * xd_scale;
    if (xd < 1e-20) continue;
    if (std::fabs(Tff - Tgas) < 5.0) Tff = Tgas;
    tau_d += xd * dust_kpla(v, Tff, b);
  }
  const double tau = tau_d * lj;
  return (tau < 1.0) ? 1.0 : 1.0 / (tau * tau);
}

// shared prefactors of the exchange balance (templates over the parameter
// type: instantiated only by a network that uses the binned dust physics,
// which then has d2g and Tcmb)
template <class P>
KROME_HD inline void dust_pre(const DustView& v, double Tgas, const double* n,
                              double& pre, double& lj, double& xd_scale,
                              const P& p) {
  const double ntot = ntot_(n);
  const double rhogas = rho_(n);
  xd_scale = v.xds_fix > 0.0 ? v.xds_fix : rhogas * p.d2g;
  const double vgas = std::sqrt(8.0 * constants::boltzmann_erg / constants::pi
                                / constants::p_mass * Tgas);
  pre = 0.5 * 0.5 * vgas * constants::boltzmann_erg * ntot;
  const double rr = std::fmax(rhogas, 1e-40);
  const double mu = rr / std::fmax(ntot, 1e-40) / constants::p_mass;
  lj = std::sqrt(constants::pi * constants::boltzmann_erg * Tgas / rr
                 / constants::p_mass / constants::gravity / mu);
}

// absorbed radiation term intJ [erg cm^-2 s^-1] for bin b (legacy
// get_int_JQabs): Phi is the photon flux per bin; legacy pi J corresponds
// to Phi E / 4 (photon flux <-> intensity conversion)
KROME_HD inline double dust_intJ(const DustView& v, int b, const double* Phi) {
  if (Phi == nullptr || v.B == 0) return 0.0;
  double s = 0.0;
  for (int i = 0; i < v.B; ++i)
    s += Phi[i] * v.emid[i] * v.qmid[static_cast<size_t>(b) * v.B + i];
  return 0.25 * s * constants::eV_to_erg;
}

// dust cooling from the temperatures Tds (legacy cooling_dust)
template <class P>
KROME_HD inline double dust_cooling_state(const DustView& v, double Tgas,
                                          const double* n, const P& p,
                                          const double* Tds,
                                          const double* Phi) {
  double pre, lj, xds;
  dust_pre(v, Tgas, n, pre, lj, xds, p);
  const double be = dust_besc(v, Tgas, xds, lj, Tds);
  double cool = 0.0;
  for (int b = 0; b < v.nb; ++b) {
    const double intCMB = dust_intBB(v, b, p.Tcmb);
    cool += (dust_intBB(v, b, Tds[b]) - intCMB - dust_intJ(v, b, Phi)) * be
        * v.xw[b] * xds * dust_a2_(v, b);
  }
  return 4.0 * constants::pi * cool;
}

// ---- dust.evolve_Td: the grains' energy equation
// heat capacity of one grain of bin b at Td [erg/K] (Draine & Li 2001
// eqs. 9-11): graphite (N - 2) k [f2'(T/863 K) + 2 f2'(T/2504 K)],
// silicate (N - 2) k [2 f2'(T/500 K) + f3'(T/1500 K)], with f_n(x) =
// (1/n) int_0^1 y^n dy / (exp(y/x) - 1), as in the paper (the grains are
// destroyed long before any high-temperature limit matters). N: the
// grain's atoms, from the bulk density and a mean atomic mass (graphite
// 12.011 amu; silicate 24.60 amu, olivine MgFeSiO4)
KROME_HD inline double dust_fnprime_(int nn, double x) {
  // d f_n/dx = (1/n) x^n int_0^{1/x} u^(n+1) e^u / (e^u - 1)^2 du
  const double U = std::fmin(1.0 / x, 60.0);
  constexpr double gx[8] = {-0.9602898564975363, -0.7966664774136267,
                            -0.5255324099163290, -0.1834346424956498,
                            0.1834346424956498, 0.5255324099163290,
                            0.7966664774136267, 0.9602898564975363};
  constexpr double gw[8] = {0.1012285362903763, 0.2223810344533745,
                            0.3137066458778873, 0.3626837833783620,
                            0.3626837833783620, 0.3137066458778873,
                            0.2223810344533745, 0.1012285362903763};
  const int np = 16;  // composite Gauss-Legendre panels over [0, U]
  const double h = U / np;
  double s = 0.0;
  for (int k = 0; k < np; ++k)
    for (int q = 0; q < 8; ++q) {
      const double u = h * (k + 0.5 + 0.5 * gx[q]);
      const double em = std::expm1(u);
      s += 0.5 * h * gw[q] * std::pow(u, nn + 1) * (em + 1.0) / (em * em);
    }
  return std::pow(x, nn) * s / nn;
}
KROME_HD inline double dust_heat_capacity(const DustView& v, int b,
                                          double T) {
  const double rho0 = v.isC[b] ? 2.25 : 3.13;          // g cm^-3
  const double mat = (v.isC[b] ? 12.011 : 24.60) * constants::amu;
  const double N = 4.0 / 3.0 * constants::pi * dust_a_(v, b) * dust_a2_(v, b)
      * rho0 / mat;
  const double f = v.isC[b]
      ? dust_fnprime_(2, T / 863.0) + 2.0 * dust_fnprime_(2, T / 2504.0)
      : 2.0 * dust_fnprime_(2, T / 500.0) + dust_fnprime_(3, T / 1500.0);
  // (N - 2) floored at one atom: a grain evaporated to a few atoms
  return std::fmax(N - 2.0, 1.0) * constants::boltzmann_erg * f;
}

// dust.evolve_Td: dTd_b/dt from each grain's energy equation,
//   C_b(Td) dTd/dt = 4 pi a^2 [pre (Tgas - Td) - (intBB(Td) - intBB(Tcmb)
//                     - intJ) besc],
// gas collisions in, net radiation out: Td relaxes to the balance from
// any start (legacy -usedTdust tracked the balance through dTgas/dt
// instead, valid only when started on it)
template <class P>
KROME_HD inline void dust_dTdust(const DustView& v, double Tgas,
                                 const double* n, const P& p,
                                 const double* Tds, const double* Phi,
                                 double* out) {
  double pre, lj, xds;
  dust_pre(v, Tgas, n, pre, lj, xds, p);
  const double be = dust_besc(v, Tgas, xds, lj, Tds);
  for (int b = 0; b < v.nb; ++b) {
    const double net = pre * (Tgas - Tds[b])
        - (dust_intBB(v, b, Tds[b]) - dust_intBB(v, b, p.Tcmb)
           - dust_intJ(v, b, Phi)) * be;
    // a grain evaporated to nothing has no heat capacity: its Td stays
    out[b] = dust_a_(v, b) > 0.0
        ? 4.0 * constants::pi * dust_a2_(v, b) * net
              / dust_heat_capacity(v, b, Tds[b])
        : 0.0;
  }
}

// one bin's radiative-collisional balance at a given escape factor be
// (legacy compute_Tdust: bisection to 1e-8 K, then a Newton polish to
// roundoff on the same piecewise-linear table balance, so Td is a smooth
// function of the gas state for the difference-quotient Jacobian)
KROME_HD inline double dust_solve_bin_(const DustView& v, int b, double Tgas,
                                       double pre, double be, double intCMB,
                                       double iJ) {
  double Td1 = v.TbbMin, Td2 = v.TbbMax - 1.0, Tdmid = 0.5 * (Td1 + Td2);
  for (;;) {
    const double f1 = (dust_intBB(v, b, Td1) - intCMB - iJ) * be
        - pre * (Tgas - Td1);
    Tdmid = 0.5 * (Td1 + Td2);
    const double fmid = (dust_intBB(v, b, Tdmid) - intCMB - iJ) * be
        - pre * (Tgas - Tdmid);
    if (f1 * fmid < 0.0) Td2 = Tdmid;
    else Td1 = Tdmid;
    if (std::fabs(Td1 - Td2) < 1e-8) break;
  }
  for (int it = 0; it < 4; ++it) {
    const double g = (dust_intBB(v, b, Tdmid) - intCMB - iJ) * be
        - pre * (Tgas - Tdmid);
    const double dg = dust_intBB_slope_(v, b, Tdmid) * be + pre;
    if (!(dg > 0.0)) break;
    const double step = g / dg;
    Tdmid = std::fmin(std::fmax(Tdmid - step, v.TbbMin), v.TbbMax);
    if (std::fabs(step) <= 1e-15 * Tdmid) break;
  }
  return Tdmid;
}

// the per-bin dust temperatures Td_out[nb] of the radiative-collisional
// balance at (Tgas, n, field Phi): a function of these arguments only
// (no history). In optically thick gas the escape factor depends on Td
// (the dust Planck opacity), so the balance is a fixed point: iterated
// from Td = Tgas to a relative change below 1e-12 (it contracts by ~3
// per pass; at most 100 passes). Legacy took ONE step of this iteration
// per call, from the previous call's Td, so its Td depended on the call
// history. Thin gas (escape factor 1) needs a single pass. Bins with
// negligible dust take Td = Tgas.
template <class P>
KROME_HD inline void dust_balance(const DustView& v, double Tgas,
                                  const double* n, const P& p,
                                  const double* Phi, double* Td_out) {
  const double ntot = ntot_(n);
  double pre, lj, xd_scale;
  dust_pre(v, Tgas, n, pre, lj, xd_scale, p);
  for (int b = 0; b < v.nb; ++b) Td_out[b] = Tgas;
  for (int pass = 0; pass < 100; ++pass) {
    const double be = dust_besc(v, Tgas, xd_scale, lj, Td_out);
    double change = 0.0;
    for (int b = 0; b < v.nb; ++b) {
      double t = Tgas;
      if (!(v.xw[b] * xd_scale < 1e-30 * ntot))  // else negligible dust
        t = dust_solve_bin_(v, b, Tgas, pre, be, dust_intBB(v, b, p.Tcmb),
                            dust_intJ(v, b, Phi));
      change = std::fmax(change, std::fabs(t - Td_out[b]) / t);
      Td_out[b] = t;
    }
    if (be >= 1.0 && dust_besc(v, Tgas, xd_scale, lj, Td_out) >= 1.0) break;
    if (change < 1e-12) break;
  }
}

// dust.evolve state <-> radii: the solver integrates N_b, the
// key-element atoms locked in bin b per cm^3, N = n_d (4 pi/3) a^3 rho0 /
// m, so that gas + grains is a LINEAR invariant (BDF keeps it to
// roundoff; with the radius as state it held only to the tolerance)
KROME_HD inline double dust_key_mass_(const DustView& v, int b) {
  const int key = v.isC[b] ? DUST_KEY_C : DUST_KEY_SI;
  return key >= 0 ? species_mass_g[key] : 0.0;
}
// atoms per unit a^3 in bin b at the count scale xds
KROME_HD inline double dust_atoms_per_a3(const DustView& v, int b,
                                         double xds) {
  const double m = dust_key_mass_(v, b);
  return m > 0.0
      ? v.xw[b] * xds * 4.0 / 3.0 * constants::pi * (v.isC[b] ? 2.25 : 3.13)
            / m
      : 0.0;
}
KROME_HD inline double dust_radius_of_atoms(const DustView& v, int b,
                                            double N, double xds) {
  const double k = dust_atoms_per_a3(v, b, xds);
  return k > 0.0 ? std::cbrt(std::fmax(N, 0.0) / k) : v.a[b];
}
KROME_HD inline double dust_atoms_of_radius(const DustView& v, int b,
                                            double ab, double xds) {
  return dust_atoms_per_a3(v, b, xds) * ab * ab * ab;
}

// ---- dust.evolve: the grain radii's rates of change
// per bin b of type t (carbon: key element C, rho0 2.25 g cm^-3;
// silicate: Si, 3.13), with m the key element's mass and
// da = (m / rho0)^(1/3) the thickness of one atomic layer:
//   GROWTH  +S(Tg, Td) v n_key m / (4 rho0), S of Leitch-Devlin &
//           Williams 1985 (legacy krome_dust_stick), v = sqrt(8 kB Tg /
//           (pi m)) the key element's thermal speed (legacy: the
//           proton's)
//   SPUTTER -n_p y(Tg) m / (4 rho0) above 1e5 K, y the legacy fit of
//           the Nozawa et al. 2006 yields times the thermal speed, n_p
//           the heavy particles (ntot - n_e): the sputtered mass flux
//           over the grain surface (legacy mixed the initial and current
//           sizes, fixed rho0 = 2.3 and counted the electrons)
//   EVAP    -nu0 da exp(-E0 / Td), nu0 = 1e12 s^-1, E0 = 7.2 eV
//           (carbon, Lenzuni et al. 1995), 4.66 eV (silicate, Nozawa et
//           al. 2006): the Polanyi-Wigner layer removal of Grassi et al.
//           2017 sec. 2.3, size independent (legacy: nu0 m /(4 pi a^2
//           rho0) exp(-E/Td), one atom per grain)
//   KELVIN  (optional) the evaporation times exp(2 sigma Omega / (a kB
//           Td)), the curvature (Kelvin) enhancement of a small grain's
//           vapour pressure: Omega = m / rho0 the volume per key atom,
//           sigma the surface energy, 1400 (carbon) and 800 (silicon)
//           erg cm^-2 (Nozawa et al. 2003 Table 2), a floored at one
//           layer da (the continuum law has no meaning below it)
// all times s(a) = a^2 / (a^2 + da^2): a grain smaller than one layer
// neither grows nor shrinks, so an evaporated bin stays empty (no
// growth from zero-size seeds) and the rates stay smooth. The atoms
// go to / come from the key element: dn_key -= n_d 4 pi a^2 rho0 da/dt
// / m, which conserves its mass exactly.
constexpr double DUST_NU0_EVAP = 1e12;                              // s^-1
constexpr double DUST_E0_C_K = 7.2 / constants::boltzmann_eV;       // K
constexpr double DUST_E0_SI_K = 4.66 / constants::boltzmann_eV;     // K
constexpr double DUST_SIGMA_C = 1400.0;                             // erg cm^-2
constexpr double DUST_SIGMA_SI = 800.0;                             // erg cm^-2
template <class P>
KROME_HD inline void dust_dadt(const DustView& v, double Tgas,
                               const double* n, const P& p,
                               const double* Tds, double* da, double* dn) {
  (void)p;
  double ntot_heavy = ntot_(n);
#ifdef KROME_HAS_E
  ntot_heavy -= n[IDX_E];
#endif
  ntot_heavy = std::fmax(ntot_heavy, 0.0);
  double ysp = 0.0;
  if constexpr (EvolveConfig::dust_sputter) {
    if (Tgas >= 1e5) {
      const double lT = std::log10(Tgas);
      ysp = std::pow(10.0, std::exp(0.39280689440337335 * lT)
                               / (19.746828032397993
                                  - 4.0031865960055839 * lT)
                           + 7.8081665612858187);
    }
  }
  double xds = v.xds_fix;
  if (!(xds > 0.0)) {
    double pre, lj;
    dust_pre(v, Tgas, n, pre, lj, xds, p);
  }
  for (int b = 0; b < v.nb; ++b) {
    const int key = v.isC[b] ? DUST_KEY_C : DUST_KEY_SI;
    da[b] = 0.0;
    if (key < 0) continue;  // a type without its key element: frozen
    const double m = species_mass_g[key];
    const double rho0 = v.isC[b] ? 2.25 : 3.13;
    const double lay = std::cbrt(m / rho0);
    const double ab = dust_a_(v, b);
    const double sw = ab * ab / (ab * ab + lay * lay);
    double r = 0.0;
    if constexpr (EvolveConfig::dust_growth) {
      const double Td = Tds[b];
      const double S = std::fmin(1.9e-2 * Tgas * (1.7e-3 * Td + 0.4)
                                 * std::exp(-7e-3 * Tgas), 1.0);
      const double vth = std::sqrt(8.0 * constants::boltzmann_erg * Tgas
                                   / (constants::pi * m));
      r += S * vth * std::fmax(n[key], 0.0) * m / (4.0 * rho0);
    }
    if constexpr (EvolveConfig::dust_sputter)
      r -= ntot_heavy * ysp * m / (4.0 * rho0);
    if constexpr (EvolveConfig::dust_evap) {
      double ex = -(v.isC[b] ? DUST_E0_C_K : DUST_E0_SI_K) / Tds[b];
      if constexpr (EvolveConfig::dust_kelvin)
        ex += 2.0 * (v.isC[b] ? DUST_SIGMA_C : DUST_SIGMA_SI) * (m / rho0)
            / (std::fmax(ab, lay) * constants::boltzmann_erg * Tds[b]);
      r -= DUST_NU0_EVAP * lay * std::exp(ex);
    }
    da[b] = sw * r;
    // the key element: n_d 4 pi a^2 rho0 / m per unit da/dt
    dn[key] -= v.xw[b] * xds * 4.0 * constants::pi * ab * ab * rho0 / m
        * da[b];
  }
}

// gas-grain collisional exchange [erg cm^-3 s^-1] at the state Td (the
// gas cooling with dust.evolve_Td: energy conserving; it equals the
// radiative form dust_cooling_state only at the balance)
template <class P>
KROME_HD inline double dust_cooling_exchange(const DustView& v, double Tgas,
                                             const double* n, const P& p,
                                             const double* Tds) {
  double pre, lj, xds;
  dust_pre(v, Tgas, n, pre, lj, xds, p);
  double cool = 0.0;
  for (int b = 0; b < v.nb; ++b)
    cool += pre * (Tgas - Tds[b]) * v.xw[b] * xds * dust_a2_(v, b);
  return 4.0 * constants::pi * cool;
}

struct DustBins {
  // configuration (legacy krome_init_dust arguments and -dust option)
  int nd = 10;                    // bins per grain type
  bool useSi = true, useC = true;  // legacy -dust=...,Si,C
  // the MRN size distribution (Mathis, Rumpl & Nordsieck 1977): radii
  // alow..aup [cm], dn/da ~ a^phi (phi: the slope, not a photon flux)
  double alow = 5e-7, aup = 2.5e-5, phi = -3.5;
  int nT = 100000;                // intBB Tbb samples (legacy dust_nT 1e5)
  double TbbMin = 1.0, TbbMax = 1e4;

  // derived (nb = nd * ntypes bins)
  int nb = 0;
  std::vector<double> asize, asize2, xw;  // xdust = xw * rho_gas * d2g
  std::vector<int> isC;                   // grain type per bin (1 carbon)
  std::vector<double> aw;                 // area weights xw * a^2 (for the
                                          // dust-chemistry hook averages)
  std::vector<double> Td;                 // per-bin Tdust (persistent;
                                          // 30 K default as legacy)
  // grain radii [cm]: the bins' current sizes (init(): the reference
  // sizes asize). dust.evolve evolves them: solve() starts from and
  // writes back these, solve_cells uses per-cell rows (CellsOpts::adust).
  // The cross sections (cooling, opacity, the energy equation, the
  // surface-chemistry area weights) and the heat capacity use them; Qabs
  // stays that of the reference sizes (as legacy)
  std::vector<double> a;
  // during a solve: the radii of the state (the ODE slice), and the fixed
  // grain count scale rho_gas d2g of the solve start (dust.evolve: the
  // counts do not follow the gas depleted into the grains); nullptr / < 0
  // outside
  const double* acur = nullptr;
  double xds_fix = -1.0;
  std::vector<double> awcur;  // area weights xw a^2 at the current radii
  std::vector<double> qb;        // Qabs per bin on the data E grid
  std::vector<double> qmid;      // Qabs per bin at photobin midpoints
  std::vector<double> emid;      // photobin midpoints [eV] (interp_qabs)
  // intBB tables [bin][nT]: read-only after init(), shared by the copies
  // (solve_cells gives every cell its own bins for its Td; copying these
  // ~3 nb nT doubles per cell would cost ~48 MB each with the defaults)
  struct IntBB { std::vector<double> ibb, ibb_sig; };
  std::shared_ptr<const IntBB> tab_;
  double TbbMult = 0.0;

  static constexpr double PI_ = constants::pi;
  static constexpr double KB_EV = constants::boltzmann_eV;
  static constexpr double H_EV = constants::planck_eV;
  static constexpr double CL_ = constants::clight;
  static constexpr double SB_ = constants::stefboltz_erg;
  static constexpr double EV2ERG_ = constants::eV_to_erg;

  static double planckBB_(double x, double T) {   // legacy planckBB [eV]
    const double xexp = x / KB_EV / T;
    if (xexp >= 3e2 || x <= 1e-10) return 0.0;
    return 2.0 * x * x * x / (H_EV * H_EV) / (CL_ * CL_)
        / (std::exp(xexp) - 1.0);
  }

  // per-bin energy grid: the Si and C data files sample the same layout
  // but differ in the last printed digit, so each bin keeps its own
  // type's grid (legacy reused whichever type loaded last — a documented
  // tiny deviation)
  int nE = 0;
  std::vector<const double*> eg;  // per-bin energy grid pointer

  void init() {
    const int ntypes = (useSi ? 1 : 0) + (useC ? 1 : 0);
    nb = nd * ntypes;
    asize.assign(nb, 0.0);
    asize2.assign(nb, 0.0);
    xw.assign(nb, 0.0);
    aw.assign(nb, 0.0);
    awcur.assign(nb, 0.0);
    a.assign(nb, 0.0);
    isC.assign(nb, 0);
    Td.assign(nb, 30.0);  // legacy default
    nE = QABS_SI_NE;      // Si and C grids have the same length
    eg.assign(nb, nullptr);
    qb.assign(static_cast<size_t>(nb) * nE, 0.0);
    TbbMult = (nT - 1) / (TbbMax - TbbMin);

    const double p1 = phi + 1.0, p4 = phi + 4.0;
    // key-element fractions (10^solar, Asplund+09), normalized
    double fSi = useSi ? std::pow(10.0, 7.51 - 12.0) : 0.0;
    double fC = useC ? std::pow(10.0, 8.43 - 12.0) : 0.0;
    const double fsum = fSi + fC;
    fSi /= fsum;
    fC /= fsum;
    const double la = std::log10(alow), lb = std::log10(aup);
    int type = 0;
    for (int t = 0; t < 2; ++t) {
      const bool isSi = (t == 0);
      if ((isSi && !useSi) || (!isSi && !useC)) continue;
      const double rho0 = isSi ? 3.13 : 2.25;   // g/cm3, Zhukovska
      const double frac = isSi ? fSi : fC;
      // xdust = myc * (a+^p1 - a^p1)/p1 with
      // myc = frac * rho_gas*d2g / ((4pi/3) rho0 (aup^p4-alow^p4)/p4);
      // xw stores everything but rho_gas*d2g
      const double myc = frac / (4.0 * PI_ / 3.0 * rho0
          * (std::pow(aup, p4) - std::pow(alow, p4)) / p4);
      const int na = isSi ? QABS_SI_NA : QABS_C_NA;
      const double* qa = isSi ? qabs_Si_a : qabs_C_a;
      const double* qq = isSi ? qabs_Si_Q : qabs_C_Q;
      const double* qE = isSi ? qabs_Si_E : qabs_C_E;
      for (int i = 0; i < nd; ++i) {
        const int b = type * nd + i;
        eg[b] = qE;
        isC[b] = isSi ? 0 : 1;
        const double a1 = std::pow(10.0, i * (lb - la) / nd + la);
        const double a2 = std::pow(10.0, (i + 1) * (lb - la) / nd + la);
        asize[b] = std::pow((std::pow(a2, p4) - std::pow(a1, p4))
                            / (std::pow(a2, p1) - std::pow(a1, p1))
                            * (p1 / p4), 1.0 / 3.0);
        asize2[b] = asize[b] * asize[b];
        xw[b] = myc * (std::pow(a2, p1) - std::pow(a1, p1)) / p1;
        aw[b] = xw[b] * asize2[b];
        a[b] = asize[b];
        awcur[b] = aw[b];
        // Qabs interpolated linearly in grain size (legacy)
        int idx = na - 1;
        for (int j = 0; j < na; ++j)
          if (asize[b] < qa[j]) { idx = j; break; }
        if (idx == 0) idx = 1;
        const double wA = (asize[b] - qa[idx - 1]) / (qa[idx] - qa[idx - 1]);
        for (int j = 0; j < nE; ++j)
          qb[static_cast<size_t>(b) * nE + j] =
              (qq[idx * nE + j] - qq[(idx - 1) * nE + j]) * wA
              + qq[(idx - 1) * nE + j];
      }
      ++type;
    }
    // intBB tables: trapezoid over each bin's data E grid at every Tbb
    // (legacy dust_init_intBB); Planck weights shared per grid
    auto tb = std::make_shared<IntBB>();
    std::vector<double>& ibb = tb->ibb;
    std::vector<double>& ibb_sig = tb->ibb_sig;
    ibb.assign(static_cast<size_t>(nb) * nT, 0.0);
    ibb_sig.assign(static_cast<size_t>(nb) * nT, 0.0);
    std::vector<double> wB(nE);
    const double pre = PI_ / H_EV * EV2ERG_;
    for (int i = 0; i < nT; ++i) {
      const double Tbb = i / TbbMult + TbbMin;
      const double* wgrid = nullptr;
      for (int b = 0; b < nb; ++b) {
        if (eg[b] != wgrid) {  // recompute trapezoid weights per grid
          wgrid = eg[b];
          for (int j = 0; j < nE; ++j) wB[j] = 0.0;
          for (int j = 1; j < nE; ++j) {
            const double dE = wgrid[j] - wgrid[j - 1];
            wB[j] += 0.5 * dE * planckBB_(wgrid[j], Tbb);
            wB[j - 1] += 0.5 * dE * planckBB_(wgrid[j - 1], Tbb);
          }
        }
        double s = 0.0;
        const double* q = &qb[static_cast<size_t>(b) * nE];
        for (int j = 0; j < nE; ++j) s += wB[j] * q[j];
        ibb[static_cast<size_t>(b) * nT + i] = s * pre;
        ibb_sig[static_cast<size_t>(b) * nT + i] =
            s * pre / (SB_ * Tbb * Tbb * Tbb * Tbb);
      }
    }
    tab_ = tb;
  }

  // the population as a DustView (its own vectors; acur / xds_fix as set)
  DustView view() const {
    DustView v;
    v.nb = nb;
    v.nT = nT;
    v.TbbMin = TbbMin;
    v.TbbMax = TbbMax;
    v.TbbMult = TbbMult;
    v.xw = xw.data();
    v.isC = isC.data();
    v.a = a.data();
    v.acur = acur;
    v.xds_fix = xds_fix;
    if (tab_) {
      v.ibb = tab_->ibb.data();
      v.ibb_sig = tab_->ibb_sig.data();
    }
    if (!qmid.empty()) {
      v.B = static_cast<int>(emid.size());
      v.qmid = qmid.data();
      v.emid = emid.data();
    }
    return v;
  }

  // the methods: the DustView functions on this population (see there)
  double get_intBB(int b, double T) const { return dust_intBB(view(), b, T); }
  // the radii the physics sees: av (a solve's state slice), or nullptr =
  // the population's own a[]; refreshes the area weights of the hook
  void set_radii(const double* av) {
    acur = av;
    dust_area_weights(view(), awcur.data());
  }

  // dust-chemistry coupling: install this population into the generated
  // DustBinsHook of a Params (p.hk.dust), so table_Tdust and the surface
  // rates (sticking, evaporation) see area-weighted bin averages. Td may
  // be overridden with the ODE-state slice (dust.evolve_Td). The solver
  // installs the hook in its own copy of the parameters; call these
  // directly only when evaluating rates() outside a solve.
  template <class P>
  void hook_install(P& p, const double* Tds = nullptr) const {
    p.hk.dust.nb = nb;
    p.hk.dust.Td = Tds != nullptr ? Tds : Td.data();
    p.hk.dust.aw = awcur.data();
    p.hk.dust.isC = isC.data();
  }
  template <class P>
  static void hook_clear(P& p) {
    p.hk.dust = DustBinsHook();
  }

  template <class P>
  double cooling_state(double Tgas, const double* n, const P& p,
                       const double* Tds, const double* Phi = nullptr)
      const {
    return dust_cooling_state(view(), Tgas, n, p, Tds, Phi);
  }
  template <class P>
  void dTdust(double Tgas, const double* n, const P& p,
              const double* Tds, const double* Phi, double* out) const {
    dust_dTdust(view(), Tgas, n, p, Tds, Phi, out);
  }

  // interpolate Qabs of every bin onto photobin midpoints (legacy
  // interp_qabs); call once after init when radiation is co-evolved
  template <class Rad>
  void interp_qabs(const Rad& rad) {
    const int B = rad.nbins();
    qmid.assign(static_cast<size_t>(nb) * B, 0.0);
    emid.resize(B);
    for (int i = 0; i < B; ++i)
      emid[i] = 0.5 * (rad.Elo[i] + rad.Ehi[i]);
    for (int b = 0; b < nb; ++b) {
      const double* q = &qb[static_cast<size_t>(b) * nE];
      const double* Egrid = eg[b];
      for (int i = 0; i < B; ++i) {
        const double E = emid[i];
        double v = q[nE - 1];
        for (int j = 1; j < nE; ++j)
          if (E < Egrid[j]) {
            v = (E - Egrid[j - 1]) / (Egrid[j] - Egrid[j - 1])
                * (q[j] - q[j - 1]) + q[j - 1];
            break;
          }
        qmid[static_cast<size_t>(b) * B + i] = v;
      }
    }
  }
  template <class P>
  void balance(double Tgas, const double* n, const P& p, const double* Phi,
               double* Td_out) const {
    dust_balance(view(), Tgas, n, p, Phi, Td_out);
  }
  // the balance Td (stored in Td[]) and the total dust cooling [erg cm^-3
  // s^-1] (legacy compute_Tdust; negative = the dust heats the gas)
  template <class P>
  double cooling(double Tgas, const double* n, const P& p,
                 const double* Phi = nullptr) {
    balance(Tgas, n, p, Phi, Td.data());
    return cooling_state(Tgas, n, p, Td.data(), Phi);
  }
  double heat_capacity(int b, double T) const {
    return dust_heat_capacity(view(), b, T);
  }
  double radius_of_atoms(int b, double N, double xds) const {
    return dust_radius_of_atoms(view(), b, N, xds);
  }
  double atoms_of_radius(int b, double ab, double xds) const {
    return dust_atoms_of_radius(view(), b, ab, xds);
  }
  template <class P>
  void dadt(double Tgas, const double* n, const P& p, const double* Tds,
            double* da, double* dn) const {
    dust_dadt(view(), Tgas, n, p, Tds, da, dn);
  }
};

} // namespace krome
