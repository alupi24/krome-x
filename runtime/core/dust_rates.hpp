// KROME-X dust module: recombination on dust grains (Weingartner & Draine
// 2001 fits) and freezeout/evaporation rates, ported verbatim from legacy
// KROME (src/krome_subs.f90, src/krome_grfuncs.f90).
//
// Include the generated net_<name>.hpp BEFORE this header. Each function is
// compiled only when the network uses it (KROME_EXTERN_* defines), so
// networks without e.g. an electron species do not instantiate the
// recombination fits.
//
// Conventions: p.G0 is the local Habing field (legacy GHabing), p.d2g the
// dust-to-gas mass ratio (0.01 = Milky Way; the WD01 fits are normalized to
// MW dust, hence the d2g/0.01 scaling), p.Tdust the dust temperature.
#pragma once
#include <cmath>

#define KROME_DUST_RATES_HPP_INCLUDED 1

namespace krome {

// sqrt(8 * boltzmann_erg / pi), legacy pre_kvgas_sqrt (a literal: no
// constexpr sqrt in C++17; tests/test_constants.py checks it against
// kromex/constants.py)
constexpr double PRE_KVGAS_SQRT = 1.87504501503247e-08;

#if defined(KROME_EXTERN_H_recombination_on_dust) || \
    defined(KROME_EXTERN_He_recombination_on_dust) || \
    defined(KROME_EXTERN_C_recombination_on_dust) || \
    defined(KROME_EXTERN_Si_recombination_on_dust) || \
    defined(KROME_EXTERN_O_recombination_on_dust) || \
    defined(KROME_EXTERN_Fe_recombination_on_dust)
// Weingartner & Draine 2001 fit, k = pre / (1 + c1 psi^p1 (1 + c2 Tgas^q
// psi^(p2a + p2b ln Tgas))), psi = G0 sqrt(Tgas)/n_e.
// Limits (p1 + p2a > 0 for every species): psi -> 0 gives k -> pre (grains
// maximally negatively charged, fastest recombination); psi -> inf (n_e ->
// 0) gives k -> 0. NOTE, deliberate deviation from legacy KROME: legacy
// returned 0 for G0 <= 0, switching the channel off discontinuously; here
// the psi -> 0 limit of the fit is returned instead.
KROME_HD static inline double wd01_fit_(double pre, double c1, double p1, double c2,
                               double q, double p2a, double p2b,
                               double Tgas, const double* n,
                               const Params& p) {
  if (n[IDX_E] < 1e-20) return 0.0;  // psi -> inf limit
  const double scale = pre * (p.d2g / 1e-2);
  const double psi = p.G0 * std::sqrt(Tgas) / n[IDX_E];
  if (psi <= 0.0) return scale;      // psi -> 0 limit
  return scale / (1.0 + c1 * std::pow(psi, p1)
                  * (1.0 + c2 * std::pow(Tgas, q)
                     * std::pow(psi, p2a + p2b * std::log(Tgas))));
}
#endif

#ifdef KROME_EXTERN_H_recombination_on_dust
KROME_HD inline double H_recombination_on_dust(double Tgas, const double* n,
                                      const Params& p) {
  return wd01_fit_(1.225e-13, 8.074e-6, 1.378, 5.087e2, 0.01586,
                   -0.4723, -1.102e-5, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_He_recombination_on_dust
KROME_HD inline double He_recombination_on_dust(double Tgas, const double* n,
                                       const Params& p) {
  return wd01_fit_(5.572e-14, 3.185e-7, 1.512, 5.115e3, 3.903e-7,
                   -0.4956, -5.494e-7, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_C_recombination_on_dust
KROME_HD inline double C_recombination_on_dust(double Tgas, const double* n,
                                      const Params& p) {
  return wd01_fit_(4.558e-13, 6.089e-3, 1.128, 4.331e2, 0.04845,
                   -0.8120, -1.333e-4, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_Si_recombination_on_dust
KROME_HD inline double Si_recombination_on_dust(double Tgas, const double* n,
                                       const Params& p) {
  return wd01_fit_(2.166e-14, 5.678e-8, 1.874, 4.375e4, 1.635e-6,
                   -0.8964, -7.538e-5, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_O_recombination_on_dust
KROME_HD inline double O_recombination_on_dust(double Tgas, const double* n,
                                      const Params& p) {
  // legacy: 0.25 * k_H (Weingartner & Draine 2001 scaling)
  return 0.25 * wd01_fit_(1.225e-13, 8.074e-6, 1.378, 5.087e2, 0.01586,
                          -0.4723, -1.102e-5, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_Fe_recombination_on_dust
KROME_HD inline double Fe_recombination_on_dust(double Tgas, const double* n,
                                       const Params& p) {
  // Fe+ (the KROME branch of the galaxy_ism networks; not in legacy
  // KROME's krome_subs.f90): the same fit with the Fe coefficients
  return wd01_fit_(1.701e-14, 9.554e-8, 1.851, 5.763e4, 4.116e-8,
                   -0.9456, -2.198e-5, Tgas, n, p);
}
#endif

#ifdef KROME_EXTERN_krate_stickSi
// freezeout on an MRN silicate population (legacy krate_stickSi defaults:
// amin 5e-7 cm, amax 2.5e-5 cm, p=-3.5, rho0=3 g/cm3) with the cell's d2g
// (legacy hard-coded d2g = 1e-2 here, whatever the dust-to-gas ratio)
KROME_HD inline double krate_stickSi(int sp, double Tgas, const double* n,
                            const Params& p) {
  const double amin = 5e-7, amax = 2.5e-5, pexp = -3.5, rho0 = 3.0;
  const double p3 = pexp + 3.0, p4 = pexp + 4.0;
  // rho_ is the weighted helper: correct in both representations
  const double rhod = rho_(n) * p.d2g;
  const double Tg100 = Tgas * 1e-2;
  const double vpre = PRE_KVGAS_SQRT * std::sqrt(Tgas)
      / std::sqrt(species_mass_g[sp]);
  // dust-chemistry coupling: with a binned population installed
  // (DustBinsHook), the grain area moment AND the sticking factor come
  // from the bins as area-weighted averages over the per-bin Td
  if (p.hk.dust.nb > 0) {
    const DustBinsHook& h = p.hk.dust;
    const double rhod_b = rho_(n) * p.d2g;
    double stw = 0.0;
    for (int b = 0; b < h.nb; ++b) {
      const double Td100 = h.Td[b] * 1e-2;
      const double st = 1.0 / (1.0 + 0.4 * std::sqrt(Tg100 + Td100)
                               + 0.2 * Tg100 + 0.08 * Tg100 * Tg100);
      stw += h.aw[b] * st;
    }
    // Sum_b n_d,b pi a_b^2 stick_b = pi (rho_gas d2g) Sum_b aw_b stick_b
    return vpre * constants::pi * rhod_b * stw;
  }
  // effective MRN population: live Tdust (table when dust.tables is
  // generated, else the p.Tdust parameter)
  const double Td100 = table_Tdust(Tgas, n, p) * 1e-2;
  const double stick = 1.0 / (1.0 + 0.4 * std::sqrt(Tg100 + Td100)
                              + 0.2 * Tg100 + 0.08 * Tg100 * Tg100);
  return vpre * rhod / (4.0 / 3.0 * rho0) * p4 / p3
         * (std::pow(amax, p3) - std::pow(amin, p3))
         / (std::pow(amax, p4) - std::pow(amin, p4)) * stick;
}
#endif

#ifdef KROME_EXTERN_krate_evaporation
// thermal evaporation, nu0 * exp(-Ebind_bare / Tdust)
KROME_HD inline double krate_evaporation(int sp, double Tgas, const double* n,
                                const Params& p) {
  const double Eb = species_Ebind_bare[sp];
  // dust-chemistry coupling: with bins installed, the area-weighted
  // average of the per-bin RATE (not the rate at a mean Td — the
  // exponential makes hot small grains dominate desorption)
  if (p.hk.dust.nb > 0) {
    const DustBinsHook& h = p.hk.dust;
    double s = 0.0, w = 0.0;
    for (int b = 0; b < h.nb; ++b) {
      s += h.aw[b] * std::exp(-Eb / h.Td[b]);
      w += h.aw[b];
    }
    return 1e12 * s / w;
  }
  // live dust temperature (see krate_stickSi)
  return 1e12 * std::exp(-Eb / table_Tdust(Tgas, n, p));
}
#endif

// gas-dust collisional thermal exchange (Hollenbach & McKee 1979; legacy
// cooling_dust / dustCool with fact = 0.5 mean-collider factor) for the
// effective MRN population: 2 kB fact v_gas (Tgas - Tdust) ntot Sum n_d a^2,
// with Sum n_d a^2 from the total dust mass rho_(n) * d2g and the same MRN
// defaults as krate_stickSi (amin 5e-7, amax 2.5e-5 cm, p = -3.5,
// rho0 = 3 g/cm3). Positive = gas cools (Tgas > Tdust), negative = dust
// heats the gas, exactly as legacy. Tdust is resolved by the caller
// (table_Tdust).
template <class P>
KROME_HD inline double cooling_dustExchange(double Tgas, double Tdust,
                                   const double* n, const P& p) {
  const double amin = 5e-7, amax = 2.5e-5, pexp = -3.5, rho0 = 3.0;
  const double p3 = pexp + 3.0, p4 = pexp + 4.0;
  const double rhod = rho_(n) * p.d2g;
  // Sum n_d a^2 = rho_d * 3/(4 pi rho0) * (p4/p3) * da^p3 / da^p4
  const double na2 = rhod * 3.0 / (4.0 * constants::pi * rho0) * (p4 / p3)
      * (std::pow(amax, p3) - std::pow(amin, p3))
      / (std::pow(amax, p4) - std::pow(amin, p4));
  // vgas = sqrt(8 kB Tgas / (pi m_p)), legacy sqrt(kvgas_erg * Tgas)
  const double vgas = PRE_KVGAS_SQRT * std::sqrt(Tgas / constants::p_mass);
  return 2.0 * constants::boltzmann_erg * 0.5 * vgas * (Tgas - Tdust)
      * ntot_(n) * na2;
}

#if defined(KROME_EXTERN_krate_H2dust)
// H2 formation efficiency on grains (legacy H2_eps_Si / H2_eps_C, Cazaux &
// Spaans 2009 style), and the Hollenbach & McKee 1979 sticking of H
KROME_HD inline double h2dust_eps_Si_(double Tg, double Td) {
  const double Ec = 1.5e4, Es = -1e3, Ep = 7e2, apc = 1.7e-10;
  const double sq = 1.0 + std::sqrt((Ec - Es) / (Ep - Es));
  const double func = 2.0 * std::exp(-(Ep - Es) / (Ep + Tg)) / (sq * sq);
  const double e = 1.0 / (1.0 + 16.0 * Td / (Ec - Es) * std::exp(-Ep / Td)
                          * std::exp(4e9 * apc * std::sqrt(Ep - Es))) + func;
  return std::fmin(1.0, std::fmax(e, 0.0));
}
KROME_HD inline double h2dust_eps_C_(double Tg, double Td) {
  const double Ep = 8e2, Ec = 7e3, Es = 2e2;
  const double sq = 1.0 + std::sqrt((Ec - Es) / (Ep - Es));
  const double TH = 4.0 / (sq * sq) * std::exp(-(Ep - Es) / (Ep + Tg));
  const double e = (1.0 - TH) / (1.0 + 0.25 * sq * sq * std::exp(-Es / Td));
  return std::fmin(1.0, std::fmax(e, 0.0));
}
KROME_HD inline double h2dust_stick_(double Tg, double Td) {
  return 1.0 / (1.0 + 0.4 * std::sqrt((Tg + Td) * 1e-2) + 0.2 * Tg * 1e-2
                + 0.08 * (Tg * 1e-2) * (Tg * 1e-2));
}

// H2 formation on dust (legacy -dustOptions=H2, krome_H2_dust): the rate
// of H + {H} -> H2, nH2dust / n_H = 0.5 v_H sum_b n_d,b pi a_b^2 eps S. With
// a binned population installed (DustBinsHook) the bins, their Td and
// grain types; else the effective MRN silicate population of krate_stickSi
// (with the real d2g) at the live Tdust. Legacy used the first type's bin
// sizes for every type; here each bin has its own.
KROME_HD inline double krate_H2dust(double Tgas, const double* n, const Params& p) {
  const double vH = PRE_KVGAS_SQRT * std::sqrt(Tgas / constants::p_mass);
  if (p.hk.dust.nb > 0) {
    const DustBinsHook& h = p.hk.dust;
    double s = 0.0;
    for (int b = 0; b < h.nb; ++b) {
      const double Td = h.Td[b];
      const double eps = h.isC[b] ? h2dust_eps_C_(Tgas, Td)
                                  : h2dust_eps_Si_(Tgas, Td);
      s += h.aw[b] * eps * h2dust_stick_(Tgas, Td);
    }
    // n_d,b a_b^2 = aw_b rho_gas d2g
    return 0.5 * vH * constants::pi * rho_(n) * p.d2g * s;
  }
  const double amin = 5e-7, amax = 2.5e-5, pexp = -3.5, rho0 = 3.0;
  const double p3 = pexp + 3.0, p4 = pexp + 4.0;
  // sum n_d a^2 of the MRN population holding rho_gas d2g
  const double na2 = rho_(n) * p.d2g * 3.0 / (4.0 * constants::pi * rho0)
      * (p4 / p3) * (std::pow(amax, p3) - std::pow(amin, p3))
      / (std::pow(amax, p4) - std::pow(amin, p4));
  const double Td = table_Tdust(Tgas, n, p);
  return 0.5 * vH * constants::pi * na2 * h2dust_eps_Si_(Tgas, Td)
      * h2dust_stick_(Tgas, Td);
}
#endif

#ifdef KROME_EXTERN_krate_H2dust_Jura
// H2 formation on dust, constant form (legacy -useDustH2const, Jura 1975 +
// Gnedin): nH2dust = 3.5e-17 n_H n_Hnuclei (d2g/0.01) clump, the rate of
// H + {H} -> H2. Legacy's dust2gas_ratio is relative to the MW value; KROME-X's
// d2g is absolute (0.01 = MW)
KROME_HD inline double krate_H2dust_Jura(double Tgas, const double* n,
                                const Params& p) {
  (void)Tgas;
  return 3.5e-17 * Hnuclei_(n) * (p.d2g / 1e-2) * p.clump;
}
#endif

#ifdef KROME_EXTERN_krate_nonthermal_evaporation
// photo + cosmic-ray desorption (Hollenbach+ 2009), legacy constants:
// F_not = 1e8 desorbing photons cm^-2 s^-1, site separation 3e-8 cm
KROME_HD inline double krate_nonthermal_evaporation(int sp, double Gnot, double Av,
                                           double crflux, double yield,
                                           double Tgas, const double* n,
                                           const Params& p) {
  (void)Tgas; (void)n; (void)p;
  const double crnot = 1.3e-17, Fnot = 1e8, ap2 = 3e-8 * 3e-8;
  const double f70 = 3.16e-19 * crflux / crnot;
  const double kevap70 = 1e12 * std::exp(-species_Ebind_bare[sp] / 70.0);
  return Gnot * Fnot * ap2 * yield * std::exp(-1.8 * Av) + f70 * kevap70;
}
#endif

#ifdef KROME_EXTERN_sqimass
// 1 / sqrt(m) of a species [g^-1/2] (legacy get_imass_sqrt)
KROME_HD inline double sqimass(int sp, double Tgas, const double* n,
                      const Params& p) {
  (void)Tgas; (void)n; (void)p;
  return 1.0 / std::sqrt(species_mass_g[sp]);
}
#endif

} // namespace krome
