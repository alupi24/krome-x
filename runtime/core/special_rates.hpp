// KROME-X miscellaneous rate functions ported verbatim from legacy
// KROME (src/krome_subs.f90). Include the generated net_<name>.hpp
// BEFORE this header (uses Hnuclei_).
#pragma once
#include <cmath>

#define KROME_SPECIAL_RATES_HPP_INCLUDED 1

namespace krome {

#ifdef KROME_EXTERN_dissH2_Martin96
// H2 collisional dissociation: collision-induced dissociation +
// dissociative tunneling, density-dependent interpolation between the
// low- and high-density limits (Martin, Schwarz & Mandy 1996), legacy
// dissH2_Martin96 with its exact fit coefficients. n_H is the total
// hydrogen nuclei density.
KROME_HD inline double dissH2_Martin96(double Tgas, const double* n,
                              const Params& p) {
  (void)p;
  static constexpr double m[2][21] = {
      // collision-induced dissociation
      {-178.4239, -68.42243, 43.20243, -4.633167, 69.70086, 40870.38,
       -23705.70, 128.8953, -53.91334, 5.315517, -19.73427, 16780.95,
       -25786.11, 14.82123, -4.890915, 0.4749030, -133.8283, -1.164408,
       0.8227443, 0.5864073, -2.056313},
      // dissociative tunneling
      {-142.7664, 42.70741, -2.027365, -0.2582097, 21.36094, 27535.31,
       -21467.79, 60.34928, -27.43096, 2.676150, -11.28215, 14254.55,
       -23125.20, 9.305564, -2.464009, 0.1985955, 743.0600, -1.174242,
       0.7502286, 0.2358848, 2.937507}};
  const double n_H = Hnuclei_(n);
  const double logT = std::log10(Tgas);
  const double invT = 1.0 / Tgas;
  const double logT2 = logT * logT;
  const double logT3 = logT2 * logT;
  double k_CID = 0.0;
  for (int i = 0; i < 2; ++i) {
    const double* c = m[i];
    const double logk_h1 = c[0] + c[1] * logT + c[2] * logT2
        + c[3] * logT3 + c[4] * std::log10(1.0 + c[5] * invT);
    const double logk_h2 = c[6] * invT;
    const double logk_l1 = c[7] + c[8] * logT + c[9] * logT2
        + c[10] * std::log10(1.0 + c[11] * invT);
    const double logk_l2 = c[12] * invT;
    const double logn_c1 = c[13] + c[14] * logT + c[15] * logT2
        + c[16] * invT;
    const double logn_c2 = c[17] + logn_c1;
    const double pp = c[18] + c[19] * std::exp(-Tgas / 1.850e3)
        + c[20] * std::exp(-Tgas / 4.40e2);
    const double n_c1 = std::pow(10.0, logn_c1);
    const double n_c2 = std::pow(10.0, logn_c2);
    const double logk_CID = logk_h1
        - (logk_h1 - logk_l1) / (1.0 + std::pow(n_H / n_c1, pp))
        + logk_h2
        - (logk_h2 - logk_l2) / (1.0 + std::pow(n_H / n_c2, pp));
    k_CID += std::pow(10.0, logk_CID);
  }
  return k_CID;
}
#endif

#ifdef KROME_EXTERN_elec_recomb_ST93
// H+ recombination in the early universe including the effect of the
// radiation field (Sasaki & Takahara 1993), legacy elec_recomb_ST93. The
// cosmology comes from the parameters, where legacy hard-coded h = 0.704
// and Omega0 = 1: h = hubble0 / (100 km s^-1 Mpc^-1), Omega0 = omega0,
// z = zredshift.
KROME_HD inline double elec_recomb_ST93(double nabund, double nelec, double ntot,
                               double nucleiH, double Trad, double Tgas,
                               const double* n, const Params& p) {
  (void)Tgas; (void)n;
  const double hubble = p.hubble0 / 3.2407792894443648e-18;
  const double Omega0 = p.omega0;
  const double z = p.zredshift;
  const double al = 8.227;
  const double ak = 22.06 / (hubble * (1.0 + z)
                             * std::sqrt(1.0 + Omega0 * z));
  const double rc2 = 8.76e-11 * std::pow(1.0 + z, -0.58);
  const double r2c = std::pow(1.80e10 * Trad, 1.5)
      * std::exp(-3.9472e4 / Trad) * rc2;
  const double a0 = nucleiH * rc2;
  const double b0 = ak * al * nucleiH;
  const double c0 = ak * rc2 * nucleiH * nucleiH;
  const double d0 = r2c * std::exp(-1.18416e5 / Trad);
  const double e0 = ak * r2c * nucleiH;
  const double a1 = -d0 * (1.0 + b0);
  const double b1 = d0 * (1.0 + 2.0 * b0);
  const double c1 = a0 + b0 * (a0 - d0);
  const double d1 = -a0 * b0;
  const double e1 = a0 * c0;
  const double f1 = 1.0 + b0 + e0;
  const double g1 = -(b0 + e0);
  const double h1 = c0;
  const double r = nabund / ntot;
  const double v = ntot * (a1 + b1 * r + c1 * r * r + d1 * r * r * r
                           + e1 * r * r * r * r)
      / (f1 + g1 * r + h1 * r * r);
  return v / (nabund * nelec);
}
#endif

#if defined(KROME_EXTERN_recfast_H_rec) || defined(KROME_EXTERN_recfast_H_ion) || \
    defined(KROME_EXTERN_recfast_He_rec) || defined(KROME_EXTERN_recfast_He_ion)
// Cosmological recombination of H and He I as RECFAST 1.5.2 (Seager,
// Sasselov & Scott 1999, 2000; Wong, Moss & Scott 2008; Rubino-Martin et
// al. 2010; recfast.for, subroutine ION, with Hswitch = 1 and Heswitch =
// 6), written as rate coefficients: each net RECFAST rate
// C (n_e n_X+ alpha - n_X beta e^-E/kT) becomes a recombination
// X+ + e- -> X with k = C alpha and a photoionization X -> X+ + e- with
// k = C beta e^-E/kT (the CMB in detailed balance at the matter
// temperature, as RECFAST). The escape factors C need the Hubble rate
// H(z) = hubble0 (1 + z) sqrt(1 + omega0 z), the matter-dominated form of
// the EXPANSION cooling and of the earlyUniverse driver (RECFAST's own
// H(z) has radiation, curvature and Lambda), z = zredshift. Constants:
// RECFAST's atomic data (wavenumbers, two-photon and Einstein rates, HeI
// cross sections, the fits and fudges) converted to cgs; the fundamental
// constants are krome::constants.
struct RecfastRates_ {
  double aH, bH, CH;        // H: alpha_B, beta e^(-h nu_a / kT), C factor
  double aHe, bHe, CHe;     // HeI singlet
  double aTr, bTr, CTr;     // HeI triplets
};
// RECFAST 1.5.2 data [cgs]
constexpr double RF_L_H_ion = 1.096787737e5;      // cm^-1
constexpr double RF_L_H_alpha = 8.225916453e4;    // cm^-1 (2p, 2s average)
constexpr double RF_L_He1_ion = 1.98310772e5;     // Drake (1993)
constexpr double RF_L_He_2s = 1.66277434e5;
constexpr double RF_L_He_2p = 1.71134891e5;
constexpr double RF_L_He_2Pt = 1.690871466e5;     // Drake & Morton (2007)
constexpr double RF_L_He_2St = 1.5985597526e5;
constexpr double RF_L_He2St_ion = 3.8454693845e4;
constexpr double RF_Lambda = 8.2245809;           // H 2s-1s two-photon [s^-1]
constexpr double RF_Lambda_He = 51.3;             // HeI 2s-1s (Dalgarno)
constexpr double RF_A2P_s = 1.798287e9;           // Morton, Wu & Drake (2006)
constexpr double RF_A2P_t = 177.58;               // Lach & Pachucki (2001)
constexpr double RF_sigma_He_2Ps = 1.436289e-18;  // cm^2, Hummer & Storey (1998)
constexpr double RF_sigma_He_2Pt = 1.484872e-18;
constexpr double RF_not4 = 3.9715;                // m_He / m_H
constexpr double RF_fu = 1.125;                   // H fudge (Hswitch = 1)
constexpr double RF_b_He = 0.86;                  // HeI fudge
KROME_HD inline RecfastRates_ recfast_rates_(double nH0, double nHe0,
                                             double Tgas, const Params& p) {
  RecfastRates_ r;
  const double T = Tgas;
  const double z = p.zredshift;
  const double Hz = p.hubble0 * (1.0 + z) * std::sqrt(1.0 + p.omega0 * z);
  const double kB = constants::boltzmann_erg, hck = constants::hc_k;
  const double CR = 2.0 * constants::pi * constants::e_mass * kB
      / (constants::planck_erg * constants::planck_erg);  // cm^-2 K^-1
  const double crT15 = std::pow(CR * T, 1.5);
  // --- H: Pequignot, Petitjean & Boisson case-B alpha, the Peebles C with
  // the fudge fu and the double-Gaussian K correction (Hswitch = 1)
  const double t4 = T * 1e-4;
  r.aH = 1e-13 * 4.309 * std::pow(t4, -0.6166) / (1.0 + 0.6703 *
                                                  std::pow(t4, 0.53));
  const double RupH = r.aH * crT15
      * std::exp(-hck * (RF_L_H_ion - RF_L_H_alpha) / T);
  r.bH = RupH * std::exp(-hck * RF_L_H_alpha / T);
  const double lz = std::log(1.0 + z);
  const double K = 1.0 / (RF_L_H_alpha * RF_L_H_alpha * RF_L_H_alpha)
      / (8.0 * constants::pi * Hz)
      * (1.0 - 0.14 * std::exp(-((lz - 7.28) / 0.18) * ((lz - 7.28) / 0.18))
         + 0.079 * std::exp(-((lz - 6.73) / 0.33) * ((lz - 6.73) / 0.33)));
  const double KL = K * RF_Lambda * nH0;
  r.CH = RF_fu * (1.0 + KL) / (1.0 + KL + RF_fu * K * RupH * nH0);
  // --- HeI singlet: Verner & Ferland type alpha (the SSS fit), the
  // Sobolev escape probability of 2^1P-1^1S with the H continuum opacity
  // (Kholupenko, Ivanchik & Varshalovich 2007, b_He fudge)
  const double sq0 = std::sqrt(T / std::pow(10.0, 0.477121));
  const double sq1 = std::sqrt(T / std::pow(10.0, 5.114));
  r.aHe = std::pow(10.0, -10.744) / (sq0 * std::pow(1.0 + sq0, 1.0 - 0.711))
      / std::pow(1.0 + sq1, 1.0 + 0.711);
  const double RupHe = 4.0 * r.aHe * crT15
      * std::exp(-hck * (RF_L_He1_ion - RF_L_He_2s) / T);
  r.bHe = RupHe * std::exp(-hck * RF_L_He_2s / T);
  const double Bf = hck * (RF_L_He_2p - RF_L_He_2s) / T;
  const double HeBoltz = std::exp(std::fmin(Bf, 680.0));
  const double tau_s = RF_A2P_s / (RF_L_He_2p * RF_L_He_2p * RF_L_He_2p)
      / (8.0 * constants::pi) * 3.0 * nHe0 / Hz;
  const double p_s = tau_s > 1e-8 ? -std::expm1(-tau_s) / tau_s : 1.0;
  const double mH = constants::p_mass + constants::e_mass;
  const double c = constants::clight;
  // the H continuum opacity: gamma ~ n_He0 / n_H0 (none without H0)
  double AHcon = 0.0;
  if (nH0 > 0.0) {
    const double Dop = c * RF_L_He_2p
        * std::sqrt(2.0 * kB * T / (mH * RF_not4 * c * c));
    const double g = 3.0 * RF_A2P_s * nHe0 / nH0
        / (std::sqrt(constants::pi) * RF_sigma_He_2Ps * 8.0 * constants::pi
           * Dop * RF_L_He_2p * RF_L_He_2p);
    AHcon = RF_A2P_s / (1.0 + 0.36 * std::pow(g, RF_b_He));
  }
  // K_He n_He0 = 1 / (3 (A p + A_Hcon)): finite as n_He0 -> 0
  const double KHn = 1.0 / (3.0 * (RF_A2P_s * p_s + AHcon));
  r.CHe = (1.0 + KHn * RF_Lambda_He * HeBoltz)
      / (1.0 + KHn * (RF_Lambda_He + RupHe) * HeBoltz);
  // --- HeI triplets (Heswitch >= 3), with the H continuum opacity
  r.aTr = std::pow(10.0, -10.306) / (sq0 * std::pow(1.0 + sq0, 1.0 - 0.761))
      / std::pow(1.0 + sq1, 1.0 + 0.761);
  const double RupTr = r.aTr * std::exp(-hck * RF_L_He2St_ion / T) * crT15
      * (4.0 / 3.0);
  r.bTr = 3.0 * RupTr * std::exp(-hck * RF_L_He_2St / T);
  const double tau_t = RF_A2P_t * nHe0 * 3.0
      / (8.0 * constants::pi * Hz * RF_L_He_2Pt * RF_L_He_2Pt * RF_L_He_2Pt);
  const double p_t = tau_t > 1e-8 ? -std::expm1(-tau_t) / tau_t : 1.0;
  double AHt = 0.0;
  if (nH0 > 0.0) {
    const double Dop = c * RF_L_He_2Pt
        * std::sqrt(2.0 * kB * T / (mH * RF_not4 * c * c));
    const double g = 3.0 * RF_A2P_t * nHe0 / nH0
        / (std::sqrt(constants::pi) * RF_sigma_He_2Pt * 8.0 * constants::pi
           * Dop * RF_L_He_2Pt * RF_L_He_2Pt);
    AHt = RF_A2P_t / (1.0 + 0.66 * std::pow(g, 0.9)) / 3.0;
  }
  // C = cf / (Rup + cf), cf = (A p + A_Hcon) e^(-hc(2^3P - 2^3S)/kT), as
  // 1 / (1 + Rup / cf): the two Boltzmann factors combined (cf alone
  // underflows below ~18 K, and 0 / 0)
  const double ru = r.aTr * crT15 * (4.0 / 3.0) / (RF_A2P_t * p_t + AHt)
      * std::exp(-hck * (RF_L_He2St_ion - (RF_L_He_2Pt - RF_L_He_2St)) / T);
  r.CTr = 1.0 / (1.0 + ru);
  return r;
}
#endif
#ifdef KROME_EXTERN_recfast_H_rec
// H+ + e- -> H, RECFAST: C alpha_B [cm^3 s^-1]; nH0 = n(H)
KROME_HD inline double recfast_H_rec(double nH0, double Tgas,
                                     const double* n, const Params& p) {
  (void)n;
  const RecfastRates_ r = recfast_rates_(nH0, 0.0, Tgas, p);
  return r.CH * r.aH;
}
#endif
#ifdef KROME_EXTERN_recfast_H_ion
// H -> H+ + e- by the CMB, RECFAST: C beta e^(-h nu_a / kT) [s^-1]
KROME_HD inline double recfast_H_ion(double nH0, double Tgas,
                                     const double* n, const Params& p) {
  (void)n;
  const RecfastRates_ r = recfast_rates_(nH0, 0.0, Tgas, p);
  return r.CH * r.bH;
}
#endif
#ifdef KROME_EXTERN_recfast_He_rec
// He+ + e- -> He, RECFAST: singlet C alpha + triplet C alpha [cm^3 s^-1];
// nHe0 = n(He), nH0 = n(H)
KROME_HD inline double recfast_He_rec(double nHe0, double nH0, double Tgas,
                                      const double* n, const Params& p) {
  (void)n;
  const RecfastRates_ r = recfast_rates_(nH0, nHe0, Tgas, p);
  return r.CHe * r.aHe + r.CTr * r.aTr;
}
#endif
#ifdef KROME_EXTERN_recfast_He_ion
// He -> He+ + e- by the CMB, RECFAST: singlet + triplet [s^-1]
KROME_HD inline double recfast_He_ion(double nHe0, double nH0, double Tgas,
                                      const double* n, const Params& p) {
  (void)n;
  const RecfastRates_ r = recfast_rates_(nH0, nHe0, Tgas, p);
  return r.CHe * r.bHe + r.CTr * r.bTr;
}
#endif

} // namespace krome
