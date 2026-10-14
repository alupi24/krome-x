// KROME-X standard cooling module, ported verbatim from legacy KROME
// (src/krome_cooling.f90):
//   cooling_Atomic  Cen 1992: collisional ionization and excitation
//                   [erg cm^-3 s^-1]; the recombination cooling of
//                   thermal.cooling = ["ATOMIC"] is generated from the network's
//                   rates (cooling_recomb_)
//   cooling_H2      Glover & Abel 2008 low-density colliders + Hollenbach &
//                   McKee 1979 / Glover 2015 high-density limit
// Terms whose species are absent from the network compile away via the
// generated KROME_HAS_* flags (the analog of legacy #IFKROME_hasX).
// Include the generated net_<name>.hpp BEFORE this header, and this header
// BEFORE solve_cell.hpp to enable the standard-term hooks there.
#pragma once
#include <cmath>

#define KROME_COOLING_HPP_INCLUDED 1

namespace krome {

KROME_HD inline double cool_sigmoid_(double x, double x0, double s) {
  return 10.0 / (10.0 + std::exp(-s * (x - x0)));  // legacy 'sigmoid'
}

KROME_HD inline double wCool_(double logT, double logTmin, double logTmax) {
  const double x = (logT - logTmin) / (logTmax - logTmin);
  const double w = std::pow(10.0, 200.0 * (cool_sigmoid_(x, -0.2, 50.0)
                            * cool_sigmoid_(-x, -1.2, 50.0) - 1.0));
  return (w < 1e-199) ? 0.0 : w;
}

#ifdef KROME_HAS_E
KROME_HD inline double cooling_Atomic(double Tgas, const double* n) {
  const double temp = std::fmax(Tgas, 10.0);
  const double T5 = temp / 1e5;
  const double ne = n[IDX_E];
  (void)T5; (void)ne;  // unused when the network has none of the species
  double cool = 0.0;
  // collisional ionization: H, He, He+, He(2S)
#ifdef KROME_HAS_H
  cool += 1.27e-21 * std::sqrt(temp) / (1.0 + std::sqrt(T5))
          * std::exp(-1.578091e5 / temp) * ne * n[IDX_H];
#endif
#ifdef KROME_HAS_He
  cool += 9.38e-22 * std::sqrt(temp) / (1.0 + std::sqrt(T5))
          * std::exp(-2.853354e5 / temp) * ne * n[IDX_He];
#endif
#ifdef KROME_HAS_Hej
  cool += 4.95e-22 * std::sqrt(temp) / (1.0 + std::sqrt(T5))
          * std::exp(-6.31515e5 / temp) * ne * n[IDX_Hej];
  cool += 5.01e-27 * std::pow(temp, -0.1687) / (1.0 + std::sqrt(T5))
          * std::exp(-5.5338e4 / temp) * ne * ne * n[IDX_Hej];
#endif
  // recombination: not here; cooling_recomb_ (generated) takes it from
  // the network's own recombination rates (Cen 1992's fits for H+, He+,
  // He++ and the He+ dielectronic recombination until 2.0.0)
  // collisional excitation: H(all n), He+(n=2)
#ifdef KROME_HAS_H
  cool += 7.5e-19 / (1.0 + std::sqrt(T5)) * std::exp(-1.18348e5 / temp)
          * ne * n[IDX_H];
#endif
#ifdef KROME_HAS_Hej
  cool += 9.1e-27 * std::pow(temp, -0.1687) / (1.0 + std::sqrt(T5))
          * std::exp(-1.3179e4 / temp) * ne * ne * n[IDX_Hej];
  cool += 5.54e-17 * std::pow(temp, -0.397) / (1.0 + std::sqrt(T5))
          * std::exp(-4.73638e5 / temp) * ne * n[IDX_Hej];
#endif
  return std::fmax(cool, 0.0);
}
#else
KROME_HD inline double cooling_Atomic(double, const double*) { return 0.0; }
#endif

#ifdef KROME_HAS_H2_ANY
KROME_HD inline double cooling_H2(double Tgas, const double* n) {
  const double temp = Tgas;
  const double T3 = temp * 1e-3;
  const double logt3 = std::log10(T3);
  const double logt = std::log10(temp);
  const double logt32 = logt3 * logt3, logt33 = logt32 * logt3,
               logt34 = logt33 * logt3, logt35 = logt34 * logt3,
               logt36 = logt35 * logt3, logt37 = logt36 * logt3,
               logt38 = logt37 * logt3;
  const double w14 = wCool_(logt, 1.0, 4.0);
  const double w24 = wCool_(logt, 2.0, 4.0);
  (void)w14; (void)w24;  // unused without the species below
  double cool = 0.0;

#ifdef KROME_HAS_H
  double fH2H;  // H2-H (Glover & Abel 2008)
  if (temp <= 1e2)
    fH2H = std::pow(10.0, -16.818342 + 37.383713 * logt3
           + 58.145166 * logt32 + 48.656103 * logt33
           + 20.159831 * logt34 + 3.8479610 * logt35) * n[IDX_H];
  else if (temp <= 1e3)
    fH2H = std::pow(10.0, -24.311209 + 3.5692468 * logt3
           - 11.332860 * logt32 - 27.850082 * logt33
           - 21.328264 * logt34 - 4.2519023 * logt35) * n[IDX_H];
  else if (temp <= 6e3)
    fH2H = std::pow(10.0, -24.311209 + 4.6450521 * logt3
           - 3.7209846 * logt32 + 5.9369081 * logt33
           - 5.5108049 * logt34 + 1.5538288 * logt35) * n[IDX_H];
  else
    fH2H = 1.862314467912518e-22 * wCool_(logt, 1.0, std::log10(6e3))
           * n[IDX_H];
  cool += fH2H;
#endif

#ifdef KROME_HAS_Hj
  double fH2Hp;  // H2-H+
  if (temp > 1e1 && temp <= 1e4)
    fH2Hp = std::pow(10.0, -22.089523 + 1.5714711 * logt3
            + 0.015391166 * logt32 - 0.23619985 * logt33
            - 0.51002221 * logt34 + 0.32168730 * logt35) * n[IDX_Hj];
  else
    fH2Hp = 1.182509139382060e-21 * n[IDX_Hj] * w14;
  cool += fH2Hp;
#endif

  // H2-H2
  cool += w24 * std::pow(10.0, -23.962112 + 2.09433740 * logt3
          - 0.77151436 * logt32 + 0.43693353 * logt33
          - 0.14913216 * logt34 - 0.033638326 * logt35) * nH2_(n);

#ifdef KROME_HAS_E
  double fH2e;  // H2-e
  if (temp <= 5e2)
    fH2e = std::pow(10.0, std::fmin(-21.928796 + 16.815730 * logt3
           + 96.743155 * logt32 + 343.19180 * logt33
           + 734.71651 * logt34 + 983.67576 * logt35
           + 801.81247 * logt36 + 364.14446 * logt37
           + 70.609154 * logt38, 30.0)) * n[IDX_E];
  else
    fH2e = std::pow(10.0, -22.921189 + 1.6802758 * logt3
           + 0.93310622 * logt32 + 4.0406627 * logt33
           - 4.7274036 * logt34 - 8.8077017 * logt35
           + 8.9167183 * logt36 + 6.4380698 * logt37
           - 6.3701156 * logt38) * n[IDX_E];
  cool += fH2e * w24;
#endif

#ifdef KROME_HAS_He
  double fH2He;  // H2-He
  if (temp > 1e1 && temp <= 1e4)
    fH2He = std::pow(10.0, -23.689237 + 2.1892372 * logt3
            - 0.81520438 * logt32 + 0.29036281 * logt33
            - 0.16596184 * logt34 + 0.19191375 * logt35) * n[IDX_He];
  else
    fH2He = 1.002560385050777e-22 * n[IDX_He] * w14;
  cool += fH2He;
#endif

  if (cool <= 0.0) return 0.0;

  // high-density (LTE) limit: HM79 + GP98 below 2e3 K, Glover 2015 above
  double HDL;
  if (temp < 2e3) {
    const double HDLR = 9.5e-22 * std::pow(T3, 3.76)
                        / (1.0 + 0.12 * std::pow(T3, 2.1))
                        * std::exp(-std::pow(0.13 / T3, 3.0))
                        + 3e-24 * std::exp(-0.51 / T3);
    const double HDLV = 6.7e-19 * std::exp(-5.86 / T3)
                        + 1.6e-18 * std::exp(-11.7 / T3);
    HDL = HDLR + HDLV;
  } else if (temp <= 1e4) {
    HDL = std::pow(10.0, -20.584225 + 5.0194035 * logt3
          - 1.5738805 * logt32 - 4.7155769 * logt33
          + 2.4714161 * logt34 + 5.4710750 * logt35
          - 3.9467356 * logt36 - 2.2148338 * logt37
          + 1.8161874 * logt38);
  } else {
    const double dump14 = 1.0
        / (1.0 + std::exp(std::fmin((temp - 3e4) * 2e-4, 3e2)));
    HDL = 5.531333679406485e-19 * dump14;
  }
  if (HDL == 0.0) return 0.0;
  return nH2_(n) / (1.0 / HDL + 1.0 / cool);
}
#else
KROME_HD inline double cooling_H2(double, const double*) { return 0.0; }
#endif

#if defined(KROME_HAS_HD) && defined(KROME_HAS_H)
// HD cooling, Lipovka, Nunez-Lopez & Avila-Reese 2005 (Table 1 polynomial
// in log10 Tgas and log10 n_H), ported verbatim from legacy KROME.
// Extrapolated limits as legacy: Tgas capped at 1e4 K, n_H clamped to
// [1e-2, 1e10]; returns 0 below the CMB temperature p.Tcmb.
template <class P>
KROME_HD inline double cooling_HD(double Tgas, const double* n, const P& p) {
  if (Tgas < p.Tcmb) return 0.0;
  static constexpr double c[5][5] = {
      {-42.56788, 0.92433, 0.54962, -0.07676, 0.00275},
      {21.93385, 0.77952, -1.06447, 0.11864, -0.00366},
      {-10.19097, -0.54263, 0.62343, -0.07366, 0.002514},
      {2.19906, 0.11711, -0.13768, 0.01759, -0.00066631},
      {-0.17334, -0.00835, 0.0106, -0.001482, 0.00006192}};
  const double T = std::fmin(Tgas, 1e4);
  const double dd = std::fmin(std::fmax(n[IDX_H], 1e-2), 1e10);
  const double logT = std::log10(T);
  const double lognH = std::log10(dd);
  double logW = 0.0;
  for (int j = 0; j < 5; ++j) {
    const double lHj = std::pow(lognH, j);
    for (int i = 0; i < 5; ++i)
      logW += c[i][j] * std::pow(logT, i) * lHj;
  }
  return std::pow(10.0, logW) * n[IDX_HD];
}
#else
template <class P>
KROME_HD inline double cooling_HD(double, const double*, const P&) {
  return 0.0;
}
#endif

// Compton cooling off the CMB, Cen 1992 (legacy cooling_compton); the
// redshift comes from p.zredshift (legacy phys_zredshift @common).
#ifdef KROME_HAS_E
template <class P>
KROME_HD inline double cooling_Compton(double Tgas, const double* n,
                              const P& p) {
  const double zp1 = 1.0 + p.zredshift;
  return 5.65e-36 * zp1 * zp1 * zp1 * zp1
      * (Tgas - 2.73 * zp1) * n[IDX_E];  // erg cm^-3 s^-1
}
#else
template <class P>
KROME_HD inline double cooling_Compton(double, const double*, const P&) {
  return 0.0;
}
#endif

// adiabatic expansion cooling, Galli & Palla 1998 (legacy
// cooling_expansion): 3 ntot kB Tgas H0 (1+z) sqrt(omega0 z + 1).
// p.hubble0 [s^-1], p.omega0 and p.zredshift are host-supplied
// (legacy hubble0/omega0/phys_zredshift @common variables).
template <class P>
KROME_HD inline double cooling_expansion(double Tgas, const double* n, const P& p) {
  return 3.0 * ntot_(n) * constants::boltzmann_erg * Tgas * p.hubble0
      * (1.0 + p.zredshift)
      * std::sqrt(p.omega0 * p.zredshift + 1.0);  // erg cm^-3 s^-1
}

// H2 collision-induced emission cooling, fit to Ripamonti & Abel 2004
// (legacy cooling_CIE, ported verbatim): valid 100 K - 1e6 K, cubic
// extrapolation above 1e5 K, RA04 opacity correction via tauCIE.
#ifdef KROME_HAS_H2_ANY
template <class P>
KROME_HD inline double cooling_CIE(double Tgas, const double* n, const P& p) {
  // avoid division by zero in the opacity term (tauCIE -> 0)
  if (nH2_(n) < 1e-12) return 0.0;
  if (Tgas < p.Tcmb) return 0.0;
  const double x = std::log10(Tgas);
  double logcool = -1e99;  // zero outside the fit boundaries
  if (x > 2.0 && x < 2.95) {
    static constexpr double a[6] = {-30.3314216559651, 19.0004016698518,
                                -17.1507937874082, 9.49499574218739,
                                -2.54768404538229, 0.265382965410969};
    logcool = a[0] + x * (a[1] + x * (a[2] + x * (a[3]
              + x * (a[4] + x * a[5]))));
  } else if (x >= 2.95 && x < 5.0) {
    static constexpr double b[6] = {-180.992524120965, 168.471004362887,
                                -67.499549702687, 13.5075841245848,
                                -1.31983368963974, 0.0500087685129987};
    logcool = b[0] + x * (b[1] + x * (b[2] + x * (b[3]
              + x * (b[4] + x * b[5]))));
  } else if (x >= 5.0) {
    logcool = 3.0 * x - 21.2968837223113;  // cubic extrapolation
  }
  const double tauCIE = std::pow(nH2_(n) * 1.4285714e-16, 2.8);  // 1/7e15
  const double cool = constants::p_mass * std::pow(10.0, logcool);  // p_mass
  return cool * std::fmin(1.0, (1.0 - std::exp(-tauCIE)) / tauCIE)
      * nH2_(n) * ntot_(n);  // erg cm^-3 s^-1
}
#else
template <class P>
KROME_HD inline double cooling_CIE(double, const double*, const P&) {
  return 0.0;
}
#endif

// Planck mean opacity fit (Lenzuni+ 1996; legacy kpla: density-only fit,
// valid T < 3e3 K, zero below rho = 1e-12 g/cm3, capped at 0.5 g/cm3)
KROME_HD inline double kpla_(const double* n) {
  const double rho = rho_(n);
  if (rho < 1e-12) return 0.0;
  const double y = std::log10(std::fmin(rho, 0.5));
  return std::pow(10.0, 1.000042 * y + 2.14989);  // cm^2/g
}

// continuum cooling (legacy cooling_Continuum): 4 sigma_SB T^4 kgas rho
// times the tau^-2 escape fraction over a Jeans length
KROME_HD inline double cooling_Continuum(double Tgas, const double* n) {
  const double rho = std::fmax(rho_(n), 1e-40);
  const double kgas = kpla_(n);
  // legacy get_jeans_length: mu = rho/(ntot m_p), constants as legacy
  const double mu = rho / std::fmax(ntot_(n), 1e-40) / constants::p_mass;
  const double lj = std::sqrt(constants::pi * constants::boltzmann_erg * Tgas
                              / rho / constants::p_mass / constants::gravity / mu);
  const double tau = lj * kgas * rho + 1e-40;
  const double beta = std::fmin(1.0, 1.0 / (tau * tau));
  return 4.0 * constants::stefboltz_erg * Tgas * Tgas * Tgas * Tgas
      * kgas * rho * beta;  // erg cm^-3 s^-1
}

} // namespace krome
