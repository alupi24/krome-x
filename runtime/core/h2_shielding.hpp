// KROME-X shielding / Lyman-Werner module, ported verbatim from legacy
// KROME (src/krome_phfuncs.f90, krome_photo.f90, krome_subs.f90).
// Include the generated net_<name>.hpp BEFORE this header.
#pragma once
#include <cmath>

namespace krome {

#ifdef KROME_EXTERN_fselfH2
// H2 self-shielding factor, Draine & Bertoldi 1996 (eq. 37):
// N = H2 column density [cm^-2], b = Doppler broadening [cm s^-1]
KROME_HD inline double fselfH2(double N, double b, double Tgas, const double* n,
                      const Params& p) {
  (void)Tgas; (void)n; (void)p;
  const double x = N * 2e-15;   // N / 5e14
  const double b5 = b * 1e-5;
  return 0.965 / ((1.0 + x / b5) * (1.0 + x / b5))
         + 0.035 / std::sqrt(1.0 + x)
           * std::exp(std::fmax(-8.5e-4 * std::sqrt(1.0 + x), -250.0));
}
#endif

#ifdef KROME_EXTERN_fselfH2_WG11
// H2 self-shielding, Wolcott-Green & Haiman 2011 (DB96 with exponent 1.1
// and thermal Doppler b from Tgas; legacy uses the H mass here, kept for
// bit-compatibility): N = H2 column density [cm^-2]
KROME_HD inline double fselfH2_WG11(double N, double Tgas, const double* n,
                           const Params& p) {
  (void)n; (void)p;
  const double H_mass = constants::p_mass + constants::e_mass;  // legacy p+e
  const double b5 = std::sqrt(constants::boltzmann_erg * Tgas / H_mass) * 1e-5;
  const double x = N * 2e-15;
  return 0.965 / std::pow(1.0 + x / b5, 1.1)
         + 0.035 / std::sqrt(1.0 + x)
           * std::exp(-8.5e-4 * std::sqrt(1.0 + x));
}
#endif

#ifdef KROME_EXTERN_fselfH2_R14
// H2 self-shielding, Richings, Schaye & Oppenheimer 2014 (II), T-dependent:
// N = H2 column density [cm^-2], bturb = turbulent Doppler b [cm s^-1]
// (Richings' default is 7.1 km/s). NOTE: legacy KROME's calc_H2shieldR14
// has `Tgas>=3d3 .or. Tgas<4d3` (always true), so its T>=4000 K branch was
// dead code; this port implements the intended piecewise from the paper.
KROME_HD inline double fselfH2_R14(double N, double bturb, double Tgas,
                          const double* n, const Params& p) {
  (void)n; (void)p;
  const double H_mass = constants::p_mass + constants::e_mass;  // legacy p+e
  const double btherm2 = constants::boltzmann_erg * Tgas / H_mass;
  const double b5 = std::sqrt(btherm2 + bturb * bturb) * 1e-5;
  const double omega = 0.013 * std::pow(1.0 + std::pow(Tgas / 2.7e3, 1.3),
                                        1.0 / 1.3)
                       * std::exp(-std::pow(Tgas / 3.9e3, 14.6));
  double alpha, Ncrit;
  if (Tgas < 3e3) {
    alpha = 1.4;
    Ncrit = 1.3 * (1.0 + std::pow(Tgas / 6e2, 0.8));
  } else if (Tgas < 4e3) {
    alpha = std::pow(Tgas / 4.5e3, -0.8);
    Ncrit = std::pow(Tgas / 4.76e3, -3.8);
  } else {
    alpha = 1.1;
    Ncrit = 2.0;
  }
  const double x = N * 1e-14 / Ncrit;
  return (1.0 - omega) / std::pow(1.0 + x / b5, alpha)
             * std::exp(-5e-7 * (1.0 + x))
         + omega / std::sqrt(1.0 + x)
             * std::exp(-8.5e-4 * std::sqrt(1.0 + x));
}
#endif

#ifdef KROME_EXTERN_H2_sigmaLW_R14
// Richings+14 (I, eq. 2.7) Solomon process recast as an effective flat
// cross section: Gamma = 7.5e-11 s^-1 * n_gamma(12.24-13.51 eV)/2.256e-4
// = sigma0 * Phi(12.24-13.51 eV) with sigma0 = 7.5e-11/(2.256e-4 * c).
// Use as photo(xsec = H2_sigmaLW_R14(E_eV)); the band limits are built in.
KROME_HD inline double H2_sigmaLW_R14(double E_eV) {
  const double sigma0 = 7.5e-11 / (2.256e-4 * constants::clight);  // cm^2
  return (E_eV >= 12.24 && E_eV <= 13.51) ? sigma0 : 0.0;
}
#endif

#ifdef KROME_EXTERN_H2_solomonLW
// Solomon-process H2 photodissociation rate from the radiation intensity
// at 12.87 eV (flux_eV in eV cm^-2 s^-1 Hz^-1 sr^-1, as legacy):
// k = 1.38e9 * J(12.87 eV) [converted to erg]
KROME_HD inline double H2_solomonLW(double flux_eV, double Tgas, const double* n,
                           const Params& p) {
  (void)Tgas; (void)n; (void)p;
  const double eV_to_erg = constants::eV_to_erg;
  return 1.38e9 * flux_eV * eV_to_erg;
}
#endif

#ifdef KROME_EXTERN_H2_sigmaLW
// effective H2 Lyman-Werner cross section [cm^2] vs photon energy,
// ortho/para-weighted sum of line + continuum fits (legacy H2_sigmaLW).
// NOTE: the ortho:para ratio is fixed at the legacy default 3:1 here; the
// xsec is baked into the sigma_rb table at init, so a runtime-variable
// ratio needs a table rebuild (future option).
KROME_HD inline double H2_sigmaLW(double E_eV) {
  const double ortho_para = 3.0;
  double sL0 = 0.0, sW0 = 0.0, sL1 = 0.0, sW1 = 0.0;
  if (E_eV > 14.675 && E_eV < 16.820)
    sL0 = 1e-18 * std::pow(10.0, 15.1289 - 1.05139 * E_eV);
  else if (E_eV > 16.820 && E_eV < 17.6)
    sL0 = 1e-18 * std::pow(10.0, -31.41 + 1.8042e-2 * E_eV * E_eV * E_eV
                                  - 4.2339e-5 * std::pow(E_eV, 5.0));
  if (E_eV > 14.675 && E_eV < 17.7)
    sW0 = 1e-18 * std::pow(10.0, 13.5311 - 0.9182618 * E_eV);
  if (E_eV > 14.159 && E_eV < 15.302)
    sL1 = 1e-18 * std::pow(10.0, 12.0218406 - 0.819429 * E_eV);
  else if (E_eV > 15.302 && E_eV < 17.2)
    sL1 = 1e-18 * std::pow(10.0, 16.04644 - 1.082438 * E_eV);
  if (E_eV > 14.159 && E_eV < 17.2)
    sW1 = 1e-18 * std::pow(10.0, 12.87367 - 0.85088597 * E_eV);
  const double fact = 1.0 / (ortho_para + 1.0);
  return fact * (sL0 + sW0) + (1.0 - fact) * (sL1 + sW1);
}
#endif

} // namespace krome
