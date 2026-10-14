// KROME-X radiation helpers (legacy krome_set_photoBin_* / opacity API),
// network-independent. KROME-X keeps the photon field as the caller's per-cell
// array Phi[nbins] (photon flux per bin [cm^-2 s^-1]); these functions
// build binnings and fields for it. Legacy stored the specific intensity
// J_nu [eV s^-1 cm^-2 Hz^-1 sr^-1] per bin and sampled spectra at the bin
// midpoint; here every field is the bin integral
//   Phi_b = 4 pi / h  int_{Elo}^{Ehi} J_nu(E) d ln E,
// analytically where possible (Draine, J21), else by composite
// Gauss-Legendre quadrature in ln E. A J constant across a bin maps to
// Phi_b = 4 pi J ln(Ehi/Elo) / h (legacy photoBinJ semantics).
// Physical constants from krome::constants (kromex/constants.py).
#pragma once
#include <cmath>
#include <vector>
#include "data_radiation.hpp"

namespace krome {
namespace spectra {

constexpr double H_EV = constants::planck_eV;
constexpr double C_CM = constants::clight;
constexpr double KB_EV = constants::boltzmann_eV;
constexpr double PI4 = 4.0 * constants::pi;

// ---- binnings (legacy krome_set_photobinE_lin / _log; the mid energy is
// not stored in KROME-X: the bins are their edges)
inline void bins_lin(int B, double lower, double upper, double* Elo,
                     double* Ehi) {
  const double dE = (upper - lower) / B;
  for (int b = 0; b < B; ++b) {
    Elo[b] = lower + b * dE;
    Ehi[b] = lower + (b + 1) * dE;
  }
}
inline void bins_log(int B, double lower, double upper, double* Elo,
                     double* Ehi) {
  const double l0 = std::log10(lower), dl = (std::log10(upper) - l0) / B;
  for (int b = 0; b < B; ++b) {
    Elo[b] = std::pow(10.0, l0 + b * dl);
    Ehi[b] = std::pow(10.0, l0 + (b + 1) * dl);
  }
}

// ---- J (constant across the bin, legacy photoBinJ) <-> Phi
inline double J_to_phi(double J, double Elo, double Ehi) {
  return PI4 * J * std::log(Ehi / Elo) / H_EV;
}
inline double phi_to_J(double phi, double Elo, double Ehi) {
  return phi * H_EV / (PI4 * std::log(Ehi / Elo));
}

// composite 4-point Gauss-Legendre of f(E) d ln E over [lo, hi], with
// sub-intervals of at most dmax in ln E
template <class F>
inline double int_dlnE_(F f, double lo, double hi, double dmax = 0.02) {
  if (!(hi > lo)) return 0.0;
  static const double x[4] = {-0.8611363115940526, -0.3399810435848563,
                              0.3399810435848563, 0.8611363115940526};
  static const double w[4] = {0.3478548451374538, 0.6521451548625461,
                              0.6521451548625461, 0.3478548451374538};
  const double a = std::log(lo), L = std::log(hi) - a;
  const int m = std::max(1, static_cast<int>(std::ceil(L / dmax)));
  const double h = L / m;
  double s = 0.0;
  for (int i = 0; i < m; ++i) {
    const double c = a + (i + 0.5) * h;
    for (int q = 0; q < 4; ++q) s += w[q] * f(std::exp(c + 0.5 * h * x[q]));
  }
  return s * 0.5 * h;
}

// ---- spectra: phi[B] for bins [Elo, Ehi]

// blackbody specific intensity B_nu(T) [eV s^-1 cm^-2 Hz^-1 sr^-1] (legacy
// planckBB: 0 beyond E/kT = 300), no dilution
inline double planck_J(double E, double Tbb) {
  const double x = E / (KB_EV * Tbb);
  if (!(x < 300.0) || !(E > 1e-10)) return 0.0;
  return 2.0 * E * E * E / (H_EV * H_EV * C_CM * C_CM) / std::expm1(x);
}
inline void blackbody(int B, const double* Elo, const double* Ehi,
                      double Tbb, double* phi) {
  for (int b = 0; b < B; ++b)
    phi[b] = PI4 / H_EV * int_dlnE_(
        [&](double E) { return planck_J(E, Tbb); }, Elo[b], Ehi[b]);
}

// legacy krome_set_photoBin_BBlog_auto range: from the Rayleigh-Jeans
// energy where B_nu = eps (1e-6) to the energy above the peak where B_nu
// drops to eps (bisection between the legacy Wien-lambda estimate and
// 290 kT); returns 0, or -1 when no such range exists
inline int blackbody_range(double Tbb, double* lower, double* upper,
                           double eps) {
  *lower = H_EV * C_CM * std::sqrt(0.5 * eps / (Tbb * KB_EV));
  double x0 = Tbb / 2.8977721e-1 * C_CM * H_EV;
  double x1 = 290.0 * KB_EV * Tbb;
  if (planck_J(x0, Tbb) < eps || planck_J(x1, Tbb) > eps) return -1;
  double xm = 0.5 * (x0 + x1);
  for (int it = 0; it < 400; ++it) {
    xm = 0.5 * (x0 + x1);
    const double f = planck_J(xm, Tbb) - eps;
    if (std::fabs(f) < 1e-9 * eps) break;
    (f > 0.0 ? x0 : x1) = xm;
  }
  *upper = xm;
  return 0;
}

// Draine 1978 ISRF, photon intensity F(E) = 1.658e6 E - 2.152e5 E^2 +
// 6.919e3 E^3 [cm^-2 s^-1 sr^-1 eV^-1] for 5 < E < 13.6 eV (legacy
// draineLin/Log/Custom): Phi_b = 4 pi int F dE over the overlap, exactly
inline void draine(int B, const double* Elo, const double* Ehi,
                   double* phi) {
  auto P = [](double E) {
    return 1.658e6 * E * E / 2.0 - 2.152e5 * E * E * E / 3.0
        + 6.919e3 * E * E * E * E / 4.0;
  };
  for (int b = 0; b < B; ++b) {
    const double lo = std::fmax(Elo[b], 5.0), hi = std::fmin(Ehi[b], 13.6);
    phi[b] = hi > lo ? PI4 * (P(hi) - P(lo)) : 0.0;
  }
}

// J21 power law: J = J21 1e-21 erg s^-1 cm^-2 Hz^-1 sr^-1 (13.6 eV / E)^alpha
// (legacy J21lin/log: J21 = 1, alpha = 1, no cutoff), exactly
inline void j21(int B, const double* Elo, const double* Ehi, double J21,
                double alpha, double* phi) {
  const double J0 = J21 * 6.2415e-10;  // legacy 1e-21 erg in eV
  for (int b = 0; b < B; ++b) {
    const double I = alpha == 0.0
        ? std::log(Ehi[b] / Elo[b])
        : (std::pow(13.6 / Elo[b], alpha) - std::pow(13.6 / Ehi[b], alpha))
            / alpha;
    phi[b] = PI4 / H_EV * J0 * I;
  }
}

// int over [e0, e1] of J d ln E for J interpolated between the nodes
// (E0, J0), (E1, J1): linear in E, or a power law (log-log; 0 when a node
// is not positive), exactly
inline double seg_dlnE_(double e0, double e1, double E0, double J0,
                        double E1, double J1, bool loglog) {
  if (!(e1 > e0)) return 0.0;
  if (!loglog) {
    const double s = (J1 - J0) / (E1 - E0);
    return (J0 - s * E0) * std::log(e1 / e0) + s * (e1 - e0);
  }
  if (!(J0 > 0.0) || !(J1 > 0.0)) return 0.0;
  const double beta = std::log(J1 / J0) / std::log(E1 / E0);
  const double L = std::log(e1 / e0);
  if (std::fabs(beta * L) < 1e-8)  // expm1 series limit
    return J0 * std::pow(e0 / E0, beta) * L * (1.0 + 0.5 * beta * L);
  return J0 * std::pow(e0 / E0, beta) * std::expm1(beta * L) / beta;
}

// tabulated J(E) [eV s^-1 cm^-2 Hz^-1 sr^-1] at nE energies E [eV]
// (ascending), linear or log-log interpolation between the nodes, zero
// outside the table, integrated exactly over each bin (legacy
// krome_load_photoBin_file_2col averaged J linearly over the bin)
inline void table(int B, const double* Elo, const double* Ehi, int nE,
                  const double* Et, const double* Jt, bool loglog,
                  double* phi) {
  for (int b = 0; b < B; ++b) {
    double s = 0.0;
    for (int i = 0; i < nE - 1; ++i) {
      const double e0 = std::fmax(Elo[b], Et[i]);
      const double e1 = std::fmin(Ehi[b], Et[i + 1]);
      s += seg_dlnE_(e0, e1, Et[i], Jt[i], Et[i + 1], Jt[i + 1], loglog);
    }
    phi[b] = PI4 / H_EV * s;
  }
}

// Haardt & Madau 2012 UV background at redshift z (legacy
// krome_set_photoBin_HMlog / HMCustom: bilinear in (z, log10 E) of log10 J,
// there sampled at the bin midpoint). Here: log10 J linear in z between the
// two bracketing redshifts, then a power law between the energy nodes,
// integrated exactly over each bin; zero outside the tabulated
// 0.1237 eV - 49.97 MeV. Returns 0, or -1 for z outside [0, 15.66]
inline int hm2012(int B, const double* Elo, const double* Ehi, double z,
                  double* phi) {
  if (z < hm12_z[0] || z > hm12_z[HM12_NZ - 1]) return -1;
  int i = 0;
  while (i < HM12_NZ - 2 && hm12_z[i + 1] <= z) ++i;
  const double fz = (z - hm12_z[i]) / (hm12_z[i + 1] - hm12_z[i]);
  const double* r0 = hm12_logJ + static_cast<size_t>(i) * HM12_NE;
  const double* r1 = r0 + HM12_NE;
  std::vector<double> E(HM12_NE), J(HM12_NE);
  for (int j = 0; j < HM12_NE; ++j) {
    E[j] = std::pow(10.0, hm12_logE[j]);
    J[j] = std::pow(10.0, r0[j] + fz * (r1[j] - r0[j]));
  }
  table(B, Elo, Ehi, HM12_NE, E.data(), J.data(), true, phi);
  return 0;
}

// ---- dust opacity per bin

// bin-mean dust opacity [cm^2 g^-1] from a table kappa(E) (E [eV]
// ascending): trapezoid of the linearly interpolated kappa over each bin
// (legacy load_opacity_table, with robust node handling). Returns 0, or
// -1 when a bin lies outside the table (legacy stopped)
inline int dust_kappa_bins(int nE, const double* Et, const double* kt,
                           int B, const double* Elo, const double* Ehi,
                           double* kbin) {
  auto kap = [&](double E) {
    int i = 0;
    while (i < nE - 2 && Et[i + 1] < E) ++i;
    const double f = (E - Et[i]) / (Et[i + 1] - Et[i]);
    return kt[i] + f * (kt[i + 1] - kt[i]);
  };
  for (int b = 0; b < B; ++b) {
    if (Elo[b] < Et[0] || Ehi[b] > Et[nE - 1]) return -1;
    double s = 0.0, e0 = Elo[b], k0 = kap(e0);
    for (int i = 0; i < nE; ++i) {
      if (Et[i] <= Elo[b] || Et[i] >= Ehi[b]) continue;
      const double k1 = kt[i];
      s += 0.5 * (k0 + k1) * (Et[i] - e0);
      e0 = Et[i];
      k0 = k1;
    }
    s += 0.5 * (k0 + kap(Ehi[b])) * (Ehi[b] - e0);
    kbin[b] = s / (Ehi[b] - Elo[b]);
  }
  return 0;
}

// ... from the embedded Draine R_V = 3.1 kappa_abs (legacy
// krome_load_average_kabs, kabs_draine_Rv31.dat, micron -> eV)
inline int dust_kappa_draine(int B, const double* Elo, const double* Ehi,
                             double* kbin) {
  std::vector<double> E(KABS_DRAINE_N), k(KABS_DRAINE_N);
  for (int i = 0; i < KABS_DRAINE_N; ++i) {  // descending lambda: ascending E
    E[i] = H_EV * C_CM / (kabs_draine_lam[i] * 1e-4);
    k[i] = kabs_draine_kap[i];
  }
  return dust_kappa_bins(KABS_DRAINE_N, E.data(), k.data(), B, Elo, Ehi,
                         kbin);
}

// ---- G0 / Av finder (legacy krome_find_G0_Av / estimate_G0_Av): fit
//   ln(Phi_b / Phi_Draine,b) = ln G0 - kappa_b d2g 1.8e21 m_p Av / xH
// over the bins where the Draine field is nonzero, by least squares, with
// the legacy priors (ln G0 clamped to [-7, 7] and Av refitted; Av < 0 ->
// Av = 0, ln G0 = the mean clamped to [-7, 7]). xH: H mass fraction (Hnuclei m_p / rho). Legacy
// cached d2g, composition and binning from its first call; here every
// call uses its arguments. Returns 0, or -1 with fewer than 2 bins.
inline int find_G0_Av(int B, const double* Elo, const double* Ehi,
                      const double* phi, const double* kbin, double d2g,
                      double xH, double* G0, double* Av) {
  std::vector<double> pd(B);
  draine(B, Elo, Ehi, pd.data());
  std::vector<double> xs, ys;
  for (int b = 0; b < B; ++b)
    if (pd[b] > 0.0) {
      xs.push_back(-kbin[b] * 1.8e21 * constants::p_mass * d2g / xH);
      ys.push_back(std::log(phi[b] + 1e-200) - std::log(pd[b]));
    }
  const int m = static_cast<int>(xs.size());
  if (m < 2) return -1;
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (int i = 0; i < m; ++i) {
    sx += xs[i]; sy += ys[i]; sxx += xs[i] * xs[i]; sxy += xs[i] * ys[i];
  }
  const double den = m * sxx - sx * sx;
  double a = den != 0.0 ? (m * sxy - sx * sy) / den : 0.0;
  double lg = (sy - a * sx) / m;
  if (lg < -7.0 || lg > 7.0) {
    lg = std::fmin(std::fmax(lg, -7.0), 7.0);
    double num = 0.0;
    for (int i = 0; i < m; ++i) num += (ys[i] - lg) * xs[i];
    a = sxx > 0.0 ? num / sxx : 0.0;
  }
  if (a < 0.0) {  // Av = 0: the best ln G0 is the mean, within the prior
    a = 0.0;       // (legacy took the bare mean, escaping [-7, 7])
    lg = std::fmin(std::fmax(sy / m, -7.0), 7.0);
  }
  *Av = a;
  *G0 = std::exp(lg);
  return 0;
}

} // namespace spectra
} // namespace krome
