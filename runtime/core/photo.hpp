// KROME-X photo module: photon bins, bin-averaged cross sections, photo
// rates from the current radiation field, and per-bin opacity.
//
// Include the generated net_<name>.hpp BEFORE this header (it provides
// NPHOTO, photo_xsec[], photo_reaction[], photo_reactant[], photo_Emin/Emax).
//
// Conventions:
//   Phi_b   photon flux in bin b [photons cm^-2 s^-1], flat in E inside bin
//   sigma_rb = (1/w_b) * Int_{bin b  intersect [Emin_r, Emax_r]} sigma_r(E) dE
//   k_r     = sum_b sigma_rb * Phi_b                       [s^-1]
//   kabs_b  = sum_r sigma_rb * n[reactant_r]               [cm^-1]
//   (photon number density in bin: Phi_b / c)
// Radiative channels (spec 5.4): a photo rate may carry a shielding
// factor s_r (photo(shield=...)), k_r = s_r sum_b sigma_rb Phi_b, which
// also scales its photoheating; a J_LW rate k = khat J_LW, J_LW =
// sum_b lww_b Phi_b, acts on the bins as the cross section khat lww_b.
// Each event absorbs 1/pdiss photons, so the photon removal is
//   kabs_b = sum_r s_r sigma_rb n_r / pdiss_r + sum_q khat_q lww_b n_q / pdiss_q
#pragma once
#include <cmath>
#include <vector>

#define KROME_PHOTO_HPP_INCLUDED 1

namespace krome {

constexpr double CLIGHT = constants::clight;  // cm s^-1

struct PhotoSums;

// raw-pointer view of a PhotoModule's tables, in host or device memory: the
// photo kernels below take it, so they run wherever the data lives (the
// PhotoModule methods pass their own view; a device back-end passes one of
// device copies)
struct PhotoView {
  int B = 0;
  double Esec = 100.0;
  const double* sigma = nullptr;    // NPHOTO x B
  const double* sigmaE = nullptr;
  const double* sigmaEs = nullptr;
  const double* sigmaCs = nullptr;
  const double* sigmaPs = nullptr;
  const double* habw = nullptr;     // B
  const double* Elo = nullptr;
  const double* Ehi = nullptr;
  const double* lww = nullptr;
};

// the precomputed bin sums of a static field (PhotoSums) as a view
struct PhotoSumsView {
  const double* s = nullptr;   // [NPHOTO]
  const double* se = nullptr;
  const double* ss = nullptr;
  const double* sc = nullptr;
  const double* sp = nullptr;
  double hab = 0.0;
};

struct Radiation {
  std::vector<double> Elo, Ehi;  // bin edges [eV]
  std::vector<double> Phi;       // photon flux per bin [cm^-2 s^-1]
  // Doppler broadening of the tabulated (file) cross sections, applied
  // once when the bins are built (legacy -useBroadening with the Tgas of
  // krome_set_photobinE_*): temperature [K] and turbulent velocity
  // [cm s^-1]; 0 and 0 = none
  double broad_T = 0.0;
  double broad_vturb = 0.0;
  int nbins() const { return static_cast<int>(Phi.size()); }

  // B logarithmic bins between emin and emax [eV], zero field
  static Radiation log_bins(int B, double emin, double emax) {
    Radiation r;
    double f = std::pow(emax / emin, 1.0 / B);
    double e = emin;
    for (int b = 0; b < B; ++b, e *= f) {
      r.Elo.push_back(e);
      r.Ehi.push_back(e * f);
      r.Phi.push_back(0.0);
    }
    return r;
  }
};

struct PhotoModule {
  int B = 0;
  double Esec = 100.0;          // photoelectron energy [eV] above which the
                                // SvS85 secondary split applies (below: heat)
  std::vector<double> sigma;    // NPHOTO x B, bin-averaged [cm^2]
  std::vector<double> sigmaE;   // NPHOTO x B, sigma*(E - Eth) avg [cm^2 eV]
  std::vector<double> sigmaEs;  // as sigmaE, fast-photoelectron part only
  std::vector<double> sigmaCs;  // sigma count, fast part only [cm^2]
  std::vector<double> sigmaPs;  // sigma * E(photon), fast part [cm^2 eV]
  std::vector<double> habw;     // per-bin weights of the Habing integral
  std::vector<double> Elo, Ehi; // bin edges [eV] (photo_J_ for pumping)
  std::vector<double> lww;      // J_LW = sum_b lww_b Phi_b (photo_J_LW_)

  // Habing flux G [Habing units] of the current field: energy flux in
  // 6-13.6 eV over 1.6e-3 erg cm^-2 s^-1 (legacy GHabing_thin, but with
  // exact fractional bin overlap instead of the bin-midpoint test)
  double habing(const double* Phi, const PhotoSums* ps = nullptr) const;

  // the per-reaction bin sums of a field held fixed during a solve
  // (PhotoSums), computed once instead of on every rhs evaluation
  void sums(const double* Phi, PhotoSums& ps) const;

  void init(const Radiation& rad, int quad_pts = 64) {
    B = rad.nbins();
    Elo = rad.Elo;
    Ehi = rad.Ehi;
    const size_t NB = static_cast<size_t>(NPHOTO) * B;
    sigma.assign(NB, 0.0);
    sigmaE.assign(NB, 0.0);
    sigmaEs.assign(NB, 0.0);
    sigmaCs.assign(NB, 0.0);
    sigmaPs.assign(NB, 0.0);
    habw.assign(B, 0.0);
    lww.assign(B, 0.0);
    for (int b = 0; b < B; ++b) {
      const double a = std::fmax(rad.Elo[b], LW_EMIN);
      const double c = std::fmin(rad.Ehi[b], LW_EMAX);
      if (c > a)  // photo_Jbin_(b) ln(c/a) / ln(LW_EMAX/LW_EMIN) per Phi_b
        lww[b] = constants::planck_eV * std::log(c / a)
            / (4.0 * constants::pi * std::log(rad.Ehi[b] / rad.Elo[b])
               * std::log(LW_EMAX / LW_EMIN));
    }
    for (int b = 0; b < B; ++b) {
      const double w = rad.Ehi[b] - rad.Elo[b];
      const double a = std::fmax(rad.Elo[b], 6.0);
      const double c = std::fmin(rad.Ehi[b], 13.6);
      if (c > a)  // Int E * (Phi/w) dE over the overlap, in Habing units
        habw[b] = (c * c - a * a) / (2.0 * w) * constants::eV_to_erg / 1.6e-3;
    }
    for (int r = 0; r < NPHOTO; ++r) {
      if (photo_binned[r] >= 0) {  // lines given per bin (h2pd): no
        // photoelectron (Solomon: no deposit), sigmaE* stay 0
        photo_bin_xsec(r, B, rad.Elo.data(), rad.Ehi.data(),
                       sigma.data() + static_cast<size_t>(r) * B);
        continue;
      }
      // Doppler broadening (tabulated cross sections): each E sees the
      // mean of the table over [E - w/2, E + w/2], w = E b / c, b^2 =
      // 2 kB T / m + v_turb^2 (legacy's box width; legacy widened each
      // table segment instead, which inflated the integrated cross
      // section, and dropped segments overlapping a bin's left edge)
      const double b2 = (rad.broad_T > 0.0
                             ? 2.0 * constants::boltzmann_erg * rad.broad_T
                                   / species_mass_g[photo_reactant[r]]
                             : 0.0) + rad.broad_vturb * rad.broad_vturb;
      const double dsh = photo_xsec_int[r] != nullptr
          ? std::sqrt(b2) / CLIGHT : 0.0;
      auto xs = [&](double E) {
        if (dsh > 0.0) {
          const double wd = dsh * E;
          return photo_xsec_int[r](E - 0.5 * wd, E + 0.5 * wd) / wd;
        }
        return photo_xsec[r](E);
      };
      // the broadened cross section reaches w/2 beyond the table's ends
      const double Elo_r = photo_Emin[r] * (1.0 - 0.5 * dsh);
      const double Ehi_r = photo_Emax[r] * (1.0 + 0.5 * dsh);
      for (int b = 0; b < B; ++b) {
        const double w = rad.Ehi[b] - rad.Elo[b];
        const double a = std::fmax(rad.Elo[b], Elo_r);
        const double c = std::fmin(rad.Ehi[b], Ehi_r);
        if (c <= a) continue;
        double s = 0.0, se = 0.0, ss = 0.0, sc = 0.0, sp = 0.0;
        const double dE = (c - a) / quad_pts;  // midpoint rule over [a, c]
        for (int q = 0; q < quad_pts; ++q) {
          const double E = a + (q + 0.5) * dE;
          const double x = xs(E) * dE;
          // photoelectron energy (>= 0: a broadened threshold)
          const double Ee = std::fmax(E - photo_Emin[r], 0.0);
          s += x;
          se += x * Ee;
          if (Ee > Esec) {
            ss += x * Ee;
            sc += x;
            sp += x * E;
          }
        }
        const size_t i = static_cast<size_t>(r) * B + b;
        sigma[i] = s / w;
        sigmaE[i] = se / w;
        sigmaEs[i] = ss / w;
        sigmaCs[i] = sc / w;
        sigmaPs[i] = sp / w;
      }
    }
  }

  // raw-pointer view of the tables (valid while this module lives)
  PhotoView view() const {
    PhotoView v;
    v.B = B;
    v.Esec = Esec;
    v.sigma = sigma.data();
    v.sigmaE = sigmaE.data();
    v.sigmaEs = sigmaEs.data();
    v.sigmaCs = sigmaCs.data();
    v.sigmaPs = sigmaPs.data();
    v.habw = habw.data();
    v.Elo = Elo.data();
    v.Ehi = Ehi.data();
    v.lww = lww.data();
    return v;
  }

  // add photo rate coefficients into k[] (rates() left them at 0);
  // sh: the shielding factors (photo_shields), nullptr for none
  // ps: the precomputed sums of this field (PhotoModule::sums), or nullptr
  void rates(const double* Phi, double* k, const double* sh = nullptr,
             const PhotoSums* ps = nullptr) const;

  // photoheating rate [erg cm^-3 s^-1]: each absorption deposits the
  // photoelectron energy (E - Eth), except in the Solomon process
  // (photo_deposit false: line pumping, heated by CHEM / PHOTOAV). With
  // radiation.secondary_ionization, secondary() splits the deposit instead.
  double heating(const double* Phi, const double* n,
                 const double* sh = nullptr,
                 const PhotoSums* ps = nullptr) const;

#ifdef KROME_HAS_SECION
  // secondary ionization + heat split (photo_secondary)
  double secondary(const double* Phi, const double* n, double* dn_add,
                   const double* sh = nullptr,
                   const PhotoSums* ps = nullptr) const;
#endif

  // per-bin photon removal coefficient [cm^-1] (photo_absorption)
  void absorption(const double* n, const double* sh, const double* kh,
                  double* kabs) const;
};

// per-reaction bin sums of a field held fixed during a solve: s_r =
// sum_b sigma_rb Phi_b and the energy-weighted variants (the rates,
// photoheating, secondary ionization) plus the Habing flux. Each sum is
// accumulated in the same order as the direct evaluation, so results are
// bitwise identical.
struct PhotoSums {
  std::vector<double> s, se, ss, sc, sp;  // [NPHOTO]
  double hab = 0.0;
  PhotoSumsView view() const {
    PhotoSumsView v;
    v.s = s.data();
    v.se = se.data();
    v.ss = ss.data();
    v.sc = sc.data();
    v.sp = sp.data();
    v.hab = hab;
    return v;
  }
};

// ---- the photo kernels (host and device) --------------------------------

// Habing flux G [Habing units]: energy flux in 6-13.6 eV over
// 1.6e-3 erg cm^-2 s^-1 (ps: the precomputed value, or nullptr)
KROME_HD inline double photo_habing(const PhotoView& v, const double* Phi,
                                    const PhotoSumsView* ps = nullptr) {
  if (ps != nullptr) return ps->hab;
  double g = 0.0;
#pragma omp simd reduction(+:g)
  for (int b = 0; b < v.B; ++b) g += v.habw[b] * Phi[b];
  return g;
}

// the per-reaction bin sums of a field into the arrays s .. sp [NPHOTO]
KROME_HD inline void photo_sums(const PhotoView& v, const double* Phi,
                                double* s_, double* se_, double* ss_,
                                double* sc_, double* sp_) {
  for (int r = 0; r < NPHOTO; ++r) {
    double s = 0.0, se = 0.0, ss = 0.0, sc = 0.0, sp = 0.0;
#pragma omp simd reduction(+:s, se, ss, sc, sp)
    for (int b = 0; b < v.B; ++b) {
      const size_t i = static_cast<size_t>(r) * v.B + b;
      s += v.sigma[i] * Phi[b];
      se += v.sigmaE[i] * Phi[b];
      ss += v.sigmaEs[i] * Phi[b];
      sc += v.sigmaCs[i] * Phi[b];
      sp += v.sigmaPs[i] * Phi[b];
    }
    s_[r] = s;
    se_[r] = se;
    ss_[r] = ss;
    sc_[r] = sc;
    sp_[r] = sp;
  }
}

// photo rate coefficients into k[] (shields sh or nullptr)
KROME_HD inline void photo_rates(const PhotoView& v, const double* Phi,
                                 double* k, const double* sh = nullptr,
                                 const PhotoSumsView* ps = nullptr) {
  for (int r = 0; r < NPHOTO; ++r) {
    double s;
    if (ps != nullptr) {
      s = ps->s[r];
    } else {
      s = 0.0;
#pragma omp simd reduction(+:s)
      for (int b = 0; b < v.B; ++b)
        s += v.sigma[static_cast<size_t>(r) * v.B + b] * Phi[b];
    }
    // += : several channels may feed one reaction (h2pd); rates() left
    // the photo reactions at 0
    k[photo_reaction[r]] += sh != nullptr ? sh[r] * s : s;
  }
}

// photoheating rate [erg cm^-3 s^-1] (see PhotoModule::heating)
KROME_HD inline double photo_heating(const PhotoView& v, const double* Phi,
                                     const double* n,
                                     const double* sh = nullptr,
                                     const PhotoSumsView* ps = nullptr) {
  const double eV_to_erg = constants::eV_to_erg;
  double h = 0.0;
  for (int r = 0; r < NPHOTO; ++r) {
    double s;
    if (ps != nullptr) {
      s = ps->se[r];
    } else {
      s = 0.0;
#pragma omp simd reduction(+:s)
      for (int b = 0; b < v.B; ++b)
        s += v.sigmaE[static_cast<size_t>(r) * v.B + b] * Phi[b];
    }
    if (photo_deposit[r])
      h += s * n[photo_reactant[r]] * (sh != nullptr ? sh[r] : 1.0);
  }
  return h * eV_to_erg;
}

#ifdef KROME_HAS_SECION
// the four secondary-ionization sums of reaction r (precomputed or not)
KROME_HD inline void photo_sums_secondary_(const PhotoView& v,
                                           const double* Phi, int r,
                                           const PhotoSumsView* ps,
                                           double& se, double& ss,
                                           double& sc, double& sp) {
  if (ps != nullptr) {
    se = ps->se[r];
    ss = ps->ss[r];
    sc = ps->sc[r];
    sp = ps->sp[r];
    return;
  }
  se = ss = sc = sp = 0.0;
#pragma omp simd reduction(+:se, ss, sc, sp)
  for (int b = 0; b < v.B; ++b) {
    const size_t i = static_cast<size_t>(r) * v.B + b;
    se += v.sigmaE[i] * Phi[b];
    ss += v.sigmaEs[i] * Phi[b];
    sc += v.sigmaCs[i] * Phi[b];
    sp += v.sigmaPs[i] * Phi[b];
  }
}

// Secondary ionization and heating by fast photoelectrons, following
// Quadri et al. 2025 (arXiv:2505.05556, Appendix B): the Dalgarno+99
// results refit in the Shull & van Steenberg 1985 form, valid for an
// H/H2/He mixture over the full ionized-fraction range.
//   x_e = n_e / n_Htot ;  x_eff = 1.83 x_e / (1 + 0.83 x_e)
//   phi_i(E) = (E/Eth_i - 1) * A_i(x_e, r),  r = n_H2/n_H
// NOTE: the paper's eq. B.7 has 1.89 n_H2/n_H in the phi_H denominator
// while its carbon example (B.15) uses 1.82 — we follow the fit (1.89).
//
// volumetric secondary ionization rate of a fit species i:
//   R_i = A_i * sum_j n_j sum_b sigma_jb(fast) Phi_b (E/Eth_i - 1)
//       = A_i * (P/Eth_i - N)
// with P, N the photon-energy-weighted and counting deposit sums.
//
// secondary ionization + heat split: adds the secondary ionization
// fluxes [cm^-3 s^-1] into dn_add and returns the net photoheating rate
// [erg cm^-3 s^-1]. Slow photoelectrons (< Esec) heat fully; fast ones
// heat with the mixture efficiency eta (paper eqs. B.18-B.21). The
// heating deposit uses the exact photoelectron energy (E - Eth), which
// the paper approximates as E in the X-ray regime.
// the heating efficiency of the fast (secondary-ionizing) part of the
// photoelectron deposit for the mixture (Quadri+25 eqs. B.18-B.21)
KROME_HD inline double photo_heat_eta_(const double* n) {
  const double xe = std::fmin(std::fmax(
      n[IDX_E] / std::fmax(Hnuclei_(n), 1e-300), 0.0), 1.0);
  const double xeff = 1.83 * xe / (1.0 + 0.83 * xe);
#ifdef KROME_HAS_H2_ANY
  const double rr = nH2_(n) / std::fmax(n[IDX_H], 1e-300);
#else
  const double rr = 0.0;
#endif
  const double etaHeH = 1.0 - 0.8957
      * std::pow(1.0 - std::pow(xe, 0.5078), 2.9620);
  const double etaH2He = 1.0 - 0.9377
      * std::pow(1.0 - std::pow(xeff, 0.3428), 1.4169);
  return (10.0 * rr * etaH2He + etaHeH) / (10.0 * rr + 1.0);
}

KROME_HD inline double photo_secondary(const PhotoView& v, const double* Phi,
                                       const double* n, double* dn_add,
                                       const double* sh = nullptr,
                                       const PhotoSumsView* ps = nullptr) {
  const double eV_to_erg = constants::eV_to_erg;
  double D_tot = 0.0, D_hi = 0.0;  // deposits [eV cm^-3 s^-1]
  double P = 0.0, N = 0.0;         // photon-E-weighted, count (fast)
  for (int r = 0; r < NPHOTO; ++r) {
    if (!photo_deposit[r]) continue;  // Solomon: no photoelectron
    double se, ss, sc, sp;
    photo_sums_secondary_(v, Phi, r, ps, se, ss, sc, sp);
    const double nr = n[photo_reactant[r]] * (sh != nullptr ? sh[r] : 1.0);
    D_tot += se * nr;
    D_hi += ss * nr;
    N += sc * nr;
    P += sp * nr;
  }
  const double xe = std::fmin(std::fmax(
      n[IDX_E] / std::fmax(Hnuclei_(n), 1e-300), 0.0), 1.0);
  const double xeff = 1.83 * xe / (1.0 + 0.83 * xe);
  const double nH = n[IDX_H];
#ifdef KROME_HAS_H2_ANY
  const double rr = nH2_(n) / std::fmax(nH, 1e-300);
#else
  const double rr = 0.0;
#endif
  const double AH = 0.3443
      * std::pow(1.0 - std::pow(xe, 0.6649), 4.3928) / (1.0 + 1.89 * rr);
  const double AH2 = rr > 0.0
      ? 0.3749 * std::pow(1.0 - std::pow(xeff, 0.5805), 2.3414)
        / (1.0 + 0.53 / rr)
      : 0.0;
  const double AHe = 0.0509
      * std::pow(1.0 - std::pow(xe, 0.7883), 4.9301);
  const double SH = std::fmax(P / 13.6 - N, 0.0);
  for (int t = 0; t < N_SECION; ++t) {
    double R;
    switch (secion_kind[t]) {
      case 0: R = AH * SH; break;
      case 1: R = AH2 * std::fmax(P / 15.4 - N, 0.0); break;
      case 2: R = AHe * std::fmax(P / 24.6 - N, 0.0); break;
      default:  // metal: rescaled per-H-atom rate (paper eq. B.14)
        R = secion_ratio[t] * n[secion_target[t]]
            * AH * SH / std::fmax(nH, 1e-300);
    }
    dn_add[secion_target[t]] -= R;
    dn_add[secion_ion[t]] += R;
    dn_add[IDX_E] += R;
  }
  // heating efficiency for the mixture (paper eqs. B.18-B.21)
  const double eta = photo_heat_eta_(n);
  (void)v;
  return ((D_tot - D_hi) + eta * D_hi) * eV_to_erg;
}
#endif

// per-bin photon removal coefficient [cm^-1] of the photochemistry: the
// photo reactions (shielding sh, nullptr for none) and the J_LW rates
// (unit rates kh = lw_unit_rates, nullptr for none), each event absorbing
// 1/pdiss photons
KROME_HD inline void photo_absorption(const PhotoView& v, const double* n,
                                      const double* sh, const double* kh,
                                      double* kabs) {
  for (int b = 0; b < v.B; ++b) kabs[b] = 0.0;
  for (int r = 0; r < NPHOTO; ++r) {
    const double f = n[photo_reactant[r]]
        * (sh != nullptr ? sh[r] : 1.0) / photo_pdiss[r];
    for (int b = 0; b < v.B; ++b)
      kabs[b] += v.sigma[static_cast<size_t>(r) * v.B + b] * f;
  }
  if (kh != nullptr)
    for (int q = 0; q < NLW; ++q) {
      const double f = n[lw_reactant[q]] * kh[q] / lw_pdiss[q];
      for (int b = 0; b < v.B; ++b) kabs[b] += v.lww[b] * f;
    }
}

// ---- PhotoModule: the host-side methods, over its own view --------------

inline void PhotoModule::sums(const double* Phi, PhotoSums& ps) const {
  ps.s.assign(NPHOTO, 0.0);
  ps.se.assign(NPHOTO, 0.0);
  ps.ss.assign(NPHOTO, 0.0);
  ps.sc.assign(NPHOTO, 0.0);
  ps.sp.assign(NPHOTO, 0.0);
  photo_sums(view(), Phi, ps.s.data(), ps.se.data(), ps.ss.data(),
             ps.sc.data(), ps.sp.data());
  ps.hab = habing(Phi);
}

inline double PhotoModule::habing(const double* Phi, const PhotoSums* ps)
    const {
  if (ps != nullptr) return ps->hab;
  return photo_habing(view(), Phi);
}

inline void PhotoModule::rates(const double* Phi, double* k, const double* sh,
                               const PhotoSums* ps) const {
  if (ps != nullptr) {
    const PhotoSumsView pv = ps->view();
    photo_rates(view(), Phi, k, sh, &pv);
  } else {
    photo_rates(view(), Phi, k, sh);
  }
}

inline double PhotoModule::heating(const double* Phi, const double* n,
                                   const double* sh, const PhotoSums* ps)
    const {
  if (ps != nullptr) {
    const PhotoSumsView pv = ps->view();
    return photo_heating(view(), Phi, n, sh, &pv);
  }
  return photo_heating(view(), Phi, n, sh);
}

#ifdef KROME_HAS_SECION
inline double PhotoModule::secondary(const double* Phi, const double* n,
                                     double* dn_add, const double* sh,
                                     const PhotoSums* ps) const {
  if (ps != nullptr) {
    const PhotoSumsView pv = ps->view();
    return photo_secondary(view(), Phi, n, dn_add, sh, &pv);
  }
  return photo_secondary(view(), Phi, n, dn_add, sh);
}
#endif

inline void PhotoModule::absorption(const double* n, const double* sh,
                                    const double* kh, double* kabs) const {
  photo_absorption(view(), n, sh, kh, kabs);
}

} // namespace krome
