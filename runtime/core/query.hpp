// KROME-X query / diagnostic functions (legacy krome_get_*), per cell.
// Implementation header: the functions are declared and documented in
// krome.hpp. Everything is evaluated exactly as the solver
// sees it: the temperature cap [Tcmb, 1e9 K], the photo field and binned
// dust of the SolverOpts, the legacy-helper hooks (n_global = n).
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace krome {

// the solver's view of one cell, for evaluations outside a solve: context
// with fixed Tgas (the query layout), the hooks (Params::hk) installed in
// its own copy of the parameters
struct QueryCell {
  SolveCtx c;
  int flag = 0;
  QueryCell(const double* n, double Tgas, const Params& p,
            const SolverOpts& o) {
    SolverOpts oq = o;
    oq.evolve_T = false;
    oq.conserve = 0;  // no solve: no post-solve correction
    flag = init_ctx_(c, Tgas, p, oq);
    // the queries evaluate n alone: with dust.evolve the grains keep
    // the population's radii (DustBins::a), not state slots
    c.ia = -1;
    c.A = 0;
    if (flag == 0) install_hooks_(c, n, o);
  }
  double T(const double* n) const { return state_T_(n, &c); }
};

// rate coefficients k[NREA] (generated + photo rates); 0 or a solve()
// input flag (-102/-103/-105)
inline int get_rates(const double* n, double Tgas, const Params& p,
                     const SolverOpts& o, double* k) {
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  double* kk = q.c.tc.k;
  eval_rates_(q.T(n), n, &q.c, kk);
  for (int r = 0; r < NREA; ++r) k[r] = kk[r];
  return 0;
}

// reaction fluxes k * prod(n) [cm^-3 s^-1], fl[NREA]
inline int get_fluxes(const double* n, double Tgas, const Params& p,
                      const SolverOpts& o, double* fl) {
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  eval_rates_(q.T(n), n, &q.c, q.c.tc.k);
  fluxes(n, q.c.tc.k, fl);
  return 0;
}

// the dust temperature of the rates: table_Tdust at the rate temperature
// (dust.tables, the Tdust parameter, or the binned-dust hook, installed
// and balanced by eval_rates_ as in a solve)
inline int get_Tdust(const double* n, double Tgas, const Params& p,
                     const SolverOpts& o, double* Td) {
#if defined(KROME_HAS_DUST_TABLES) || defined(KROME_HAS_PAR_Tdust)
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  const double T = q.T(n);
  eval_rates_(T, n, &q.c, q.c.tc.k);
  *Td = table_Tdust(T, n, q.c.tc.p);
  return 0;
#else
  (void)n; (void)Tgas; (void)p; (void)o; (void)Td;
  return -1;
#endif
}

// the thermal balance [erg cm^-3 s^-1]: totals, the active terms in the
// heat_term_names / cool_term_names order (heats[NHEAT], cools[NCOOL],
// each may be nullptr), and dTgas/dt [K s^-1]
inline int get_thermal(const double* n, double Tgas, const Params& p,
                       const SolverOpts& o, double* heat, double* cool,
                       double* heats, double* cools, double* dTdt) {
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  eval_rates_(q.T(n), n, &q.c, q.c.tc.k);
  // the thermal terms at the thermal temperature, as in a solve
  const double T = thermal_T_(Tgas);
  double H[NHEAT_ALL], C[NCOOL_ALL];
  thermal_terms_(T, n, q.c.tc.k, &q.c.tc, q.c.dnx, nullptr, H, C);
  double hs = 0.0, cs = 0.0;
  int ih = 0, ic = 0;
  for (int i = 0; i < NHEAT_ALL; ++i) {
    hs += H[i];
    if (heat_active_[i] && heats != nullptr) heats[ih] = H[i];
    ih += heat_active_[i] ? 1 : 0;
  }
  for (int i = 0; i < NCOOL_ALL; ++i) {
    cs += C[i];
    if (cool_active_[i] && cools != nullptr) cools[ic] = C[i];
    ic += cool_active_[i] ? 1 : 0;
  }
  if (heat != nullptr) *heat = hs;
  if (cool != nullptr) *cool = cs;
  if (dTdt != nullptr)
    *dTdt = (gamma_ad_(T, n, q.c.tc.p) - 1.0) * (hs - cs)
        / (KB_ERG * ntot_(n));
  return 0;
}

// names of the active thermal terms, in the order of get_thermal
inline const char* heating_name(int i) {
  for (int a = 0, m = 0; a < NHEAT_ALL; ++a)
    if (heat_active_[a] && m++ == i) return heat_term_names[a];
  return nullptr;
}
inline const char* cooling_name(int i) {
  for (int a = 0, m = 0; a < NCOOL_ALL; ++a)
    if (cool_active_[a] && m++ == i) return cool_term_names[a];
  return nullptr;
}

// adiabatic index at (n, Tgas) (thermal.gamma), temperature capped as the
// solver's thermal terms
inline double get_gamma(const double* n, double Tgas, const Params& p) {
  return gamma_ad_(thermal_T_(Tgas), n, p);
}

// mean molecular weight rho / (ntot m_p), legacy get_mu (floored sums)
inline double get_mu(const double* n) {
  return std::fmax(rho_(n), 1e-40)
      / (std::fmax(ntot_(n), 1e-40) * constants::p_mass);
}
inline double get_rho(const double* n) { return rho_(n); }
inline double get_ntot(const double* n) { return ntot_(n); }
inline double get_Hnuclei(const double* n) { return Hnuclei_(n); }

// net charge density sum q n [e cm^-3]; the electron density that makes the
// cell neutral, max(sum over ions q n, 0) (legacy get_electrons)
inline double get_charge(const double* n) {
  double q = 0.0;
  for (int i = 0; i < NSPEC; ++i) q += species_charge[i] * n[i];
  return q;
}
inline double get_electrons_neutral(const double* n) {
  double q = get_charge(n);
#ifdef KROME_HAS_E
  q += n[IDX_E];  // all species but the electrons themselves
#endif
  return std::fmax(q, 0.0);
}

// element totals [atoms cm^-3], tot[NELEM] (opaque species excluded)
inline void get_element_totals(const double* n, double* tot) {
  for (int e = 0; e < NELEM; ++e) tot[e] = 0.0;
  for (int i = 0; i < NSPEC; ++i)
    for (int e = 0; e < NELEM; ++e)
      tot[e] += species_elements[i * NELEM + e] * n[i];
}

// free-fall time sqrt(3 pi / (32 G rho)) [s] and Jeans length
// sqrt(pi kB T / (G rho mu m_p)) [cm] (legacy get_free_fall_time,
// get_jeans_length); get_free_fall_time_rho from the mass density rho
// [g cm^-3] (legacy krome_get_free_fall_time_rho)
inline double get_free_fall_time_rho(double rho) {
  return std::sqrt(3.0 * constants::pi / (32.0 * constants::gravity
                                          * std::fmax(rho, 1e-40)));
}
inline double get_free_fall_time(const double* n) {
  return get_free_fall_time_rho(rho_(n));
}

// the column density [cm^-2] the network's column estimator
// (radiation.num2col, legacy -columnDensityMethod) gives a number density num
// [cm^-3] of the cell (n, Tgas) (legacy krome_num2col)
inline double num2col(double num, const double* n, double Tgas,
                      const Params& p) {
  return num2col_(num, std::fmin(std::fmax(Tgas, p.Tcmb), 1e9), n, p);
}
inline double get_jeans_length(const double* n, double Tgas) {
  const double rr = std::fmax(rho_(n), 1e-40);
  return std::sqrt(constants::pi * constants::boltzmann_erg * Tgas
                   / (constants::gravity * rr * get_mu(n) * constants::p_mass));
}

#ifdef KROME_PHOTO_HPP_INCLUDED
// int_0^N f(N') dN' for a column-form shielding factor: f is flat below
// 1e-10 N, then composite 4-point Gauss-Legendre in ln N'
template <class F>
inline double col_integral_(F f, double N) {
  if (!(N > 0.0)) return 0.0;
  static const double x[4] = {-0.8611363115940526, -0.3399810435848563,
                              0.3399810435848563, 0.8611363115940526};
  static const double w[4] = {0.3478548451374538, 0.6521451548625461,
                              0.6521451548625461, 0.3478548451374538};
  const double N0 = 1e-10 * N;
  const double a = std::log(N0), L = std::log(N) - a;
  const int m = 96;
  const double h = L / m;
  double s = f(N0) * N0;
  for (int i = 0; i < m; ++i) {
    const double c = a + (i + 0.5) * h;
    for (int q = 0; q < 4; ++q) {
      const double Np = std::exp(c + 0.5 * h * x[q]);
      s += 0.5 * h * w[q] * f(Np) * Np;
    }
  }
  return s;
}

// the radiative-channel factors at (n, Tgas): shields, J_LW unit rates
struct RadChannels_ {
  double T;
  double sh[NPHOTO + 1], kh[NLW + 1];
  RadChannels_(const double* n, double Tgas, const Params& p)
      : T(std::fmin(std::fmax(Tgas, p.Tcmb), 1e9)) {
    Params pq = p;  // n_global = n for the shield expressions
    pq.hk.n_global = n;
    for (int r = 0; r <= NPHOTO; ++r) sh[r] = 1.0;
    for (int q = 0; q <= NLW; ++q) kh[q] = 0.0;
    photo_shields(T, n, pq, sh);
    lw_unit_rates(T, n, pq, kh);
  }
};

// optical depth per photon bin, tau[B], with the num2col columns (legacy
// krome_get_opacity; the length is the one num2col implies):
//   photo r   sigma_rb / pdiss_r  int_0^N s_r dN'   (N = num2col(n_r))
//   J_LW q    lww_b / pdiss_q     int_0^N khat_q dN'
// The integral is the curve of growth of a column-form shielding factor
// (one fselfH2* call, its column argument swept); otherwise s(N) N,
// which for an unshielded channel is sigma N. With d2g > 0 the dust adds
// kappa_b d2g Sigma_gas, Sigma_gas = N_H rho / n_Hnuclei (N_H =
// NH_col_, as the Habing dust shielding): the dust bins of o.dust when
// configured (pi a^2 Qabs per grain), else kdust[B] [cm^2 per g of dust]
// when given. Returns 0, -1 (no photon bins) or -105 (a REQUIRED
// parameter unset)
inline int get_opacity(const double* n, double Tgas, const Params& p,
                       const SolverOpts& o, double d2g, const double* kdust,
                       double* tau) {
  if (o.photo == nullptr) return -1;
  if (params_missing(p) >= 0) return -105;
  const PhotoModule& ph = *o.photo;
  const int B = ph.B;
  RadChannels_ rc(n, Tgas, p);
  const double T = rc.T;
  for (int b = 0; b < B; ++b) tau[b] = 0.0;
  for (int r = 0; r < NPHOTO; ++r) {
    const double N = num2col_(n[photo_reactant[r]], T, n, p);
    const double I = photo_shield_colform[r]
        ? col_integral_([&](double Np) {
            return photo_shield_col(r, Np, T, n, p); }, N)
        : rc.sh[r] * N;
    for (int b = 0; b < B; ++b)
      tau[b] += ph.sigma[static_cast<size_t>(r) * B + b] * I
          / photo_pdiss[r];
  }
  for (int q = 0; q < NLW; ++q) {
    const double N = num2col_(n[lw_reactant[q]], T, n, p);
    const double I = lw_colform[q]
        ? col_integral_([&](double Np) {
            return lw_unit_rate_col(q, Np, T, n, p); }, N)
        : rc.kh[q] * N;
    for (int b = 0; b < B; ++b) tau[b] += ph.lww[b] * I / lw_pdiss[q];
  }
  if (d2g > 0.0) {
    const double nH = Hnuclei_(n);
    const double Sg = nH > 0.0 ? NH_col_(T, n, p) * rho_(n) / nH : 0.0;
    bool bins = false;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
    if (o.dust != nullptr && !o.dust->qmid.empty()) {
      const DustBins& d = *o.dust;
      for (int i = 0; i < d.nb; ++i)
        for (int b = 0; b < B; ++b)
          tau[b] += DustBins::PI_ * d.asize2[i] * d.xw[i] * d2g * Sg
              * d.qmid[static_cast<size_t>(i) * B + b];
      bins = true;
    }
#endif
    if (!bins && kdust != nullptr)
      for (int b = 0; b < B; ++b) tau[b] += kdust[b] * d2g * Sg;
  }
  return 0;
}

// photon removal coefficient per photon bin, kabs[B] [cm^-1]: the
// photochemistry (sum_r s_r sigma_rb n_r / pdiss_r + the J_LW channel,
// exactly the radiation.evolve sink), the binned dust of o.dust (sum_i pi
// a_i^2 n_d,i Qabs_i at the bin midpoint, legacy dust_opacity) and, when
// kdust is given, rho_gas kdust_b with kdust_b the dust opacity per unit
// GAS mass [cm^2 g^-1] (e.g. d2g kappa_b; legacy
// krome_get_opacity_size_d2g). Times a length it is legacy's
// krome_get_opacity_size(_d2g). With radiation.evolve the solver already
// removed the photochemical part from the field it returns.
inline int get_absorption(const double* n, double Tgas, const Params& p,
                          const SolverOpts& o, const double* kdust,
                          double* kabs) {
  if (o.photo == nullptr) return -1;
  if (params_missing(p) >= 0) return -105;
  const int B = o.photo->B;
  RadChannels_ rc(n, Tgas, p);
  o.photo->absorption(n, rc.sh, rc.kh, kabs);
#if defined(KROME_DUST_BINS_HPP_INCLUDED) && defined(KROME_HAS_PAR_d2g)
  if (o.dust != nullptr && !o.dust->qmid.empty()) {
    const DustBins& d = *o.dust;
    const double xs = rho_(n) * p.d2g;
    for (int i = 0; i < d.nb; ++i)
      for (int b = 0; b < B; ++b)
        kabs[b] += DustBins::PI_ * d.asize2[i] * d.xw[i] * xs
            * d.qmid[static_cast<size_t>(i) * B + b];
  }
#endif
  if (kdust != nullptr) {
    const double rr = rho_(n);
    for (int b = 0; b < B; ++b) kabs[b] += rr * kdust[b];
  }
  return 0;
}
#endif

#ifdef KROME_PHOTO_HPP_INCLUDED
// photoheating per photo reaction, heat[NPHOTO] [erg cm^-3 s^-1] (legacy
// krome_get_photoBin_heats, which gave it per absorber: heat[r] /
// n[photo_reactant[r]]): each reaction's photoelectron deposit with its
// shielding, and with radiation.secondary_ionization the mixture's heat share of
// its fast part; they sum to the PHOTO heating term (0 for the Solomon
// J_LW-free lines, which deposit nothing). -1 without photon bins, -103
// without a field, -105 with a REQUIRED parameter unset
inline int get_photo_heating(const double* n, double Tgas, const Params& p,
                             const SolverOpts& o, double* heat) {
  if (o.photo == nullptr) return -1;
  if (o.Phi == nullptr) return -103;
  if (params_missing(p) >= 0) return -105;
  const PhotoView v = o.photo->view();
  RadChannels_ rc(n, Tgas, p);
  const double eV_to_erg = constants::eV_to_erg;
#ifdef KROME_HAS_SECION
  const double eta = photo_heat_eta_(n);
#endif
  for (int r = 0; r < NPHOTO; ++r) {
    heat[r] = 0.0;
    if (!photo_deposit[r]) continue;
    const double nr = n[photo_reactant[r]] * rc.sh[r];
#ifdef KROME_HAS_SECION
    double se, ss, sc, sp;
    photo_sums_secondary_(v, o.Phi, r, nullptr, se, ss, sc, sp);
    heat[r] = ((se - ss) + eta * ss) * nr * eV_to_erg;
#else
    double se = 0.0;
    for (int b = 0; b < v.B; ++b)
      se += v.sigmaE[static_cast<size_t>(r) * v.B + b] * o.Phi[b];
    heat[r] = se * nr * eV_to_erg;
#endif
  }
  return 0;
}
#endif

// the field over the Draine ISRF (legacy krome_get_ratioFluxDraine: the
// integrated flux / 7.95e-8, i.e. sum_b Phi_b h / (4 pi) / 7.95e-8, the
// rate-expression function get_ratioFluxDraine). -1 without photon bins
// or field
inline int get_ratio_flux_draine(const SolverOpts& o, double* ratio) {
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (o.photo == nullptr || o.Phi == nullptr) return -1;
  double s = 0.0;
  for (int b = 0; b < o.photo->B; ++b) s += o.Phi[b];
  *ratio = s * constants::planck_eV / (4.0 * constants::pi) / 7.95e-8;
  return 0;
#else
  (void)o; (void)ratio;
  return -1;
#endif
}

// the line coolants (the coolZ.dat metals of -cooling=Z, then the
// thermal.coolant_files coolants; legacy krome_cooling<metal>): each one's line
// cooling [erg cm^-3 s^-1], cool[NCOOLANT], at the thermal temperature
// as get_thermal (the Z term and the coolant-file terms are their sums,
// before the Z_EXTENDED merge weight)
inline const char* coolant_name(int i) {
  return (i >= 0 && i < NCOOLANT) ? coolant_names[i] : nullptr;
}
inline int get_coolant_cooling(const double* n, double Tgas, const Params& p,
                               const SolverOpts& o, double* cool) {
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  eval_rates_(q.T(n), n, &q.c, q.c.tc.k);  // the hooks' field, as a solve
  const double T = thermal_T_(Tgas);
  for (int i = 0; i < NCOOLANT; ++i)
    cool[i] = coolant_eval_(i, T, n, q.c.tc.p, nullptr);
  return 0;
}

// the level populations [cm^-3] of coolant i, pop[coolant_nlev[i]], in the
// file's level numbering (energies coolant_level_E[coolant_lev0[i] + l]
// [K], degeneracies coolant_level_g[...]); levels no transition connects
// to the ground have population 0. -1 for a bad index
inline int get_coolant_populations(const double* n, double Tgas,
                                   const Params& p, const SolverOpts& o,
                                   int i, double* pop) {
  if (i < 0 || i >= NCOOLANT || pop == nullptr) return -1;
  QueryCell q(n, Tgas, p, o);
  if (q.flag != 0) return q.flag;
  eval_rates_(q.T(n), n, &q.c, q.c.tc.k);
  coolant_eval_(i, thermal_T_(Tgas), n, q.c.tc.p, pop);
  return 0;
}

// energy released per reaction event, q[NREA] [MeV] (the ; qeff=
// annotations, 0 elsewhere; legacy krome_get_qeff). Times get_fluxes it is
// legacy krome_stars_energy (in MeV cm^-3 s^-1).
inline int get_qeff(const double* n, double Tgas, const Params& p,
                    double* q) {
  if (params_missing(p) >= 0) return -105;
  reaction_qeff(std::fmin(std::fmax(Tgas, p.Tcmb), 1e9), n, p, q);
  return 0;
}

// number densities <-> mass fractions x_i = m_i n_i / rho (legacy
// krome_n2x / krome_x2n); rho [g cm^-3]
inline void n2x(const double* n, double* x) {
  const double rr = std::fmax(rho_(n), 1e-40);
  for (int i = 0; i < NSPEC; ++i) x[i] = species_mass_g[i] * n[i] / rr;
}
inline void x2n(const double* x, double rho, double* n) {
  for (int i = 0; i < NSPEC; ++i)
    n[i] = species_mass_g[i] > 0.0 ? rho * x[i] / species_mass_g[i] : 0.0;
}

// ---- legacy conveniences -------------------------------------------------

inline std::string names_header() {
  std::string s;
  for (int i = 0; i < NSPEC; ++i) {
    if (i) s += ' ';
    s += species_names[i];
  }
  return s;
}

inline int print_best_flux(const double* n, double Tgas, const Params& p,
                           int nbest, const SolverOpts& o, std::FILE* f) {
  std::vector<double> fl(NREA);
  const int flag = get_fluxes(n, Tgas, p, o, fl.data());
  if (flag != 0) return flag;
  std::vector<int> idx(NREA);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&](int a, int b) { return fl[a] > fl[b]; });
  std::fprintf(f, " ***************\n");
  for (int i = 0; i < std::min(nbest, NREA); ++i)
    std::fprintf(f, "%8d %-50s %17.8E\n", idx[i] + 1,
                 reaction_names[idx[i]], fl[idx[i]]);
  return 0;
}

inline int explore_flux(std::FILE* f, const double* n, double Tgas,
                        const Params& p, double xvar, const SolverOpts& o) {
  std::vector<double> fl(NREA);
  const int flag = get_fluxes(n, Tgas, p, o, fl.data());
  if (flag != 0) return flag;
  const double fmax = NREA > 0 ? *std::max_element(fl.begin(), fl.end())
                               : 0.0;
  const double fsum = std::accumulate(fl.begin(), fl.end(), 0.0);
  for (int i = 0; i < NREA; ++i)
    std::fprintf(f, "%8d %17.8E %17.8E %17.8E %17.8E %17.8E   %s\n", i + 1,
                 xvar, Tgas, fl[i], fl[i] / fmax, fl[i] / fsum,
                 reaction_names[i]);
  std::fprintf(f, "\n");
  return 0;
}

// the level populations of every coolant (legacy krome_popcool_dump),
// one row per level: the coolant name, the level (1-based), Tgas, the
// population [cm^-3] and the coolant's total; a blank line at the end
inline int print_coolant_populations(std::FILE* f, const double* n,
                                     double Tgas, const Params& p,
                                     const SolverOpts& o) {
  for (int i = 0; i < NCOOLANT; ++i) {
    std::vector<double> pop(static_cast<size_t>(coolant_nlev[i]));
    const int flag = get_coolant_populations(n, Tgas, p, o, i, pop.data());
    if (flag != 0) return flag;
    const double tot = std::accumulate(pop.begin(), pop.end(), 0.0);
    for (int l = 0; l < coolant_nlev[i]; ++l)
      std::fprintf(f, "%8s%5d%17.8E%17.8E%17.8E\n", coolant_names[i], l + 1,
                   Tgas, pop[static_cast<size_t>(l)], tot);
  }
  std::fprintf(f, "\n");
  return 0;
}

// one row: Tgas, the total cooling and every active cooling term
// [erg cm^-3 s^-1] in the cooling_name order (legacy krome_dump_cooling)
inline int dump_cooling(std::FILE* f, const double* n, double Tgas,
                        const Params& p, const SolverOpts& o) {
  double cool = 0.0;
  double cs[NCOOL > 0 ? NCOOL : 1];
  const int flag = get_thermal(n, Tgas, p, o, nullptr, &cool, nullptr, cs,
                               nullptr);
  if (flag != 0) return flag;
  std::fprintf(f, "%14.5E%14.5E", Tgas, cool);
  for (int i = 0; i < NCOOL; ++i) std::fprintf(f, "%14.5E", cs[i]);
  std::fprintf(f, "\n");
  return 0;
}

// the rate coefficients on npts log-spaced temperatures Tmin..Tmax at
// n = 1e-40 cm^-3 for every species (legacy krome_dump_rates): rows Tgas,
// the reaction (1-based), k; a blank line after each temperature
inline int dump_rates(std::FILE* f, double Tmin, double Tmax, int npts,
                      const Params& p, const SolverOpts& o) {
  if (!(Tmin > 0.0 && Tmax >= Tmin) || npts < 1) return -1;
  double n[NSPEC];
  for (int i = 0; i < NSPEC; ++i) n[i] = 1e-40;
  std::vector<double> k(NREA > 0 ? NREA : 1);
  const double a = std::log10(Tmin), b = std::log10(Tmax);
  for (int i = 0; i < npts; ++i) {
    const double T = std::pow(10.0, npts > 1 ? a + i * (b - a) / (npts - 1)
                                             : a);
    const int flag = get_rates(n, T, p, o, k.data());
    if (flag != 0) return flag;
    for (int r = 0; r < NREA; ++r)
      std::fprintf(f, "%17.8E%8d%17.8E\n", T, r + 1,
                   k[static_cast<size_t>(r)]);
    std::fprintf(f, "\n");
  }
  return 0;
}

// every parameter of p, 'name: value' per line (legacy
// krome_print_phys_variables)
inline void print_params(const Params& p, std::FILE* f) {
  Params pc = p;
  for (int i = 0; i < NPAR; ++i)
    std::fprintf(f, "%s: %.16E\n", param_names[i], *param_ptr(pc, i));
}

inline void print_info(const double* n, double Tgas, std::FILE* f) {
  std::fprintf(f, "%4s %-9s%11s%5s%11s\n", "#", "Name", "m (g)", "Chrg",
               "x");
  double s = 0.0;
  for (int i = 0; i < NSPEC; ++i) {
    std::fprintf(f, "%4d %-9s%11.3E%5d%11.3E\n", i + 1, species_names[i],
                 species_mass_g[i], int(species_charge[i]), n[i]);
    s += n[i];
  }
  std::fprintf(f, "%30s%11.3E\n", " sum", s);
  std::fprintf(f, "%14s%11.3E\n", "Tgas", Tgas);
}

} // namespace krome
