// KROME-X C interface implementation (see krome.h). Compiled once per
// project against the generated network (krome.hpp includes
// krome_network.hpp, written by codegen). Link with net_<name>_data.cpp
// (no external library).
#include <cmath>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include "krome.hpp"  // the network and the C++ API
#include "krome.h"

// the C flag constants are the C++ ones
static_assert(KROME_FLAG_OK == krome::FLAG_OK, "flag OK");
static_assert(KROME_FLAG_TOO_MUCH_WORK == krome::FLAG_TOO_MUCH_WORK, "flag TOO_MUCH_WORK");
static_assert(KROME_FLAG_TOO_MUCH_ACC == krome::FLAG_TOO_MUCH_ACC, "flag TOO_MUCH_ACC");
static_assert(KROME_FLAG_ERR_FAILURE == krome::FLAG_ERR_FAILURE, "flag ERR_FAILURE");
static_assert(KROME_FLAG_CONV_FAILURE == krome::FLAG_CONV_FAILURE, "flag CONV_FAILURE");
static_assert(KROME_FLAG_LSETUP_FAIL == krome::FLAG_LSETUP_FAIL, "flag LSETUP_FAIL");
static_assert(KROME_FLAG_LSOLVE_FAIL == krome::FLAG_LSOLVE_FAIL, "flag LSOLVE_FAIL");
static_assert(KROME_FLAG_RHSFUNC_FAIL == krome::FLAG_RHSFUNC_FAIL, "flag RHSFUNC_FAIL");
static_assert(KROME_FLAG_FIRST_RHSFUNC_ERR == krome::FLAG_FIRST_RHSFUNC_ERR, "flag FIRST_RHSFUNC_ERR");
static_assert(KROME_FLAG_REPTD_RHSFUNC_ERR == krome::FLAG_REPTD_RHSFUNC_ERR, "flag REPTD_RHSFUNC_ERR");
static_assert(KROME_FLAG_UNREC_RHSFUNC_ERR == krome::FLAG_UNREC_RHSFUNC_ERR, "flag UNREC_RHSFUNC_ERR");
static_assert(KROME_FLAG_NLS_FAIL == krome::FLAG_NLS_FAIL, "flag NLS_FAIL");
static_assert(KROME_FLAG_ILL_INPUT == krome::FLAG_ILL_INPUT, "flag ILL_INPUT");
static_assert(KROME_FLAG_TOO_CLOSE == krome::FLAG_TOO_CLOSE, "flag TOO_CLOSE");
static_assert(KROME_FLAG_NO_DUST == krome::FLAG_NO_DUST, "flag NO_DUST");
static_assert(KROME_FLAG_NO_FIELD == krome::FLAG_NO_FIELD, "flag NO_FIELD");
static_assert(KROME_FLAG_UNSUPPORTED == krome::FLAG_UNSUPPORTED, "flag UNSUPPORTED");
static_assert(KROME_FLAG_PARAM_UNSET == krome::FLAG_PARAM_UNSET, "flag PARAM_UNSET");
static_assert(KROME_FLAG_CONSERVE_GAS_TOTAL == krome::FLAG_CONSERVE_GAS_TOTAL, "flag CONSERVE_GAS_TOTAL");
static_assert(KROME_FLAG_SURFACE_BINS == krome::FLAG_SURFACE_BINS, "flag SURFACE_BINS");
static_assert(KROME_FLAG_NO_LAPACK == krome::FLAG_NO_LAPACK, "flag NO_LAPACK");
static_assert(KROME_FLAG_NO_KOKKOS == krome::FLAG_NO_KOKKOS, "flag NO_KOKKOS");
static_assert(KROME_FLAG_NONPOSITIVE_T == krome::FLAG_NONPOSITIVE_T,
              "flag NONPOSITIVE_T");
static_assert(KROME_FLAG_DEVICE_MEMORY == krome::FLAG_DEVICE_MEMORY,
              "flag DEVICE_MEMORY");
static_assert(KROME_FLAG_CONSERVE_UNREACHED == krome::FLAG_CONSERVE_UNREACHED,
              "flag CONSERVE_UNREACHED");
static_assert(KROME_FLAG_CONSERVE_LARGE == krome::FLAG_CONSERVE_LARGE,
              "flag CONSERVE_LARGE");
static_assert(KROME_FLAG_HOST_RECOVERED == krome::FLAG_HOST_RECOVERED,
              "flag HOST_RECOVERED");

struct krome_opts {
  krome::SolverOpts so;
  krome::CellsOpts co;
  std::string error_file = "KROME_ERROR_REPORT";  // so.error_file points here
  krome::EquilOpts eq;
  bool has_photo = false;
  double broad_T = 0.0, broad_vturb = 0.0;  // krome_opts_set_photo_broadening
  krome::Radiation rad;
  krome::PhotoModule photo;
  std::vector<double> kdraine;  // Draine R_V=3.1 kappa per bin (get_opacity)
  bool has_dust = false;
  krome::DustBins dust;
  std::vector<double> atol_sp, rtol_sp;  // so.atol_species / rtol_species
};

namespace {

// value slot of a key (nullptr if unknown); booleans go through `flag`
bool* bool_slot(krome_opts* h, const char* key) {
  if (!std::strcmp(key, "evolve_T")) return &h->so.evolve_T;
  if (!std::strcmp(key, "clip_negative")) return &h->so.clip_negative;
  if (!std::strcmp(key, "error_report")) return &h->so.error_report;
  if (!std::strcmp(key, "renormalize_x")) return &h->so.renormalize_x;
  return nullptr;
}

int* int_slot(krome_opts* h, const char* key) {
  if (!std::strcmp(key, "max_restarts")) return &h->so.max_restarts;
  if (!std::strcmp(key, "max_err_fails")) return &h->so.max_err_fails;
  if (!std::strcmp(key, "max_conv_fails")) return &h->so.max_conv_fails;
  if (!std::strcmp(key, "conserve")) return &h->so.conserve;
  return nullptr;
}

double* double_slot(krome_opts* h, const char* key) {
  if (!std::strcmp(key, "eq_t0")) return &h->eq.t0;
  if (!std::strcmp(key, "eq_tmax")) return &h->eq.t_max;
  if (!std::strcmp(key, "eq_factor")) return &h->eq.factor;
  if (!std::strcmp(key, "eq_rtol")) return &h->eq.ertol;
  if (!std::strcmp(key, "eq_atol")) return &h->eq.eatol;
  if (!std::strcmp(key, "rtol")) return &h->so.rtol;
  if (!std::strcmp(key, "atol_n")) return &h->so.atol_n;
  if (!std::strcmp(key, "atol_T")) return &h->so.atol_T;
  if (!std::strcmp(key, "atol_u")) return &h->so.atol_u;
  return nullptr;
}

void link_dust_photo(krome_opts* h) {
  if (h->has_dust && h->has_photo) h->dust.interp_qabs(h->rad);
}

// per-cell evaluation with the krome_solve_cells input conventions: calls
// f(c, n_c, Tgas_c, params_c, opts_c) for every cell; returns 0, -1 (bad
// arguments), -103 (photo bins configured, no Phi) or the first nonzero f
template <class F>
int for_cells_(int ncell, const double* n, const double* Tgas,
               const double* par, int par_stride, const double* Phi,
               const krome_opts* h, F f, bool need_phi = true) {
  if (ncell < 0 || n == nullptr) return -1;
  if (par != nullptr && par_stride != 0 && par_stride != krome::NPAR)
    return -1;
  krome::SolverOpts so = h != nullptr ? h->so : krome::SolverOpts();
  int B = 0;
  if (h != nullptr && h->has_photo) {
    if (Phi == nullptr && need_phi) return krome::FLAG_NO_FIELD;
    so.photo = &h->photo;
    B = Phi != nullptr ? h->photo.B : 0;
  }
  if (h != nullptr && h->has_dust)
    so.dust = const_cast<krome::DustBins*>(&h->dust);
  int first = 0;
  for (int c = 0; c < ncell; ++c) {
    krome::Params P;
    if (par != nullptr)
      krome::params_from_row(par + static_cast<size_t>(c) * par_stride, P);
    if (B > 0)  // read only: queries never evolve the field
      so.Phi = const_cast<double*>(Phi + static_cast<size_t>(c) * B);
    const int fl = f(c, n + static_cast<size_t>(c) * krome::NSPEC,
                     Tgas != nullptr ? Tgas[c] : 0.0, P, so);
    if (fl != 0 && first == 0) first = fl;
  }
  return first;
}

// f(file, n, Tgas, params, opts) on one cell, appending to path
template <class F>
int append_one_(const char* path, const double* n, double Tgas,
                const double* par, const double* Phi, const krome_opts* h,
                F f) {
  if (path == nullptr) return -1;
  std::FILE* fp = std::fopen(path, "a");
  if (fp == nullptr) return -1;
  const int r = for_cells_(1, n, &Tgas, par, 0, Phi, h,
                           [&](int, const double* nc, double T,
                               const krome::Params& P,
                               const krome::SolverOpts& o) {
                             return f(fp, nc, T, P, o);
                           });
  std::fclose(fp);
  return r;
}

}  // namespace

extern "C" {

const char* krome_network_name(void) { return krome::NETWORK_NAME; }
const char* krome_version(void) { return krome::KROME_VERSION; }

int krome_nspec(void) { return krome::NSPEC; }
int krome_nrea(void) { return krome::NREA; }
int krome_npar(void) { return krome::NPAR; }

const char* krome_species_name(int i) {
  return (i >= 0 && i < krome::NSPEC) ? krome::species_names[i] : nullptr;
}

int krome_species_index(const char* name) {
  for (int i = 0; i < krome::NSPEC; ++i)
    if (!std::strcmp(krome::species_names[i], name)) return i;
  return -1;
}

const char* krome_param_name(int i) {
  return (i >= 0 && i < krome::NPAR) ? krome::param_names[i] : nullptr;
}

int krome_param_index(const char* name) {
  for (int i = 0; i < krome::NPAR; ++i)
    if (!std::strcmp(krome::param_names[i], name)) return i;
  return -1;
}

const char* krome_param_required_by(int i) {
  return (i >= 0 && i < krome::NPAR) ? krome::param_required_by[i]
                                     : nullptr;
}

void krome_params_default(double* par) {
  krome::Params p;
  for (int i = 0; i < krome::NPAR; ++i) par[i] = *krome::param_ptr(p, i);
}

krome_opts* krome_opts_new(void) { return new krome_opts(); }

void krome_opts_free(krome_opts* h) { delete h; }

int krome_opts_set(krome_opts* h, const char* key, double value) {
  if (h == nullptr || key == nullptr) return -1;
  if (bool* b = bool_slot(h, key)) {
    *b = value != 0.0;
    return 0;
  }
  if (double* d = double_slot(h, key)) {
    *d = value;
    return 0;
  }
  if (int* iv = int_slot(h, key)) {
    if (value < 0.0) return -1;
    if (iv == &h->so.conserve && value > 3.0) return -1;  // modes 0-3
    *iv = static_cast<int>(value);
    return 0;
  }
  if (!std::strcmp(key, "max_steps")) {
    h->so.max_steps = static_cast<long>(value);
    return 0;
  }
  if (!std::strcmp(key, "jac")) {
    const int v = static_cast<int>(value);
    if (v < 0 || v > 2) return -1;
    h->so.jac = v == 0 ? krome::JAC_AUTO
              : (v == 1 ? krome::JAC_ANALYTIC : krome::JAC_FD);
    return 0;
  }
  if (!std::strcmp(key, "linsol")) {
    const int v = static_cast<int>(value);
    if (v < 0 || v > 3) return -1;
    const krome::LinSol ls[4] = {krome::LINSOL_AUTO, krome::LINSOL_DENSE,
                                 krome::LINSOL_LAPACK, krome::LINSOL_SPARSE};
    h->so.linsol = ls[v];
    return 0;
  }
  if (!std::strcmp(key, "path")) {
    const int v = static_cast<int>(value);
    if (v < 0 || v > 2) return -1;
    h->co.path = v == 0 ? krome::PATH_AUTO
               : v == 1 ? krome::PATH_PERCELL : krome::PATH_DEVICE;
    return 0;
  }
  if (!std::strcmp(key, "device_cache")) {
    h->co.device_cache = value != 0.0;
    return 0;
  }
  if (!std::strcmp(key, "device_memory")) {
    h->co.device_memory = static_cast<long long>(value);
    return 0;
  }
  if (!std::strcmp(key, "min_interleave")) {
    h->co.min_interleave = static_cast<int>(value);
    return 0;
  }
  if (!std::strcmp(key, "min_device")) {
    h->co.min_device = static_cast<int>(value);
    return 0;
  }
  return -1;
}

int krome_opts_get(const krome_opts* hc, const char* key, double* value) {
  if (hc == nullptr || key == nullptr || value == nullptr) return -1;
  auto* h = const_cast<krome_opts*>(hc);
  if (bool* b = bool_slot(h, key)) {
    *value = *b ? 1.0 : 0.0;
    return 0;
  }
  if (double* d = double_slot(h, key)) {
    *value = *d;
    return 0;
  }
  if (int* iv = int_slot(h, key)) {
    *value = static_cast<double>(*iv);
    return 0;
  }
  if (!std::strcmp(key, "max_steps")) {
    *value = static_cast<double>(h->so.max_steps);
    return 0;
  }
  if (!std::strcmp(key, "jac")) {
    *value = h->so.jac == krome::JAC_AUTO ? 0.0
           : (h->so.jac == krome::JAC_ANALYTIC ? 1.0 : 2.0);
    return 0;
  }
  if (!std::strcmp(key, "linsol")) {
    *value = h->so.linsol == krome::LINSOL_AUTO ? 0.0
           : h->so.linsol == krome::LINSOL_DENSE ? 1.0
           : h->so.linsol == krome::LINSOL_LAPACK ? 2.0 : 3.0;
    return 0;
  }
  if (!std::strcmp(key, "path")) {
    *value = h->co.path == krome::PATH_AUTO ? 0.0
           : h->co.path == krome::PATH_PERCELL ? 1.0 : 2.0;
    return 0;
  }
  if (!std::strcmp(key, "min_device")) {
    *value = static_cast<double>(h->co.min_device);
    return 0;
  }
  if (!std::strcmp(key, "device_cache")) {
    *value = h->co.device_cache ? 1.0 : 0.0;
    return 0;
  }
  if (!std::strcmp(key, "device_memory")) {
    *value = static_cast<double>(h->co.device_memory);
    return 0;
  }
  if (!std::strcmp(key, "min_interleave")) {
    *value = static_cast<double>(h->co.min_interleave);
    return 0;
  }
  return -1;
}

int krome_opts_set_photo(krome_opts* h, int nbins, const double* Elo,
                         const double* Ehi) {
  if (h == nullptr || nbins <= 0 || Elo == nullptr || Ehi == nullptr)
    return -1;
  h->rad = krome::Radiation();
  for (int b = 0; b < nbins; ++b) {
    if (!(Ehi[b] > Elo[b])) return -1;
    h->rad.Elo.push_back(Elo[b]);
    h->rad.Ehi.push_back(Ehi[b]);
    h->rad.Phi.push_back(0.0);
  }
  h->rad.broad_T = h->broad_T;
  h->rad.broad_vturb = h->broad_vturb;
  h->photo = krome::PhotoModule();
  h->photo.init(h->rad);
  h->has_photo = true;
  h->kdraine.assign(nbins, 0.0);
  if (krome::spectra::dust_kappa_draine(nbins, Elo, Ehi, h->kdraine.data())
      != 0)
    h->kdraine.clear();  // bins outside the table: no default dust opacity
  link_dust_photo(h);
  return 0;
}

int krome_opts_set_photo_broadening(krome_opts* h, double Tgas,
                                    double vturb) {
  if (h == nullptr || !(Tgas >= 0.0) || !(vturb >= 0.0)) return -1;
  h->broad_T = Tgas;
  h->broad_vturb = vturb;
  if (!h->has_photo) return 0;  // applied by krome_opts_set_photo
  std::vector<double> Elo = h->rad.Elo, Ehi = h->rad.Ehi;
  return krome_opts_set_photo(h, static_cast<int>(Elo.size()), Elo.data(),
                              Ehi.data());
}

int krome_opts_nbins(const krome_opts* h) {
  return (h != nullptr && h->has_photo) ? h->photo.B : 0;
}

int krome_opts_set_dustbins(krome_opts* h, int nd, int use_Si, int use_C,
                            double alow, double aup, double phi) {
  if (h == nullptr || nd <= 0 || !(use_Si || use_C) || !(aup > alow))
    return -1;
  h->dust = krome::DustBins();
  h->dust.nd = nd;
  h->dust.useSi = use_Si != 0;
  h->dust.useC = use_C != 0;
  h->dust.alow = alow;
  h->dust.aup = aup;
  h->dust.phi = phi;
  h->dust.init();
  h->has_dust = true;
  link_dust_photo(h);
  return 0;
}

int krome_solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                      const double* par, int par_stride,
                      const krome_cells* extra, int* flags,
                      const krome_opts* hc) {
  double* Phi = extra != nullptr ? extra->Phi : nullptr;
  double* Td = extra != nullptr ? extra->Td : nullptr;
  double* adust = extra != nullptr ? extra->adust : nullptr;
  if (hc == nullptr || ncell < 0 || n == nullptr || Tgas == nullptr ||
      dt == nullptr)
    return -1;
  if (par != nullptr && par_stride != 0 && par_stride != krome::NPAR)
    return -1;
  if (ncell == 0) return 0;
  // the solve works on a copy of the options, so one handle can serve
  // concurrent calls (the handle itself is only read)
  krome::SolverOpts so = hc->so;
  so.error_file = hc->error_file.c_str();
  krome::CellsOpts co = hc->co;
  co.flags = flags;
  if (hc->has_photo) {
    // bins configured, no field given
    if (Phi == nullptr) return krome::FLAG_NO_FIELD;
    so.photo = &hc->photo;
    co.Phi = Phi;
  }
  if (hc->has_dust) {
    so.dust = const_cast<krome::DustBins*>(&hc->dust);  // copied per cell
    co.Td = Td;
    co.adust = adust;
  }

  // runtime parameters: dense rows in param_names order
  const bool per_cell = par != nullptr && par_stride != 0;
  std::vector<krome::Params> P(per_cell ? ncell : 1);
  if (par != nullptr)
    for (size_t c = 0; c < P.size(); ++c)
      krome::params_from_row(par + c * par_stride, P[c]);
  return per_cell
      ? krome::solve_cells(ncell, n, Tgas, dt, P.data(), so, co)
      : krome::solve_cells(ncell, n, Tgas, dt, P[0], so, co);
}

int krome_solve_cells_x(int ncell, double* x, const double* rho,
                        double* Tgas, const double* dt, const double* par,
                        int par_stride, const krome_cells* extra,
                        int* flags, const krome_opts* hc) {
  // the arguments checked here, as krome_solve_cells does: its -1 after
  // the solve is then FLAG_TOO_MUCH_WORK (a cell out of steps), not a
  // rejection, and the other cells' fractions must still be written
  if (hc == nullptr || ncell < 0 || x == nullptr || rho == nullptr
      || Tgas == nullptr || dt == nullptr)
    return -1;
  if (par != nullptr && par_stride != 0 && par_stride != krome::NPAR)
    return -1;
  if (ncell == 0) return 0;
  std::vector<double> n(static_cast<size_t>(ncell) * krome::NSPEC);
  std::vector<double> xin(static_cast<size_t>(ncell));
  for (int c = 0; c < ncell; ++c) {
    const double* xc = x + static_cast<size_t>(c) * krome::NSPEC;
    xin[c] = krome::x_sum_(xc);
    krome::x_to_n_(xc, rho[c], n.data() + static_cast<size_t>(c) * krome::NSPEC);
  }
  std::vector<int> fl(static_cast<size_t>(ncell), 0);
  int* f = flags != nullptr ? flags : fl.data();
  const int r = krome_solve_cells(ncell, n.data(), Tgas, dt, par, par_stride,
                                  extra, f, hc);
  // rejected before solving (flags not set)
  if (r == krome::FLAG_NO_FIELD) return r;
  for (int c = 0; c < ncell; ++c)
    if (f[c] == 0)  // with renormalize_x rescaled to the input sum
      krome::n_to_x_out_(n.data() + static_cast<size_t>(c) * krome::NSPEC,
                         rho[c], xin[c], hc->so,
                         x + static_cast<size_t>(c) * krome::NSPEC);
  return r;
}

const char* krome_flag_message(int flag) {
  return krome::flag_message(flag);
}

int krome_opts_ndust(const krome_opts* h) {
  return h != nullptr && h->has_dust ? h->dust.nb : 0;
}

int krome_opts_set_dust_Td(krome_opts* h, const double* Td) {
  if (h == nullptr || !h->has_dust || Td == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b) {
    if (!(Td[b] > 0.0)) return -1;
  }
  for (int b = 0; b < h->dust.nb; ++b) h->dust.Td[b] = Td[b];
  return 0;
}

int krome_opts_get_dust_Td(const krome_opts* h, double* Td) {
  if (h == nullptr || !h->has_dust || Td == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b) Td[b] = h->dust.Td[b];
  return 0;
}

namespace {
// krome_opts_set_atol_species / rtol_species: copy v[NSPEC] into the
// handle (NULL: the scalar for every species); NaN, and a zero atol, are
// rejected
int set_species_tol_(krome_opts* h, const double* v, bool atol) {
  if (h == nullptr) return -1;
  std::vector<double>& buf = atol ? h->atol_sp : h->rtol_sp;
  const double*& slot = atol ? h->so.atol_species : h->so.rtol_species;
  if (v == nullptr) {
    slot = nullptr;
    return 0;
  }
  for (int i = 0; i < krome::NSPEC; ++i)
    if (std::isnan(v[i]) || (atol && v[i] == 0.0)) return -1;
  buf.assign(v, v + krome::NSPEC);
  slot = buf.data();
  return 0;
}

int get_species_tol_(const krome_opts* h, double* v, bool atol) {
  if (h == nullptr || v == nullptr) return -1;
  const double* a = atol ? h->so.atol_species : h->so.rtol_species;
  const double s = atol ? h->so.atol_n : h->so.rtol;
  for (int i = 0; i < krome::NSPEC; ++i)
    v[i] = krome::species_tol_(a, i, s);
  return 0;
}
}  // namespace

int krome_opts_set_atol_species(krome_opts* h, const double* atol) {
  return set_species_tol_(h, atol, true);
}

int krome_opts_set_rtol_species(krome_opts* h, const double* rtol) {
  return set_species_tol_(h, rtol, false);
}

int krome_opts_get_atol_species(const krome_opts* h, double* atol) {
  return get_species_tol_(h, atol, true);
}

int krome_opts_get_rtol_species(const krome_opts* h, double* rtol) {
  return get_species_tol_(h, rtol, false);
}

int krome_opts_set_dust_a(krome_opts* h, const double* a) {
  if (h == nullptr || !h->has_dust || a == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b)
    if (!(a[b] >= 0.0)) return -1;
  for (int b = 0; b < h->dust.nb; ++b) h->dust.a[b] = a[b];
  h->dust.set_radii(nullptr);
  return 0;
}

int krome_opts_get_dust_a(const krome_opts* h, double* a) {
  if (h == nullptr || !h->has_dust || a == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b) a[b] = h->dust.a[b];
  return 0;
}

int krome_opts_get_dust_counts(const krome_opts* h, double rho_gas,
                               double d2g, double* nd) {
  if (h == nullptr || !h->has_dust || nd == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b) nd[b] = h->dust.xw[b] * rho_gas * d2g;
  return 0;
}

int krome_opts_get_dust_types(const krome_opts* h, int* isC) {
  if (h == nullptr || !h->has_dust || isC == nullptr) return -1;
  for (int b = 0; b < h->dust.nb; ++b) isC[b] = h->dust.isC[b];
  return 0;
}

int krome_opts_set_error_file(krome_opts* h, const char* path) {
  if (h == nullptr || path == nullptr) return -1;
  h->error_file = path;
  return 0;
}

int krome_equilibrium_cells(int ncell, double* n, double* Tgas,
                            const double* par, int par_stride,
                            const krome_cells* extra, int* flags,
                            const krome_opts* h) {
  if (h == nullptr || Tgas == nullptr) return -1;
  double* Phi = extra != nullptr ? extra->Phi : nullptr;
  double* Td = extra != nullptr ? extra->Td : nullptr;
  double* adust = extra != nullptr ? extra->adust : nullptr;
  // per-cell states are written in place: copy them out of the const view
  // the evaluation helper gives
  std::vector<int> fl(static_cast<size_t>(std::max(ncell, 0)), 0);
  const int r = for_cells_(
      ncell, n, Tgas, par, par_stride, Phi, h,
      [&](int c, const double*, double, const krome::Params& P,
          const krome::SolverOpts& o) {
        double* nc = n + static_cast<size_t>(c) * krome::NSPEC;
        krome::SolverOpts oc = o;
        krome::DustBins dc;  // each cell's own bins (Td is per-cell state)
        double* tdc = nullptr;
        double* adc = nullptr;
        if (o.dust != nullptr) {
          dc = *o.dust;
          oc.dust = &dc;
          if (Td != nullptr) {
            tdc = Td + static_cast<size_t>(c) * dc.nb;
            for (int b = 0; b < dc.nb; ++b) dc.Td[b] = tdc[b];
          }
          if (adust != nullptr) {
            adc = adust + static_cast<size_t>(c) * dc.nb;
            for (int b = 0; b < dc.nb; ++b) dc.a[b] = adc[b];
            dc.set_radii(nullptr);
          }
        }
        fl[c] = krome::solve_equilibrium(nc, Tgas[c], P, oc, h->eq);
        if (tdc != nullptr && fl[c] >= 0)
          for (int b = 0; b < dc.nb; ++b) tdc[b] = dc.Td[b];
        if (adc != nullptr && fl[c] >= 0)
          for (int b = 0; b < dc.nb; ++b) adc[b] = dc.a[b];
        return fl[c];
      });
  if (flags != nullptr)
    for (int c = 0; c < ncell; ++c) flags[c] = fl[c];
  return r;
}

/* ---- query / diagnostic functions ------------------------------------ */

int krome_get_rates(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, const double* Phi,
                    const krome_opts* h, double* k) {
  if (Tgas == nullptr || k == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_rates(
                          nc, T, P, o, k + static_cast<size_t>(c) * krome::NREA);
                    });
}

int krome_get_fluxes(int ncell, const double* n, const double* Tgas,
                     const double* par, int par_stride, const double* Phi,
                     const krome_opts* h, double* fl) {
  if (Tgas == nullptr || fl == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_fluxes(
                          nc, T, P, o,
                          fl + static_cast<size_t>(c) * krome::NREA);
                    });
}

int krome_get_Tdust(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, const double* Phi,
                    const krome_opts* h, double* Td) {
  if (Tgas == nullptr || Td == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_Tdust(nc, T, P, o, Td + c);
                    });
}

const char* krome_names_header(void) {
  static const std::string s = krome::names_header();
  return s.c_str();
}

int krome_print_best_flux(const double* n, double Tgas, const double* par,
                          const double* Phi, const krome_opts* h,
                          int nbest) {
  const int r = for_cells_(1, n, &Tgas, par, 0, Phi, h,
                           [&](int, const double* nc, double T,
                               const krome::Params& P,
                               const krome::SolverOpts& o) {
                             return krome::print_best_flux(nc, T, P, nbest,
                                                           o, stdout);
                           });
  std::fflush(stdout);
  return r;
}

int krome_explore_flux(const char* path, const double* n, double Tgas,
                       const double* par, const double* Phi,
                       const krome_opts* h, double xvar) {
  if (path == nullptr) return -1;
  std::FILE* f = std::fopen(path, "a");
  if (f == nullptr) return -1;
  const int r = for_cells_(1, n, &Tgas, par, 0, Phi, h,
                           [&](int, const double* nc, double T,
                               const krome::Params& P,
                               const krome::SolverOpts& o) {
                             return krome::explore_flux(f, nc, T, P, xvar,
                                                        o);
                           });
  std::fclose(f);
  return r;
}


int krome_print_coolant_populations(const char* path, const double* n,
                                    double Tgas, const double* par,
                                    const double* Phi, const krome_opts* h) {
  return append_one_(path, n, Tgas, par, Phi, h,
                     [](std::FILE* f, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                       return krome::print_coolant_populations(f, nc, T, P,
                                                               o);
                     });
}

int krome_dump_cooling(const char* path, const double* n, double Tgas,
                       const double* par, const double* Phi,
                       const krome_opts* h) {
  return append_one_(path, n, Tgas, par, Phi, h,
                     [](std::FILE* f, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                       return krome::dump_cooling(f, nc, T, P, o);
                     });
}

int krome_dump_rates(const char* path, double Tmin, double Tmax, int npts,
                     const double* par, const double* Phi,
                     const krome_opts* h) {
  double n[krome::NSPEC];
  for (double& v : n) v = 1e-40;
  return append_one_(path, n, Tmin, par, Phi, h,
                     [&](std::FILE* f, const double*, double,
                         const krome::Params& P, const krome::SolverOpts& o) {
                       return krome::dump_rates(f, Tmin, Tmax, npts, P, o);
                     });
}

void krome_print_params(const double* par) {
  krome::Params P;
  if (par != nullptr)
    krome::params_from_row(par, P);
  krome::print_params(P, stdout);
  std::fflush(stdout);
}

void krome_print_info(const double* n, double Tgas) {
  if (n == nullptr) return;
  krome::print_info(n, Tgas, stdout);
  std::fflush(stdout);
}

int krome_scale_Z(int ncell, double* n, double Z) {
  if (ncell < 0 || n == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    krome::scale_Z(n + static_cast<size_t>(c) * krome::NSPEC, Z);
  return 0;
}

int krome_nheat(void) { return krome::NHEAT; }
int krome_ncool(void) { return krome::NCOOL; }
const char* krome_heating_name(int i) { return krome::heating_name(i); }
const char* krome_cooling_name(int i) { return krome::cooling_name(i); }

int krome_get_thermal(int ncell, const double* n, const double* Tgas,
                      const double* par, int par_stride, const double* Phi,
                      const krome_opts* h, double* heat, double* cool,
                      double* heats, double* cools, double* dTdt) {
  if (Tgas == nullptr) return -1;
  return for_cells_(
      ncell, n, Tgas, par, par_stride, Phi, h,
      [&](int c, const double* nc, double T, const krome::Params& P,
          const krome::SolverOpts& o) {
        const size_t cc = static_cast<size_t>(c);
        return krome::get_thermal(
            nc, T, P, o, heat ? heat + cc : nullptr,
            cool ? cool + cc : nullptr,
            heats ? heats + cc * krome::NHEAT : nullptr,
            cools ? cools + cc * krome::NCOOL : nullptr,
            dTdt ? dTdt + cc : nullptr);
      });
}

int krome_get_qeff(int ncell, const double* n, const double* Tgas,
                   const double* par, int par_stride, double* q) {
  if (Tgas == nullptr || q == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, nullptr, nullptr,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts&) {
                      return krome::get_qeff(
                          nc, T, P, q + static_cast<size_t>(c) * krome::NREA);
                    });
}

int krome_get_gamma(int ncell, const double* n, const double* Tgas,
                    const double* par, int par_stride, double* gamma) {
  if (Tgas == nullptr || gamma == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, nullptr, nullptr,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts&) {
                      if (krome::params_missing(P) >= 0)
                        return krome::FLAG_PARAM_UNSET;
                      gamma[c] = krome::get_gamma(nc, T, P);
                      return 0;
                    });
}

#define KROME_CELLWISE_(NAME, EXPR)                                        \
  int NAME(int ncell, const double* n, double* out) {                      \
    if (ncell < 0 || n == nullptr || out == nullptr) return -1;            \
    for (int c = 0; c < ncell; ++c) {                                      \
      const double* nc = n + static_cast<size_t>(c) * krome::NSPEC;        \
      out[c] = (EXPR);                                                     \
    }                                                                      \
    return 0;                                                              \
  }
KROME_CELLWISE_(krome_get_mu, krome::get_mu(nc))
KROME_CELLWISE_(krome_get_rho, krome::get_rho(nc))
KROME_CELLWISE_(krome_get_ntot, krome::get_ntot(nc))
KROME_CELLWISE_(krome_get_Hnuclei, krome::get_Hnuclei(nc))
KROME_CELLWISE_(krome_get_charge, krome::get_charge(nc))
KROME_CELLWISE_(krome_get_electrons, krome::get_electrons_neutral(nc))
KROME_CELLWISE_(krome_get_free_fall_time, krome::get_free_fall_time(nc))
#undef KROME_CELLWISE_

int krome_get_jeans_length(int ncell, const double* n, const double* Tgas,
                           double* out) {
  if (ncell < 0 || n == nullptr || Tgas == nullptr || out == nullptr)
    return -1;
  for (int c = 0; c < ncell; ++c)
    out[c] = krome::get_jeans_length(
        n + static_cast<size_t>(c) * krome::NSPEC, Tgas[c]);
  return 0;
}

int krome_get_free_fall_time_rho(int ncell, const double* rho, double* out) {
  if (ncell < 0 || rho == nullptr || out == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    out[c] = krome::get_free_fall_time_rho(rho[c]);
  return 0;
}

int krome_num2col(int ncell, const double* num, const double* n,
                  const double* Tgas, const double* par, int par_stride,
                  double* out) {
  if (num == nullptr || Tgas == nullptr || out == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, nullptr, nullptr,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts&) {
                      out[c] = krome::num2col(num[c], nc, T, P);
                      return 0;
                    });
}

int krome_ncoolant(void) { return krome::NCOOLANT; }
const char* krome_coolant_name(int i) { return krome::coolant_name(i); }
int krome_coolant_species(int i) {
  return (i >= 0 && i < krome::NCOOLANT) ? krome::coolant_species[i] : -1;
}
int krome_coolant_nlev(int i) {
  return (i >= 0 && i < krome::NCOOLANT) ? krome::coolant_nlev[i] : -1;
}
int krome_coolant_levels(int i, double* E, double* g) {
  if (i < 0 || i >= krome::NCOOLANT) return -1;
  for (int l = 0; l < krome::coolant_nlev[i]; ++l) {
    if (E != nullptr) E[l] = krome::coolant_level_E[krome::coolant_lev0[i] + l];
    if (g != nullptr) g[l] = krome::coolant_level_g[krome::coolant_lev0[i] + l];
  }
  return 0;
}

int krome_get_coolant_cooling(int ncell, const double* n, const double* Tgas,
                              const double* par, int par_stride,
                              const double* Phi, const krome_opts* h,
                              double* cool) {
  if (Tgas == nullptr || cool == nullptr) return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_coolant_cooling(
                          nc, T, P, o,
                          cool + static_cast<size_t>(c) * krome::NCOOLANT);
                    });
}

int krome_get_coolant_populations(int ncell, const double* n,
                                  const double* Tgas, const double* par,
                                  int par_stride, const double* Phi,
                                  const krome_opts* h, int i, double* pop) {
  if (Tgas == nullptr || pop == nullptr || i < 0 || i >= krome::NCOOLANT)
    return -1;
  const size_t L = static_cast<size_t>(krome::coolant_nlev[i]);
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_coolant_populations(
                          nc, T, P, o, i, pop + static_cast<size_t>(c) * L);
                    });
}

int krome_get_photo_heating(int ncell, const double* n, const double* Tgas,
                            const double* par, int par_stride,
                            const double* Phi, const krome_opts* h,
                            double* heat) {
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (Tgas == nullptr || heat == nullptr || h == nullptr || !h->has_photo)
    return -1;
  return for_cells_(ncell, n, Tgas, par, par_stride, Phi, h,
                    [&](int c, const double* nc, double T,
                        const krome::Params& P, const krome::SolverOpts& o) {
                      return krome::get_photo_heating(
                          nc, T, P, o,
                          heat + static_cast<size_t>(c) * krome::NPHOTO);
                    });
#else
  (void)ncell; (void)n; (void)Tgas; (void)par; (void)par_stride; (void)Phi;
  (void)h; (void)heat;
  return -1;
#endif
}

int krome_get_ratio_flux_draine(int ncell, const double* Phi,
                                const krome_opts* h, double* ratio) {
  if (ncell < 0 || Phi == nullptr || ratio == nullptr || h == nullptr ||
      !h->has_photo)
    return -1;
#ifdef KROME_PHOTO_HPP_INCLUDED
  krome::SolverOpts o = h->so;
  o.photo = &h->photo;
  for (int c = 0; c < ncell; ++c) {
    o.Phi = const_cast<double*>(Phi + static_cast<size_t>(c) * h->photo.B);
    const int r = krome::get_ratio_flux_draine(o, ratio + c);
    if (r != 0) return r;
  }
  return 0;
#else
  return -1;
#endif
}

int krome_get_element_totals(int ncell, const double* n, double* tot) {
  if (ncell < 0 || n == nullptr || tot == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    krome::get_element_totals(n + static_cast<size_t>(c) * krome::NSPEC,
                              tot + static_cast<size_t>(c) * krome::NELEM);
  return 0;
}

int krome_get_element_totals_x(int ncell, const double* x, double* tot) {
  if (ncell < 0 || x == nullptr || tot == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    krome::get_element_totals_x(x + static_cast<size_t>(c) * krome::NSPEC,
                                tot + static_cast<size_t>(c) * krome::NELEM);
  return 0;
}

// consistent advection: the parameter rows as Params (one per cell, or one)
static int conserve_cells_c_(int ncell, double* v, const double* ref,
                             int ref_stride, const double* par,
                             int par_stride, int* flags, bool x) {
  if (ncell < 0 || v == nullptr || ref == nullptr || ref_stride < 0
      || par_stride < 0)
    return -1;
  if (ref_stride != 0 && ref_stride != krome::NELEM) return -1;
  const bool per_cell = par != nullptr && par_stride != 0;
  std::vector<krome::Params> P(per_cell ? ncell : 1);
  if (par != nullptr)
    for (size_t c = 0; c < P.size(); ++c)
      krome::params_from_row(par + c * par_stride, P[c]);
  const int rs = ref_stride != 0 ? 1 : 0;
  return x ? krome::conserve_cells_x(ncell, v, ref, rs, P.data(),
                                     per_cell ? 1 : 0, flags)
           : krome::conserve_cells(ncell, v, ref, rs, P.data(),
                                   per_cell ? 1 : 0, flags);
}

int krome_conserve_cells(int ncell, double* n, const double* ref,
                         int ref_stride, const double* par, int par_stride,
                         int* flags) {
  return conserve_cells_c_(ncell, n, ref, ref_stride, par, par_stride, flags,
                           false);
}

int krome_conserve_cells_x(int ncell, double* x, const double* ref,
                           int ref_stride, const double* par,
                           int par_stride, int* flags) {
  return conserve_cells_c_(ncell, x, ref, ref_stride, par, par_stride, flags,
                           true);
}

int krome_n2x(int ncell, const double* n, double* x) {
  if (ncell < 0 || n == nullptr || x == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    krome::n2x(n + static_cast<size_t>(c) * krome::NSPEC,
               x + static_cast<size_t>(c) * krome::NSPEC);
  return 0;
}

int krome_x2n(int ncell, const double* x, const double* rho, double* n) {
  if (ncell < 0 || x == nullptr || rho == nullptr || n == nullptr) return -1;
  for (int c = 0; c < ncell; ++c)
    krome::x2n(x + static_cast<size_t>(c) * krome::NSPEC, rho[c],
               n + static_cast<size_t>(c) * krome::NSPEC);
  return 0;
}

/* ---- radiation --------------------------------------------------------- */

namespace {
bool bins_ok_(int B, const double* Elo, const double* Ehi) {
  if (B <= 0 || Elo == nullptr || Ehi == nullptr) return false;
  for (int b = 0; b < B; ++b)
    if (!(Ehi[b] > Elo[b]) || !(Elo[b] > 0.0)) return false;
  return true;
}
}  // namespace

int krome_bins_lin(int nbins, double lower, double upper, double* Elo,
                   double* Ehi) {
  if (nbins <= 0 || !(upper > lower) || !(lower > 0.0) || !Elo || !Ehi)
    return -1;
  krome::spectra::bins_lin(nbins, lower, upper, Elo, Ehi);
  return 0;
}
int krome_bins_log(int nbins, double lower, double upper, double* Elo,
                   double* Ehi) {
  if (nbins <= 0 || !(upper > lower) || !(lower > 0.0) || !Elo || !Ehi)
    return -1;
  krome::spectra::bins_log(nbins, lower, upper, Elo, Ehi);
  return 0;
}
int krome_J_to_phi(int nbins, const double* Elo, const double* Ehi,
                   const double* J, double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || !J || !phi) return -1;
  for (int b = 0; b < nbins; ++b)
    phi[b] = krome::spectra::J_to_phi(J[b], Elo[b], Ehi[b]);
  return 0;
}
int krome_phi_to_J(int nbins, const double* Elo, const double* Ehi,
                   const double* phi, double* J) {
  if (!bins_ok_(nbins, Elo, Ehi) || !J || !phi) return -1;
  for (int b = 0; b < nbins; ++b)
    J[b] = krome::spectra::phi_to_J(phi[b], Elo[b], Ehi[b]);
  return 0;
}
int krome_spectrum_blackbody(int nbins, const double* Elo, const double* Ehi,
                             double Tbb, double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || !(Tbb > 0.0) || !phi) return -1;
  krome::spectra::blackbody(nbins, Elo, Ehi, Tbb, phi);
  return 0;
}
int krome_blackbody_range(double Tbb, double* lower, double* upper) {
  if (!(Tbb > 0.0) || !lower || !upper) return -1;
  return krome::spectra::blackbody_range(Tbb, lower, upper);
}
int krome_spectrum_draine(int nbins, const double* Elo, const double* Ehi,
                          double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || !phi) return -1;
  krome::spectra::draine(nbins, Elo, Ehi, phi);
  return 0;
}
int krome_spectrum_j21(int nbins, const double* Elo, const double* Ehi,
                       double J21, double alpha, double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || !phi) return -1;
  krome::spectra::j21(nbins, Elo, Ehi, J21, alpha, phi);
  return 0;
}
int krome_spectrum_hm2012(int nbins, const double* Elo, const double* Ehi,
                          double z, double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || !phi) return -1;
  return krome::spectra::hm2012(nbins, Elo, Ehi, z, phi);
}
int krome_spectrum_table(int nbins, const double* Elo, const double* Ehi,
                         int nE, const double* E, const double* J,
                         int loglog, double* phi) {
  if (!bins_ok_(nbins, Elo, Ehi) || nE < 2 || !E || !J || !phi) return -1;
  for (int i = 1; i < nE; ++i)
    if (!(E[i] > E[i - 1])) return -1;
  krome::spectra::table(nbins, Elo, Ehi, nE, E, J, loglog != 0, phi);
  return 0;
}
int krome_dust_kappa_bins(int nE, const double* E, const double* kappa,
                          int nbins, const double* Elo, const double* Ehi,
                          double* kbin) {
  if (!bins_ok_(nbins, Elo, Ehi) || nE < 2 || !E || !kappa || !kbin)
    return -1;
  for (int i = 1; i < nE; ++i)
    if (!(E[i] > E[i - 1])) return -1;
  return krome::spectra::dust_kappa_bins(nE, E, kappa, nbins, Elo, Ehi,
                                         kbin);
}
int krome_dust_kappa_draine(int nbins, const double* Elo, const double* Ehi,
                            double* kbin) {
  if (!bins_ok_(nbins, Elo, Ehi) || !kbin) return -1;
  return krome::spectra::dust_kappa_draine(nbins, Elo, Ehi, kbin);
}
int krome_find_G0_Av(int nbins, const double* Elo, const double* Ehi,
                     const double* phi, const double* kbin, double d2g,
                     double xH, double* G0, double* Av) {
  if (!bins_ok_(nbins, Elo, Ehi) || !phi || !kbin || !G0 || !Av ||
      !(d2g > 0.0) || !(xH > 0.0))
    return -1;
  return krome::spectra::find_G0_Av(nbins, Elo, Ehi, phi, kbin, d2g, xH, G0,
                                    Av);
}

int krome_nphoto(void) { return krome::NPHOTO; }
int krome_photo_reaction(int r) {
  return (r >= 0 && r < krome::NPHOTO) ? photo_reaction[r] : -1;
}
int krome_photo_reactant(int r) {
  return (r >= 0 && r < krome::NPHOTO) ? photo_reactant[r] : -1;
}
double krome_photo_pdiss(int r) {
  return (r >= 0 && r < krome::NPHOTO) ? photo_pdiss[r]
                                       : std::nan("");
}
int krome_nlw(void) { return krome::NLW; }
int krome_lw_reaction(int q) {
  return (q >= 0 && q < krome::NLW) ? lw_reaction[q] : -1;
}
int krome_lw_reactant(int q) {
  return (q >= 0 && q < krome::NLW) ? lw_reactant[q] : -1;
}
double krome_lw_pdiss(int q) {
  return (q >= 0 && q < krome::NLW) ? lw_pdiss[q] : std::nan("");
}
int krome_opts_lw_weights(const krome_opts* h, double* w) {
  if (h == nullptr || !h->has_photo || w == nullptr) return -1;
  for (int b = 0; b < h->photo.B; ++b) w[b] = h->photo.lww[b];
  return 0;
}
double krome_photo_xsec(int r, double E) {
  if (r < 0 || r >= krome::NPHOTO) return std::nan("");
  return (E >= photo_Emin[r] && E <= photo_Emax[r])
      ? krome::photo_xsec[r](E) : 0.0;
}
int krome_opts_photo_sigma(const krome_opts* h, double* sigma) {
  if (h == nullptr || !h->has_photo || sigma == nullptr) return -1;
  for (size_t i = 0; i < h->photo.sigma.size(); ++i)
    sigma[i] = h->photo.sigma[i];
  return 0;
}

int krome_get_habing(int ncell, const double* Phi, const krome_opts* h,
                     double* G) {
  if (ncell < 0 || h == nullptr || !h->has_photo || !Phi || !G) return -1;
  for (int c = 0; c < ncell; ++c)
    G[c] = h->photo.habing(Phi + static_cast<size_t>(c) * h->photo.B);
  return 0;
}

int krome_get_opacity(int ncell, const double* n, const double* Tgas,
                      const double* par, int par_stride,
                      const krome_opts* h, const double* d2g,
                      const double* kdust, double* tau) {
  if (Tgas == nullptr || tau == nullptr || h == nullptr || !h->has_photo)
    return -1;
  const double* kd = kdust != nullptr ? kdust
      : (h->kdraine.empty() ? nullptr : h->kdraine.data());
  if (d2g != nullptr && !h->has_dust && kd == nullptr) return -1;
  const size_t B = static_cast<size_t>(h->photo.B);
  return for_cells_(
      ncell, n, Tgas, par, par_stride, nullptr, h,
      [&](int c, const double* nc, double T, const krome::Params& P,
          const krome::SolverOpts& o) {
        return krome::get_opacity(nc, T, P, o,
                                  d2g != nullptr ? d2g[c] : 0.0, kd,
                                  tau + c * B);
      },
      false);
}

int krome_get_absorption(int ncell, const double* n, const double* Tgas,
                         const double* par, int par_stride,
                         const krome_opts* h, const double* kdust,
                         double* kabs) {
  if (Tgas == nullptr || kabs == nullptr || h == nullptr || !h->has_photo)
    return -1;
  const size_t B = static_cast<size_t>(h->photo.B);
  return for_cells_(
      ncell, n, Tgas, par, par_stride, nullptr, h,
      [&](int c, const double* nc, double T, const krome::Params& P,
          const krome::SolverOpts& o) {
        return krome::get_absorption(nc, T, P, o, kdust, kabs + c * B);
      },
      false);
}

/* ---- species / reaction metadata ---------------------------------------- */
// (the numeric tables are unqualified: they are macros for the generated
// device-safe accessors, ::krome::<name>_arr_())

double krome_species_mass(int i) {
  return (i >= 0 && i < krome::NSPEC) ? species_mass_g[i]
                                      : std::nan("");
}
int krome_species_charge(int i) {
  return (i >= 0 && i < krome::NSPEC) ? species_charge[i] : 0;
}
int krome_nelem(void) { return krome::NELEM; }
const char* krome_element_name(int e) {
  return (e >= 0 && e < krome::NELEM) ? krome::element_names[e] : nullptr;
}
int krome_species_natoms(int i, int e) {
  if (i < 0 || i >= krome::NSPEC || e < 0 || e >= krome::NELEM) return -1;
  return species_elements[i * krome::NELEM + e];
}
const char* krome_reaction_name(int r) {
  return (r >= 0 && r < krome::NREA) ? krome::reaction_names[r] : nullptr;
}

/* ---- device backend ------------------------------------------------------ */

int krome_has_device(void) {
#ifdef KROME_USE_KOKKOS
  return 1;
#else
  return 0;
#endif
}
int krome_device_initialize(void) {
#ifdef KROME_USE_KOKKOS
  if (!Kokkos::is_initialized()) Kokkos::initialize();
  return 0;
#else
  return krome::FLAG_UNSUPPORTED;
#endif
}
int krome_device_initialize_id(int device_id) {
#ifdef KROME_USE_KOKKOS
  if (Kokkos::is_initialized()) return 1;
  Kokkos::InitializationSettings st;
  if (device_id >= 0) st.set_device_id(device_id);
  Kokkos::initialize(st);
  return 0;
#else
  (void)device_id;
  return krome::FLAG_UNSUPPORTED;
#endif
}

#ifdef KROME_USE_KOKKOS
static_assert(KROME_LAYOUT_SPECIES_MAJOR == krome::kk::LAYOUT_SPECIES_MAJOR
                  && KROME_LAYOUT_CELL_MAJOR == krome::kk::LAYOUT_CELL_MAJOR,
              "layouts");
#endif

static int solve_cells_device_c_(int ncell, double* v, const double* rho,
                                 int layout, long ld, double* Tgas,
                                 const double* dt, const double* par,
                                 int par_stride, const krome_cells* extra,
                                 int* flags, krome_device_summary* summary,
                                 const krome_opts* hc) {
  if (summary != nullptr) *summary = krome_device_summary{0, 0, 0};
  if (hc == nullptr || ncell < 0 || v == nullptr || Tgas == nullptr
      || dt == nullptr || (layout != KROME_LAYOUT_SPECIES_MAJOR
                           && layout != KROME_LAYOUT_CELL_MAJOR)
      || (ld != 0 && ld < ncell))
    return -1;
  if (par != nullptr && par_stride != 0 && par_stride != krome::NPAR)
    return -1;
  if (ncell == 0) return 0;
#ifdef KROME_USE_KOKKOS
  krome::SolverOpts so = hc->so;
  so.error_file = hc->error_file.c_str();
  krome::kk::DeviceCells a;
  a.n = v;
  a.rho = rho;
  a.layout = layout;
  a.ld = ld;
  a.Tgas = Tgas;
  a.dt = dt;
  a.par = par;
  a.par_stride = par != nullptr ? par_stride : 0;
  a.flags = flags;
  if (hc->has_photo) {
    so.photo = &hc->photo;
    a.Phi = extra != nullptr ? extra->Phi : nullptr;
  }
  if (hc->has_dust) {
    so.dust = const_cast<krome::DustBins*>(&hc->dust);
    a.Td = extra != nullptr ? extra->Td : nullptr;
    a.adust = extra != nullptr ? extra->adust : nullptr;
  }
  krome::kk::DeviceSummary sm;
  const int r = krome::kk::solve_cells_device(ncell, a, so, hc->co, &sm);
  if (summary != nullptr)
    *summary = krome_device_summary{sm.device, sm.host, sm.failed};
  return r;
#else
  (void)rho; (void)extra; (void)flags;
  return krome::FLAG_UNSUPPORTED;
#endif
}

int krome_solve_cells_device(int ncell, double* n, int layout, long ld,
                             double* Tgas, const double* dt,
                             const double* par, int par_stride,
                             const krome_cells* extra, int* flags,
                             krome_device_summary* summary,
                             const krome_opts* h) {
  return solve_cells_device_c_(ncell, n, nullptr, layout, ld, Tgas, dt, par,
                               par_stride, extra, flags, summary, h);
}

int krome_solve_cells_x_device(int ncell, double* x, const double* rho,
                               int layout, long ld, double* Tgas,
                               const double* dt, const double* par,
                               int par_stride, const krome_cells* extra,
                               int* flags, krome_device_summary* summary,
                               const krome_opts* h) {
  if (rho == nullptr) return -1;
  return solve_cells_device_c_(ncell, x, rho, layout, ld, Tgas, dt, par,
                               par_stride, extra, flags, summary, h);
}

static int conserve_cells_device_c_(int ncell, double* v, int layout,
                                    long ld, const double* ref,
                                    int ref_stride, const double* par,
                                    int par_stride, int* flags, bool x) {
  if (ncell < 0 || v == nullptr || ref == nullptr
      || (ref_stride != 0 && ref_stride != krome::NELEM)
      || (par != nullptr && par_stride != 0 && par_stride != krome::NPAR)
      || (layout != KROME_LAYOUT_SPECIES_MAJOR
          && layout != KROME_LAYOUT_CELL_MAJOR)
      || (ld != 0 && ld < ncell))
    return -1;
#ifdef KROME_USE_KOKKOS
  return krome::kk::conserve_cells_device(ncell, v, x, layout, ld, ref,
                                          ref_stride, par,
                                          par != nullptr ? par_stride : 0,
                                          flags);
#else
  (void)flags; (void)x;
  return krome::FLAG_UNSUPPORTED;
#endif
}

int krome_conserve_cells_device(int ncell, double* n, int layout, long ld,
                                const double* ref, int ref_stride,
                                const double* par, int par_stride,
                                int* flags) {
  return conserve_cells_device_c_(ncell, n, layout, ld, ref, ref_stride, par,
                                  par_stride, flags, false);
}

int krome_conserve_cells_x_device(int ncell, double* x, int layout, long ld,
                                  const double* ref, int ref_stride,
                                  const double* par, int par_stride,
                                  int* flags) {
  return conserve_cells_device_c_(ncell, x, layout, ld, ref, ref_stride, par,
                                  par_stride, flags, true);
}

int krome_device_finalize(void) {
#ifdef KROME_USE_KOKKOS
  if (Kokkos::is_initialized() && !Kokkos::is_finalized()) Kokkos::finalize();
  return 0;
#else
  return krome::FLAG_UNSUPPORTED;
#endif
}

}  // extern "C"
