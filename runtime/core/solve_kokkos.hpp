// KROME-X Kokkos backend: the cells of a solve_cells call integrated on
// a Kokkos execution space (CUDA, HIP, SYCL, or OpenMP / Serial on the
// host),
//
//   solve_cells_kokkos(ncell, n[ncell][NSPEC], Tgas[ncell], dt[ncell],
//                      params, pstride, opts, cells_opts)
//
// All cells are copied to the device in one call and integrated by one
// kernel, one thread per cell, each running its own integration to its own
// dt (core/cvode/bdf_cell.hpp: CVODE's BDF, ported, with its own step size,
// order, Newton iteration, Jacobian reuse and error test) on the host
// solver's own rhs and Jacobian (solve_cell.hpp: full_rhs_, jac_cell_ over
// a RhsCtx), so that the device integrates exactly the host's equations:
// the same state layout [n | Tgas | ln Phi | Tdust | grain atoms], the
// same JAC_AUTO (difference quotients with Tgas evolved, the analytic
// Jacobian with its photo-shielding, radiation, dust and
// secondary-ionization blocks at fixed Tgas), per-species rtol. Every
// configuration of the host is supported: fixed or evolved Tgas, static
// or co-evolved photon bins (radiation.evolve), binned dust (DUSTBINS
// cooling, the surface hook, surface.bins), dust.evolve_Td, dust.evolve,
// radiation.secondary_ionization, the dust tables. Each cell returns its own
// flag; the cells that fail are re-solved on the host (solve()) from
// their inputs.
//
// Per-cell inputs and outputs as on the host: CellsOpts::Phi rows (read;
// with radiation.evolve updated), CellsOpts::Td rows (the bins' Td: read and
// updated), CellsOpts::adust rows (dust.evolve radii: read and updated);
// without a row, every cell starts from the SolverOpts' field / DustBins
// (with radiation.evolve a single-cell call updates opts.Phi).
//
// Include the network, the runtime modules, solve_cell.hpp and
// solve_cells.hpp before this header; compile with Kokkos (and the
// matching device compiler). Kokkos must be initialized by the caller.
#pragma once
#include <algorithm>
#include <climits>
#include <cmath>
#include <exception>
#include <mutex>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <Kokkos_Core.hpp>

namespace krome {
namespace kk {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
using DView = Kokkos::View<double*, MemSpace>;
using DView2 = Kokkos::View<double**, Kokkos::LayoutRight, MemSpace>;
using IView1 = Kokkos::View<int*, MemSpace>;
using LView2 = Kokkos::View<long**, Kokkos::LayoutRight, MemSpace>;
using TView = Kokkos::View<ThermCtx*, MemSpace>;

// the static field's per-reaction photon sums of a cell (photo_sums): five
// arrays of NPHOTO, back to back in one row of Batch::psum
constexpr int NPSUM_ = 5;
#ifdef KROME_PHOTO_HPP_INCLUDED
// point v's sums (s, se, ss, sc, sp) into the row s[NPSUM_ * NPHOTO]
KOKKOS_INLINE_FUNCTION void bind_psum_(double* s, PhotoSumsView& v) {
  v.s = s;
  v.se = s + NPHOTO;
  v.ss = s + 2 * NPHOTO;
  v.sc = s + 3 * NPHOTO;
  v.sp = s + 4 * NPHOTO;
}
#endif

// what one batch needs on the device (Views: reference-counted handles,
// captured by value in the kernels)
struct Batch {
  int ncell = 0;
  int blk = NSPEC;   // equations per cell
  // the state layout (the host solve()'s, init_ctx_)
  int iT = -1;             // Tgas slot (-1: fixed)
  int iu = -1, Bev = 0;    // ln Phi slots (radiation.evolve)
  int id = -1, D = 0;      // Tdust slots (dust.evolve_Td, Tgas evolved)
  int ia = -1, A = 0;      // grain-atom slots (dust.evolve)
  int nb = 0;              // dust bins (0: no DustBins)
  bool balance = false;    // the bins' Td from the balance (no dust.evolve_Td)
  long rw = 0;             // per-cell rhs work doubles (+ 2 blk scratch)
  TView tc;      // per cell: params, k, shields, photo and dust views
  DView dt, Tfix;                        // per cell
  DView2 n0;         // [ncell][NSPEC] input densities (n_global)
  DView2 phi;        // [ncell][B] per-cell photon field (radiation.evolve: out)
  DView2 psum;       // [ncell][NPSUM_ * NPHOTO] static-field bin sums
#ifdef KROME_PHOTO_HPP_INCLUDED
  Kokkos::View<PhotoSumsView*, MemSpace> psv;  // per cell: views of psum
#endif
  DView sigma, sigmaE, sigmaEs, sigmaCs, sigmaPs, habw, elo, ehi, lww;
  // the dust population: per-cell Td and radii rows (in, and out), and its
  // tables (one copy for every cell)
  DView2 tdr, adr;   // [ncell][nb]
  DView dxw, dibb, dibbs, dqmid, demid;
  IView1 disC;
  int dnT = 0, dB = 0;
  double dTmin = 1.0, dTmax = 1e4, dTmult = 0.0;
};

// ---- device buffers -------------------------------------------------------
//
// The device arrays of a call (the Batch views, the per-cell state,
// workspace, flags and statistics) and their host staging copies. With
// CellsOpts::device_cache (default) one set is kept between calls and
// reused: allocated for the largest call so far (grown, never shrunk),
// released by release_device_cache() or at Kokkos::finalize. Only memory
// is kept: every call writes its own inputs and re-initializes every cell
// (a call never reads what a previous one left; the dust tables are
// uploaded again when the population changes). A call while the cache is
// in use (another host thread), or with device_cache off, gets buffers of
// its own for that call.

// the integrator's statistics per cell (the columns of DeviceBuffers_::
// stats, copied into SolveStats)
enum StatCol_ {
  ST_STEPS_, ST_RHS_, ST_JAC_, ST_RHS_JAC_, ST_LSETUPS_, ST_NEWTON_,
  ST_NEWTON_FAILS_, ST_ERR_FAILS_, ST_RESTARTS_,
  NSTAT_  // the number of columns
};

struct DeviceBuffers_ {
  int cap = 0;    // cells
  int blk = -1;   // equations per cell
  int B = -1;     // photon bins
  int nb = -1;    // dust bins
  int pcap = 0;   // cells of the per-cell arrays
  long rw = -1;   // per-cell work doubles of the per-cell arrays
  long wd = -1;   // per-cell integrator doubles (bdf::work_doubles: dense
                  // or sparse Newton matrices)
  Batch b;
  TView::HostMirror tch;
  DView::HostMirror dth, Tfh;
  DView2::HostMirror n0h, phh, tdh, adh;
  DView y, work;  // [cap * blk], [cap * wd]: interleaved or
  IView1 iwork;   // cell-contiguous (per-cell); iwork [cap * blk]
  DView scr;      // [cap * rw]: the contiguous layout's rhs work rows
  IView1 flags;
  LView2 stats;
  DView::HostMirror yh;
  IView1::HostMirror fh;
  LView2::HostMirror sh;
  DView atol, rtolv;
  // the dust intBB tables uploaded (their host address and size)
  const void* dtab_key = nullptr;
  long dtab_n = -1;
};

// the cache and its use counters (diagnostics: allocations since start)
inline DeviceBuffers_* device_cache_ = nullptr;
inline std::mutex device_cache_mutex_;
inline bool device_cache_hook_ = false;
inline long device_allocs_ = 0;

// free the cached device buffers (the next call allocates again); called
// by Kokkos::finalize (a finalize hook), callable any time between calls
inline void release_device_cache() {
  std::lock_guard<std::mutex> lk(device_cache_mutex_);
  if (device_cache_ != nullptr) *device_cache_ = DeviceBuffers_();
}

// the bytes a set of buffers holds on the device
inline size_t buffers_bytes_(const DeviceBuffers_& D) {
  const Batch& b = D.b;
  return b.tc.span() * sizeof(ThermCtx)
         + (b.dt.span() + b.Tfix.span() + b.n0.span() + b.phi.span()
            + b.psum.span() + b.tdr.span() + b.adr.span() + b.dxw.span()
            + b.dibb.span() + b.dibbs.span() + b.dqmid.span()
            + b.demid.span() + D.y.span() + D.scr.span() + D.work.span()
            + D.atol.span() + D.rtolv.span()) * sizeof(double)
         + (D.iwork.span() + D.flags.span() + b.disC.span()) * sizeof(int)
         + D.stats.span() * sizeof(long);
}

// cached bytes on the device (the per-cell arrays included)
inline size_t device_cache_bytes() {
  std::lock_guard<std::mutex> lk(device_cache_mutex_);
  return device_cache_ == nullptr ? 0 : buffers_bytes_(*device_cache_);
}

// tests: the kernel reports every cell c with c % device_fail_every_ == 0
// (a chunk's index) as failed (FLAG_TOO_MUCH_WORK), to exercise the host
// fallback; 0: off
inline int device_fail_every_ = 0;
// tests: an allocation for more than this many cells throws std::bad_alloc
// (a device out of memory, simulated); < 0: off
inline long alloc_fail_cells_ = -1;
inline void alloc_check_(long cap) {
  if (alloc_fail_cells_ >= 0 && cap > alloc_fail_cells_)
    throw std::bad_alloc();
}

// the first n rows of a view
template <class V>
inline auto rows1_(const V& v, int n) {
  return Kokkos::subview(v, std::make_pair(0, n));
}
template <class V>
inline auto rows2_(const V& v, int n) {
  return Kokkos::subview(v, std::make_pair(0, n), Kokkos::ALL);
}

// a host table on the device (its own size: re-allocated when it changes,
// rewritten on every call)
inline void upload_(DView& d, const std::vector<double>& h, const char* name) {
  const size_t m = h.size() > 0 ? h.size() : 1;
  if (d.extent(0) != m) {
    d = DView(std::string(name), m);
    ++device_allocs_;
  }
  if (h.empty()) return;
  Kokkos::View<const double*, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>> hv(h.data(), h.size());
  Kokkos::deep_copy(d, hv);
}

// an int table on the device (as upload_)
inline void upload_i_(IView1& d, const std::vector<int>& h, const char* name) {
  const size_t m = h.size() > 0 ? h.size() : 1;
  if (d.extent(0) != m) {
    d = IView1(std::string(name), m);
    ++device_allocs_;
  }
  if (h.empty()) return;
  Kokkos::View<const int*, Kokkos::HostSpace,
               Kokkos::MemoryTraits<Kokkos::Unmanaged>> hv(h.data(), h.size());
  Kokkos::deep_copy(d, hv);
}

// the Batch arrays for ncell cells of blk equations, B photon bins and nb
// dust bins: kept when they fit
inline void ensure_batch_(DeviceBuffers_& D, int ncell, int blk, int B,
                          int nb) {
  const bool same = blk == D.blk && B == D.B && nb == D.nb;
  if (ncell <= D.cap && same) return;
  const int cap = std::max(ncell, same ? D.cap : 0);
  // the old buffers freed first: the peak is the new ones alone
  D = DeviceBuffers_();
  alloc_check_(cap);
  Batch& b = D.b;
  b.tc = TView("krome_tc", cap);
  b.dt = DView("krome_dt", cap);
  b.Tfix = DView("krome_Tfix", cap);
  b.n0 = DView2("krome_n0", cap, NSPEC);
  b.phi = DView2();
  b.psum = DView2();
#ifdef KROME_PHOTO_HPP_INCLUDED
  b.psv = Kokkos::View<PhotoSumsView*, MemSpace>();
  if (B > 0) {
    b.phi = DView2("krome_phi", cap, B);
    b.psum = DView2("krome_psum", cap, NPSUM_ * (NPHOTO > 0 ? NPHOTO : 1));
    b.psv = Kokkos::View<PhotoSumsView*, MemSpace>("krome_psv", cap);
  }
#endif
  b.tdr = DView2();
  b.adr = DView2();
  if (nb > 0) {
    b.tdr = DView2("krome_tdr", cap, nb);
    b.adr = DView2("krome_adr", cap, nb);
  }
  D.tch = Kokkos::create_mirror_view(b.tc);
  D.dth = Kokkos::create_mirror_view(b.dt);
  D.Tfh = Kokkos::create_mirror_view(b.Tfix);
  D.n0h = Kokkos::create_mirror_view(b.n0);
  D.phh = Kokkos::create_mirror_view(b.phi);
  D.tdh = Kokkos::create_mirror_view(b.tdr);
  D.adh = Kokkos::create_mirror_view(b.adr);
  D.cap = cap;
  D.blk = blk;
  D.B = B;
  D.nb = nb;
  D.pcap = 0;  // the per-cell arrays follow the new capacity
  ++device_allocs_;
}

// the per-cell arrays at the Batch capacity, wd integrator doubles per cell
inline void ensure_percell_(DeviceBuffers_& D, long wd) {
  if (D.pcap == D.cap && D.cap > 0 && D.rw == D.b.rw && D.wd == wd) return;
  const int cap = D.cap, blk = D.blk;
  // the old arrays freed first
  D.y = DView();
  D.work = DView();
  D.iwork = IView1();
  D.scr = DView();
  D.flags = IView1();
  D.stats = LView2();
  D.yh = DView::HostMirror();
  D.fh = IView1::HostMirror();
  D.sh = LView2::HostMirror();
  D.pcap = 0;
  alloc_check_(cap);
  D.y = DView("krome_y", static_cast<size_t>(cap) * blk);
  D.work = DView("krome_bdf_work", static_cast<size_t>(cap) * wd);
  D.wd = wd;
  D.iwork = IView1("krome_bdf_iwork",
                   static_cast<size_t>(cap) * bdf::work_ints(blk));
  D.scr = DView("krome_rhs_work", static_cast<size_t>(cap) * D.b.rw);
  D.flags = IView1("krome_flags", cap);
  D.stats = LView2("krome_stats", cap, NSTAT_);
  D.yh = Kokkos::create_mirror_view(D.y);
  D.fh = Kokkos::create_mirror_view(D.flags);
  D.sh = Kokkos::create_mirror_view(D.stats);
  D.pcap = cap;
  D.rw = D.b.rw;
  ++device_allocs_;
}

// ---- the backend ---------------------------------------------------------

// the state layout of a call (the host solve()'s, init_ctx_) and its
// photon (B) and dust (nb) bins: 0, or the input flags of setup_batch_
struct Shape_ {
  int blk = NSPEC, iT = -1, iu = -1, Bev = 0, id = -1, Dd = 0, ia = -1,
      A = 0, nb = 0, B = 0;
  bool balance = false;
};
// a flag batch_shape_ returns: the options and inputs are wrong for every
// cell of the call alike (nothing is integrated)
inline bool input_flag_(int f) {
  return f == FLAG_NO_DUST || f == FLAG_NO_FIELD
         || f == FLAG_CONSERVE_GAS_TOTAL || f == FLAG_SURFACE_BINS;
}

inline int batch_shape_(const SolverOpts& o, const CellsOpts& co,
                        Shape_& sh) {
  int blk = NSPEC, iT = -1, iu = -1, Bev = 0, id = -1, Dd = 0, ia = -1,
      A = 0, nb = 0, B = 0;
  bool balance = false;
  if (o.evolve_T) iT = blk++;
#ifdef KROME_PHOTO_HPP_INCLUDED
  const PhotoModule* ph = o.photo;
  if (ph != nullptr) {
    if (o.Phi == nullptr && co.Phi == nullptr) return FLAG_NO_FIELD;
    B = ph->B;
    if constexpr (EvolveConfig::rad) {
      Bev = B;
      iu = blk;
      blk += B;
    }
  } else if (EvolveConfig::rad) {
    return FLAG_NO_FIELD;
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if constexpr (SURF_NB > 0)
    if (o.dust == nullptr || o.dust->nb != SURF_NB) return FLAG_SURFACE_BINS;
  if (TC::cool_dust_bins && o.evolve_T && o.dust == nullptr)
    return FLAG_NO_DUST;
  nb = o.dust != nullptr ? o.dust->nb : 0;
  if (EvolveConfig::Td && o.evolve_T) {
    Dd = nb;
    id = blk;
    blk += Dd;
  }
  if constexpr (EvolveConfig::dust) {
    if (o.dust == nullptr) return FLAG_NO_DUST;
    if (o.conserve == 2) return FLAG_CONSERVE_GAS_TOTAL;
    A = nb;
    ia = blk;
    blk += A;
  }
#ifdef KROME_HAS_PAR_d2g
  balance = !EvolveConfig::Td && o.dust != nullptr;
#endif
#endif
  (void)co;
  sh.blk = blk;
  sh.iT = iT;
  sh.iu = iu;
  sh.Bev = Bev;
  sh.id = id;
  sh.Dd = Dd;
  sh.ia = ia;
  sh.A = A;
  sh.nb = nb;
  sh.B = B;
  sh.balance = balance;
  return 0;
}

// the device's sparse LU: LINSOL_SPARSE, or AUTO where the host solve
// takes it (sparse_auto_; the device has no LAPACK)
inline bool percell_sparse_(const SolverOpts& o, int blk) {
  return o.linsol == LINSOL_SPARSE
         || (o.linsol == LINSOL_AUTO && sparse_auto_(blk));
}

// the device bytes per cell of a call (ensure_batch_ and ensure_percell_);
// 0 when the inputs are rejected (setup_batch_ reports them)
inline long long percell_bytes_(const SolverOpts& o, const CellsOpts& co) {
  Shape_ sh;
  if (batch_shape_(o, co, sh) != 0) return 0;
  const int blk = sh.blk;
  const long long rw = rhs_work_doubles(blk, sh.Bev, sh.nb, sh.Dd, sh.A)
                       + 2LL * blk;
  long long dbl = 2 + NSPEC + blk + rw   // dt, Tfix, n0, y, rhs work
      + bdf::work_doubles(blk, percell_sparse_(o, blk), NSPEC, SPARSE_NNZ)
      + 2LL * sh.nb;                     // Td and radius rows
  long long bytes = static_cast<long long>(sizeof(ThermCtx));
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (sh.B > 0) {
    dbl += sh.B + 5LL * (NPHOTO > 0 ? NPHOTO : 1);
    bytes += static_cast<long long>(sizeof(PhotoSumsView));
  }
#endif
  return bytes + dbl * static_cast<long long>(sizeof(double))
         + (bdf::work_ints(blk) + 1LL) * static_cast<long long>(sizeof(int))
         + NSTAT_ * static_cast<long long>(sizeof(long));
}

// the free memory of the device [bytes] (CUDA, HIP); -1: unknown (the
// host spaces and the others: no limit applied)
inline long long device_free_bytes_() {
#if defined(KOKKOS_ENABLE_CUDA)
  size_t f = 0, t = 0;
  if (cudaMemGetInfo(&f, &t) == cudaSuccess)
    return static_cast<long long>(f);
#elif defined(KOKKOS_ENABLE_HIP)
  size_t f = 0, t = 0;
  if (hipMemGetInfo(&f, &t) == hipSuccess) return static_cast<long long>(f);
#endif
  return -1;
}

// an exception from a device allocation (or a launch) out of memory
inline bool out_of_memory_(const std::exception& e) {
  if (dynamic_cast<const std::bad_alloc*>(&e) != nullptr) return true;
  const std::string w = e.what();
  for (const char* k : {"out of memory", "MemoryAllocation",
                        "failed to allocate", "insufficient memory"})
    if (w.find(k) != std::string::npos) return true;
  return false;
}

// the Batch arrays and state layout for ncell cells of the shape sh
inline void set_shape_(DeviceBuffers_& D, const Shape_& sh, int ncell) {
  const int blk = sh.blk, Bev = sh.Bev, Dd = sh.Dd, A = sh.A, nb = sh.nb;
  ensure_batch_(D, ncell, blk, sh.B, nb);
  Batch& b = D.b;
  b.ncell = ncell;
  b.blk = blk;
  b.iT = sh.iT;
  b.iu = sh.iu;
  b.Bev = Bev;
  b.id = sh.id;
  b.D = Dd;
  b.ia = sh.ia;
  b.A = A;
  b.nb = nb;
  b.balance = sh.balance;
  b.rw = rhs_work_doubles(blk, Bev, nb, Dd, A) + 2L * blk;
}

// the tables of the photo module and of the dust population on the device
// (re-uploaded when they change); B photon bins, nb dust bins
inline void upload_tables_(DeviceBuffers_& D, const SolverOpts& o, int B,
                           int nb) {
  Batch& b = D.b;
  (void)b; (void)o; (void)B; (void)nb;
#ifdef KROME_PHOTO_HPP_INCLUDED
  const PhotoModule* ph = o.photo;
  if (ph != nullptr) {
    upload_(b.sigma, ph->sigma, "krome_sigma");
    upload_(b.sigmaE, ph->sigmaE, "krome_sigmaE");
    upload_(b.sigmaEs, ph->sigmaEs, "krome_sigmaEs");
    upload_(b.sigmaCs, ph->sigmaCs, "krome_sigmaCs");
    upload_(b.sigmaPs, ph->sigmaPs, "krome_sigmaPs");
    upload_(b.habw, ph->habw, "krome_habw");
    upload_(b.elo, ph->Elo, "krome_elo");
    upload_(b.ehi, ph->Ehi, "krome_ehi");
    upload_(b.lww, ph->lww, "krome_lww");
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (nb > 0) {
    const DustBins& du = *o.dust;
    upload_(b.dxw, du.xw, "krome_dust_xw");
    upload_i_(b.disC, du.isC, "krome_dust_isC");
    // the intBB tables (nb nT each): uploaded when the population's
    // tables change
    const long ntab = static_cast<long>(nb) * du.nT;
    if (D.dtab_key != du.tab_.get() || D.dtab_n != ntab) {
      upload_(b.dibb, du.tab_->ibb, "krome_dust_ibb");
      upload_(b.dibbs, du.tab_->ibb_sig, "krome_dust_ibb_sig");
      D.dtab_key = du.tab_.get();
      D.dtab_n = ntab;
    }
    upload_(b.dqmid, du.qmid, "krome_dust_qmid");
    upload_(b.demid, du.emid, "krome_dust_emid");
    b.dnT = du.nT;
    b.dB = du.qmid.empty() ? 0 : static_cast<int>(du.emid.size());
    b.dTmin = du.TbbMin;
    b.dTmax = du.TbbMax;
    b.dTmult = du.TbbMult;
  }
#endif
}

// the per-cell ThermCtx pointers of the first ncell cells, on the device:
// the parameters' state hooks, the photo views and the static field's bin
// sums, the dust views (the contexts themselves already in b.tc)
inline void setup_ctx_(const Batch& b, int B, int ncell) {
  {
    const Batch bb = b;
    const int Bb = B;
    Kokkos::parallel_for(
        "krome_setup", Kokkos::RangePolicy<ExecSpace>(0, ncell),
        KOKKOS_LAMBDA(const int c) {
          ThermCtx& t = bb.tc(c);
          t.p.hk.n_global = &bb.n0(c, 0);
#ifdef KROME_PHOTO_HPP_INCLUDED
          if (Bb > 0) {
            t.photo_on = true;
            PhotoView v;
            v.B = Bb;
            v.sigma = bb.sigma.data();
            v.sigmaE = bb.sigmaE.data();
            v.sigmaEs = bb.sigmaEs.data();
            v.sigmaCs = bb.sigmaCs.data();
            v.sigmaPs = bb.sigmaPs.data();
            v.habw = bb.habw.data();
            v.Elo = bb.elo.data();
            v.Ehi = bb.ehi.data();
            v.lww = bb.lww.data();
            t.ph = v;
            t.Phi = &bb.phi(c, 0);
            t.p.hk.phi = t.Phi;
            t.p.hk.nbins = Bb;
            t.p.hk.elo = v.Elo;
            t.p.hk.ehi = v.Ehi;
            if (bb.Bev == 0) {
              // the static field's bin sums, once per solve (the host
              // solve_caches_)
              double* s = &bb.psum(c, 0);
              photo_sums(v, t.Phi, s, s + NPHOTO, s + 2 * NPHOTO,
                         s + 3 * NPHOTO, s + 4 * NPHOTO);
              PhotoSumsView& pv = bb.psv(c);
              bind_psum_(s, pv);
              pv.hab = photo_habing(v, t.Phi);
              t.psum = &pv;
            }
          }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
          if (bb.nb > 0) {
            DustView v;
            v.nb = bb.nb;
            v.nT = bb.dnT;
            v.TbbMin = bb.dTmin;
            v.TbbMax = bb.dTmax;
            v.TbbMult = bb.dTmult;
            v.xw = bb.dxw.data();
            v.isC = bb.disC.data();
            v.a = &bb.adr(c, 0);
            v.ibb = bb.dibb.data();
            v.ibb_sig = bb.dibbs.data();
            v.B = bb.dB;
            v.qmid = bb.dqmid.data();
            v.emid = bb.demid.data();
            t.dust_on = true;
            t.dust = v;
            t.dust_Td = &bb.tdr(c, 0);
            t.p.hk.dust.nb = bb.nb;
            t.p.hk.dust.Td = t.dust_Td;
            t.p.hk.dust.isC = v.isC;
          }
#endif
          (void)Bb;
        });
  }
}

// the per-cell device data of a call (densities, dt, fixed Tgas, ThermCtx
// with its photo and dust views, the dust rows and tables) in D.b, and the
// state layout; B: the photon bins. 0, or the host solve()'s input flags
// (init_ctx_): -102 (dust bins required), -103 (a photo module without a
// field, or radiation.evolve without one), -106 (dust.evolve with conserve
// 2), -107 (surface.bins without its bins)
inline int setup_batch_(DeviceBuffers_& D, int& B, int ncell,
                        const double* n, const double* Tgas,
                        const double* dt, const Params* P, int pstride,
                        const SolverOpts& o, const CellsOpts& co) {
  Shape_ sh;
  B = 0;
  const int sflag = batch_shape_(o, co, sh);
  if (sflag != 0) return sflag;
  B = sh.B;
  const int nb = sh.nb;
#ifdef KROME_PHOTO_HPP_INCLUDED
  const PhotoModule* ph = o.photo;
#endif
  set_shape_(D, sh, ncell);
  Batch& b = D.b;
  for (int c = 0; c < ncell; ++c) {
    D.dth(c) = dt[c];
    D.Tfh(c) = Tgas[c];
    for (int i = 0; i < NSPEC; ++i)
      D.n0h(c, i) = n[static_cast<size_t>(c) * NSPEC + i];
  }
  Kokkos::deep_copy(rows1_(b.dt, ncell), rows1_(D.dth, ncell));
  Kokkos::deep_copy(rows1_(b.Tfix, ncell), rows1_(D.Tfh, ncell));
  Kokkos::deep_copy(rows2_(b.n0, ncell), rows2_(D.n0h, ncell));

  upload_tables_(D, o, B, nb);
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (ph != nullptr) {
    for (int c = 0; c < ncell; ++c)
      for (int q = 0; q < B; ++q)
        D.phh(c, q) = co.Phi != nullptr
            ? co.Phi[static_cast<size_t>(c) * B + q] : o.Phi[q];
    Kokkos::deep_copy(rows2_(b.phi, ncell), rows2_(D.phh, ncell));
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (nb > 0) {
    const DustBins& du = *o.dust;
    for (int c = 0; c < ncell; ++c)
      for (int k = 0; k < nb; ++k) {
        const size_t r = static_cast<size_t>(c) * nb + k;
        D.tdh(c, k) = co.Td != nullptr ? co.Td[r] : du.Td[k];
        D.adh(c, k) = co.adust != nullptr ? co.adust[r] : du.a[k];
      }
    Kokkos::deep_copy(rows2_(b.tdr, ncell), rows2_(D.tdh, ncell));
    Kokkos::deep_copy(rows2_(b.adr, ncell), rows2_(D.adh, ncell));
  }
#endif
  // per-cell ThermCtx on the host, then its device-side pointers
  for (int c = 0; c < ncell; ++c) {
    ThermCtx& t = D.tch(c);
    t = ThermCtx();
    t.p = P[static_cast<size_t>(c) * pstride];
    t.p.hk = StateHooks();
    for (int r = 0; r <= NPHOTO; ++r) t.shield[r] = 1.0;
    for (int q = 0; q <= NLW; ++q) t.khat[q] = 0.0;
  }
  Kokkos::deep_copy(rows1_(b.tc, ncell), rows1_(D.tch, ncell));
  setup_ctx_(b, B, ncell);
  (void)nb;
  return 0;
}

// ---- per-cell: one BDF integration per thread --------------------
//
// Memory layout, two kernels:
// - contiguous (the default): each cell's integrator arrays (the state,
//   the Nordsieck array, the work vectors, the Newton and saved Jacobian
//   matrices, the pivots) in one block, and the rhs / Jacobian work on the
//   cell's own context (b.tc(c), its static-field sums) and on its row of
//   the rhs work (D.scr): no thread-local copies. With few warps resident
//   (255 registers) each thread's lines stay in L1.
// - interleaved (large calls): element i of cell c at c + i ncell, so the
//   32 threads of a warp reading element i read 32 consecutive doubles
//   (coalesced); bdf::Cell runs on stride ncell. Each thread works on
//   thread-local copies of its ThermCtx, of its static-field sums, of the
//   state the rhs / Jacobian read, and on a thread-local rhs work array
//   (GPU local memory: interleaved by the hardware, so these accesses
//   coalesce too), sized at compile time for at most
//   KROME_DEVICE_MAX_PHOTO_BINS co-evolved photon bins (default 32) and
//   KROME_DEVICE_MAX_DUST_BINS dust bins (default 40): a call with more
//   runs the contiguous kernel.
// The integration is the same in both (bitwise without FMA contraction).
//
// Launch: each layout's kernel has its own Kokkos::LaunchBounds<MAXT,
// MINB> (at most MAXT threads per block, at least MINB blocks per SM: on
// CUDA a register cap of 65536 / (MAXT MINB) per thread). On CUDA the
// defaults are the best of the A100 measurements (docs/gpu-percell-
// design.md, launch bounds and layout): contiguous 32,1 (more blocks
// below one wave: all SMs busy), interleaved 128,8 (64 registers: a
// large grid resident at once, where coalescing pays). Elsewhere: Kokkos'
// choice. -DKROME_PERCELL_LAUNCH_CONTIG=MAXT,MINB and
// -DKROME_PERCELL_LAUNCH_INTER=MAXT,MINB set either,
// -DKROME_PERCELL_LAUNCH=MAXT,MINB both.
//
// The layout of a call: CellsOpts::min_interleave >= 0: interleaved from
// that many cells, contiguous below. Default (< 0): contiguous; on CUDA
// interleaved above one wave of the contiguous kernel (the cells resident
// at once at 255 registers per thread: 256 threads per SM of 64K
// registers, 27 648 on an A100). -DKROME_PERCELL_INTERLEAVE=0/1 fixes the
// default layout instead (=1: the OpenMP test of the interleaved kernel).
#ifdef KROME_PERCELL_INTERLEAVE
constexpr bool PERCELL_IL_ = KROME_PERCELL_INTERLEAVE != 0;
#else
constexpr bool PERCELL_IL_ = false;
#endif
#ifdef KROME_PERCELL_LAUNCH
#define KROME_PERCELL_LC_ KROME_PERCELL_LAUNCH
#define KROME_PERCELL_LI_ KROME_PERCELL_LAUNCH
#endif
#ifdef KROME_PERCELL_LAUNCH_CONTIG
#undef KROME_PERCELL_LC_
#define KROME_PERCELL_LC_ KROME_PERCELL_LAUNCH_CONTIG
#endif
#ifdef KROME_PERCELL_LAUNCH_INTER
#undef KROME_PERCELL_LI_
#define KROME_PERCELL_LI_ KROME_PERCELL_LAUNCH_INTER
#endif
#ifdef KOKKOS_ENABLE_CUDA
#ifndef KROME_PERCELL_LC_
#define KROME_PERCELL_LC_ 32, 1
#endif
#ifndef KROME_PERCELL_LI_
#define KROME_PERCELL_LI_ 128, 8
#endif
#endif
#ifdef KROME_PERCELL_LC_
using PercellPolicyC_ = Kokkos::RangePolicy<
    ExecSpace, Kokkos::LaunchBounds<KROME_PERCELL_LC_>>;
#else
using PercellPolicyC_ = Kokkos::RangePolicy<ExecSpace>;
#endif
#ifdef KROME_PERCELL_LI_
using PercellPolicyI_ = Kokkos::RangePolicy<
    ExecSpace, Kokkos::LaunchBounds<KROME_PERCELL_LI_>>;
#else
using PercellPolicyI_ = Kokkos::RangePolicy<ExecSpace>;
#endif
template <bool IL>
using PercellPolicy_ =
    typename std::conditional<IL, PercellPolicyI_, PercellPolicyC_>::type;

// the cells the contiguous kernel holds resident at once on the current
// CUDA device (SMs x the threads of 255 registers an SM's register file
// holds); 0 when unknown or not CUDA. Read once per process (the first
// device used)
inline long percell_wave_() {
#if defined(KOKKOS_ENABLE_CUDA) && !defined(KROME_PERCELL_INTERLEAVE)
  static const long w = [] {
    int dev = 0, nsm = 0, regs = 0;
    if (cudaGetDevice(&dev) != cudaSuccess
        || cudaDeviceGetAttribute(&nsm, cudaDevAttrMultiProcessorCount, dev)
               != cudaSuccess
        || cudaDeviceGetAttribute(
               &regs, cudaDevAttrMaxRegistersPerMultiprocessor, dev)
               != cudaSuccess)
      return 0L;
    return static_cast<long>(nsm) * (regs / 256);
  }();
  return w;
#else
  return 0;
#endif
}
#ifndef KROME_DEVICE_MAX_PHOTO_BINS
#define KROME_DEVICE_MAX_PHOTO_BINS 32
#endif
#ifndef KROME_DEVICE_MAX_DUST_BINS
#define KROME_DEVICE_MAX_DUST_BINS 40
#endif
// the interleaved kernel's thread-local sizes
constexpr int BCAP_ = EvolveConfig::rad ? KROME_DEVICE_MAX_PHOTO_BINS : 0;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
constexpr int NBCAP_ = KROME_DEVICE_MAX_DUST_BINS;
#else
constexpr int NBCAP_ = 0;
#endif
constexpr int DCAP_ = EvolveConfig::Td ? NBCAP_ : 0;
constexpr int ACAP_ = EvolveConfig::dust ? NBCAP_ : 0;
constexpr int NEQCAP_ = NSPEC + 1 + BCAP_ + DCAP_ + ACAP_;
constexpr long RWCAP_ =
    rhs_work_doubles(NEQCAP_, BCAP_, NBCAP_, DCAP_, ACAP_) + 2L * NEQCAP_;

// the problem of one cell for bdf::Cell: the host's rhs and Jacobian
// (full_rhs_, jac_cell_ over the cell's RhsCtx) in its own time (no
// rescaling), with the host rhs_host_ conventions (a trial Tgas or Tdust
// <= 0, or a non-finite rhs: recoverable). Interleaved (IL): the
// integrator's vectors have stride S and the functions run on a
// thread-local copy of the state; contiguous: on the vectors themselves.
// yp, fp: the Jacobian's contiguous scratch (NEQ each)
template <bool IL>
struct CellProb_ {
  RhsCtx* c;
  long S;
  double* yp;
  double* fp;
  double yl[IL ? NEQCAP_ : 1], fl[IL ? NEQCAP_ : 1];
  KOKKOS_INLINE_FUNCTION int rhs(double, const double* y, double* f) {
    const int n = c->NEQ;
    const double* ys = y;
    double* fs = f;
    if constexpr (IL) {
      for (int i = 0; i < n; ++i) yl[i] = y[i * S];
      ys = yl;
      fs = fl;
    }
    if (!rhs_trial_ok_(ys, c)) return 1;
    full_rhs_(ys, fs, c);
    int bad = 0;
    for (int i = 0; i < n; ++i) {
      if constexpr (IL) f[i * S] = fl[i];
      if (!Kokkos::isfinite(fs[i])) bad = 1;
    }
    return bad;
  }
  // J column-major with stride S (element (i, j) at J[(i + j NEQ) S]),
  // zeroed by the caller
  KOKKOS_INLINE_FUNCTION int jac(double, const double* y, const double*,
                                 double* J, double*, double*, double*) {
    const int n = c->NEQ;
    const double* ys = y;
    if constexpr (IL) {
      for (int i = 0; i < n; ++i) yl[i] = y[i * S];
      ys = yl;
    }
    return jac_cell_(ys, DenseMat_(J, S, static_cast<long>(n) * S), c, yp,
                     fp);
  }
};

// the per-cell context of cell c (the host's init_ctx_ layout)
KOKKOS_INLINE_FUNCTION RhsCtx percell_ctx_(const Batch& b, long c) {
  RhsCtx rc;
  rc.Tfix = b.Tfix(c);
  rc.NEQ = b.blk;
  rc.iT = b.iT;
  rc.iu = b.iu;
  rc.B = b.Bev;
  rc.id = b.id;
  rc.D = b.D;
  rc.ia = b.ia;
  rc.A = b.A;
  return rc;
}

// cell c's context tc and rhs work w (the host's init_ctx_, install_hooks_
// and pack_state_): the state's slots beyond n and Tgas in yc (stride S);
// returns the Jacobian's scratch (2 blk doubles after the rhs work)
KOKKOS_INLINE_FUNCTION double* percell_begin_(const Batch& b, long c, long S,
                                              double* yc, RhsCtx& rc,
                                              ThermCtx& tc, double* w) {
  rc.th = &tc;
  double* yp = rhs_bind_(&rc, w, b.nb, b.balance);
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (rc.B > 0) {  // the evolved field starts at the cell's row
    for (int q = 0; q < rc.B; ++q) {
      rc.Phi[q] = b.phi(c, q);
      yc[(rc.iu + q) * S] = Kokkos::log(Kokkos::fmax(rc.Phi[q], 1e-300));
    }
    tc.Phi = rc.Phi;
    tc.p.hk.phi = tc.Phi;
  }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (b.nb > 0) {
    dust_area_weights(tc.dust, rc.awcur);  // the row's radii
    tc.p.hk.dust.aw = rc.awcur;
    for (int k = 0; k < rc.A; ++k) rc.arad[k] = tc.dust.a[k];
    for (int k = 0; k < rc.D; ++k) yc[(rc.id + k) * S] = b.tdr(c, k);
#ifdef KROME_HAS_PAR_d2g
    if (rc.A > 0) {  // the grain counts of the solve start
      rc.xds = rho_(&b.n0(c, 0)) * tc.p.d2g;
      tc.dust.xds_fix = rc.xds;
      tc.dust.acur = rc.arad;
      for (int k = 0; k < rc.A; ++k)
        yc[(rc.ia + k) * S] =
            dust_atoms_of_radius(tc.dust, k, tc.dust.a[k], rc.xds);
    }
#endif
  }
#endif
  return yp;
}

// after the integration of cell c (flag fl): a non-positive final Tgas is
// a failure; on success the outputs beyond n and Tgas (the evolved field,
// Td, the grain radii, the balance Td) from the final state (yc, stride S;
// yp: blk doubles of scratch). Returns the cell's flag
KOKKOS_INLINE_FUNCTION int percell_end_(const Batch& b, long c, long S,
                                        const double* yc, RhsCtx& rc,
                                        ThermCtx& tc, double* yp, int fl) {
  if (fl == 0 && rc.iT >= 0) {
    const double Tc = yc[rc.iT * S];
    if (!(Tc > 0.0 && Kokkos::isfinite(Tc))) fl = FLAG_NONPOSITIVE_T;
  }
  if (fl != 0) return fl;
  for (int i = 0; i < b.blk; ++i) yp[i] = yc[i * S];
#ifdef KROME_PHOTO_HPP_INCLUDED
  for (int q = 0; q < rc.B; ++q) b.phi(c, q) = Kokkos::exp(yp[rc.iu + q]);
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (b.nb > 0) {
    for (int k = 0; k < rc.D; ++k) b.tdr(c, k) = yp[rc.id + k];
    for (int k = 0; k < rc.A; ++k)
      b.adr(c, k) = dust_radius_of_atoms(tc.dust, k, yp[rc.ia + k], rc.xds);
#ifdef KROME_HAS_PAR_d2g
    if (rc.tdb != nullptr) {  // the final state's balance Td
      if (rc.iu >= 0) unpack_rad_(yp, &rc);
      const double* phi = nullptr;
#ifdef KROME_PHOTO_HPP_INCLUDED
      phi = tc.Phi;
#endif
      dust_balance(tc.dust, state_T_(yp, &rc), yp, tc.p, phi, &b.tdr(c, 0));
    }
#endif
  }
#endif
  (void)tc;
  return 0;
}

// cells re-solved on the host by the last per-cell call (failures on
// the device), for diagnostics
inline int percell_host_cells_ = 0;

// the per-cell arrays and the tolerances of a call (both data paths)
inline void percell_prepare_(DeviceBuffers_& D, const SolverOpts& o) {
  const Batch& b = D.b;
  const int blk = b.blk;
  const bool sparse = percell_sparse_(o, blk);
  const long wd = bdf::work_doubles(blk, sparse, NSPEC, SPARSE_NNZ);
  ensure_percell_(D, wd);
  // tolerances, the same for every cell (the host solve()'s)
  std::vector<double> at(blk), rt(blk, o.rtol);
  for (int i = 0; i < NSPEC; ++i) {
    at[i] = species_tol_(o.atol_species, i, o.atol_n);
    rt[i] = species_tol_(o.rtol_species, i, o.rtol);
  }
  if (b.iT >= 0) at[b.iT] = o.atol_T;
  for (int k = 0; k < b.Bev; ++k) at[b.iu + k] = o.atol_u;
  for (int k = 0; k < b.D; ++k) at[b.id + k] = o.atol_T;
  for (int k = 0; k < b.A; ++k) at[b.ia + k] = o.atol_n;
  upload_(D.atol, at, "krome_atol");
  upload_(D.rtolv, rt, "krome_rtolv");
}

// the per-cell kernel over the first ncell cells (their states in D.y,
// laid out as IL; percell_prepare_ done): each thread integrates its
// cell, its flag and statistics into D.flags, D.stats
template <bool IL>
inline void percell_launch_(DeviceBuffers_& D, int ncell,
                            const SolverOpts& o) {
  const Batch& b = D.b;
  const int blk = b.blk;
  const bool sparse = percell_sparse_(o, blk);
  const long wd = bdf::work_doubles(blk, sparse, NSPEC, SPARSE_NNZ);
  const long S = IL ? ncell : 1;
  const DView y = D.y, work = D.work, scr = D.scr;
  const IView1 iwork = D.iwork;
  const IView1 flags = D.flags;
  const LView2 stats = D.stats;
  const DView atol = D.atol;
  const bool varies = rtol_varies_(o);
  const DView rtolv = D.rtolv;
  const JacMode jm = o.jac != JAC_AUTO ? o.jac
                     : (o.evolve_T ? JAC_FD : JAC_ANALYTIC);
  bdf::Opts bo;
  bo.rtol = o.rtol;
  bo.max_steps = o.max_steps;
  bo.max_err_fails = o.max_err_fails;
  bo.max_conv_fails = o.max_conv_fails;
  bo.dq_jac = jm == JAC_FD;
  bo.sparse = sparse;  // (the tables: set in the kernel, sparse_tables_)
  const int max_restarts = o.max_restarts;
  const int fail_every = device_fail_every_;
  Kokkos::parallel_for(
      "krome_bdf_cells", PercellPolicy_<IL>(0, ncell),
      KOKKOS_LAMBDA(const int c) {
        int f = 0;
        bdf::Stats st;
        int nr = 0;
        const double dtc = b.dt(c);
        if (dtc > 0.0) {
          double* yc = y.data() + (IL ? c : static_cast<long>(c) * blk);
          RhsCtx rc = percell_ctx_(b, c);
          // integrate the cell with the context tc and the rhs work w
          auto run = [&](ThermCtx& tc, double* w) {
            double* yp = percell_begin_(b, c, S, yc, rc, tc, w);
            double* fp = yp + blk;
            CellProb_<IL> pr;
            pr.c = &rc;
            pr.S = S;
            pr.yp = yp;
            pr.fp = fp;
            bdf::Opts oc = bo;
            oc.atol = atol.data();
            oc.rtol_v = varies ? rtolv.data() : nullptr;
            if (oc.sparse) oc.sp = sparse_tables_();
            charge_newton_(oc);
            double* wc = work.data() + (IL ? c : static_cast<long>(c) * wd);
            int* ic = iwork.data() + (IL ? c : static_cast<long>(c) * blk);
            int fl = bdf::integrate<CellProb_<IL>, !IL>(
                pr, blk, yc, dtc, oc, max_restarts, wc, ic, &st, &nr, S);
            if (fail_every > 0 && c % fail_every == 0)
              fl = FLAG_TOO_MUCH_WORK;  // tests: device_fail_every_
            return percell_end_(b, c, S, yc, rc, tc, yp, fl);
          };
          // a runtime if, not if constexpr: nvcc rejects an extended lambda
          // that first-captures a variable in a constexpr-if context (IL is
          // a template constant, so the branch is folded anyway)
          if (IL) {
            // thread-local copies of the cell's context and work (GPU
            // local memory: interleaved by the hardware)
            ThermCtx tc = b.tc(c);
#ifdef KROME_PHOTO_HPP_INCLUDED
            double psl[NPSUM_ * (NPHOTO > 0 ? NPHOTO : 1)];
            PhotoSumsView pvl;
            if (tc.psum != nullptr) {
              const double* src = &b.psum(c, 0);
              for (int i = 0; i < NPSUM_ * NPHOTO; ++i) psl[i] = src[i];
              bind_psum_(psl, pvl);
              pvl.hab = tc.psum->hab;
              tc.psum = &pvl;
            }
#endif
            double wl[RWCAP_];
            f = run(tc, wl);
          } else {
            // the cell's own context and rhs work row in global memory
            f = run(b.tc(c), scr.data() + static_cast<long>(c) * b.rw);
          }
        }
        flags(c) = f;
        stats(c, ST_STEPS_) = st.nst;
        stats(c, ST_RHS_) = st.nfe;
        stats(c, ST_JAC_) = st.nje;
        stats(c, ST_RHS_JAC_) = st.nfeDQ;
        stats(c, ST_LSETUPS_) = st.nsetups;
        stats(c, ST_NEWTON_) = st.nni;
        stats(c, ST_NEWTON_FAILS_) = st.nnf;
        stats(c, ST_ERR_FAILS_) = st.netf;
        stats(c, ST_RESTARTS_) = nr;
      });
  Kokkos::fence();
}

template <bool IL>
inline int solve_percell_layout_(DeviceBuffers_& D, int ncell, double* n,
                                 double* Tgas, const double* dt,
                                 const Params* P, int pstride,
                                 const SolverOpts& o, const CellsOpts& co) {
  percell_prepare_(D, o);
  const Batch& b = D.b;
  const int blk = b.blk;
  const DView y = D.y;
  auto yh = D.yh;
  auto yix = [&](int c, int i) {
    return IL ? static_cast<size_t>(i) * ncell + c
              : static_cast<size_t>(c) * blk + i;
  };
  // n and Tgas; the kernel packs the other slots (field, Td, grains)
  for (int c = 0; c < ncell; ++c) {
    for (int i = 0; i < NSPEC; ++i)
      yh(yix(c, i)) = n[static_cast<size_t>(c) * NSPEC + i];
    if (b.iT >= 0) yh(yix(c, b.iT)) = Tgas[c];
  }
  const int ny = ncell * blk;
  Kokkos::deep_copy(rows1_(y, ny), rows1_(yh, ny));
  percell_launch_<IL>(D, ncell, o);
  auto fh = D.fh;
  auto sh = D.sh;
  Kokkos::deep_copy(rows1_(yh, ny), rows1_(y, ny));
  Kokkos::deep_copy(rows1_(fh, ncell), rows1_(D.flags, ncell));
  Kokkos::deep_copy(rows2_(sh, ncell), rows2_(D.stats, ncell));
  // the per-cell rows the caller gets back
  const int B = b.Bev;
  const bool phi_out = B > 0 && co.Phi != nullptr;
  const bool td_out = b.nb > 0 && co.Td != nullptr;
  const bool a_out = b.A > 0 && co.adust != nullptr;
  if (phi_out) Kokkos::deep_copy(rows2_(D.phh, ncell), rows2_(b.phi, ncell));
  if (td_out) Kokkos::deep_copy(rows2_(D.tdh, ncell), rows2_(b.tdr, ncell));
  if (a_out) Kokkos::deep_copy(rows2_(D.adh, ncell), rows2_(b.adr, ncell));
  std::vector<int> failed;
  std::vector<int> fl(ncell, 0);
  for (int c = 0; c < ncell; ++c) {
    if (fh(c) != 0) {
      failed.push_back(c);
      continue;
    }
    if (co.stats != nullptr) {
      SolveStats& s = co.stats[c];
      s = SolveStats();
      s.steps = sh(c, ST_STEPS_);
      s.rhs_evals = sh(c, ST_RHS_);
      s.jac_evals = sh(c, ST_JAC_);
      s.rhs_evals_jac = sh(c, ST_RHS_JAC_);
      s.lin_setups = sh(c, ST_LSETUPS_);
      s.newton_iters = sh(c, ST_NEWTON_);
      s.newton_fails = sh(c, ST_NEWTON_FAILS_);
      s.err_fails = sh(c, ST_ERR_FAILS_);
      s.restarts = static_cast<int>(sh(c, ST_RESTARTS_));
    }
    if (!(dt[c] > 0.0)) continue;
    double n0c[NSPEC];
    double* nc = n + static_cast<size_t>(c) * NSPEC;
    for (int i = 0; i < NSPEC; ++i) {
      n0c[i] = nc[i];
      nc[i] = yh(yix(c, i));
      if (o.clip_negative && nc[i] < 0.0) nc[i] = 0.0;
    }
    conserve_cell_(n0c, nc, P[static_cast<size_t>(c) * pstride], o.conserve);
    if (b.iT >= 0) Tgas[c] = yh(yix(c, b.iT));
    for (int q = 0; phi_out && q < B; ++q)
      co.Phi[static_cast<size_t>(c) * B + q] = D.phh(c, q);
    for (int k = 0; td_out && k < b.nb; ++k)
      co.Td[static_cast<size_t>(c) * b.nb + k] = D.tdh(c, k);
    for (int k = 0; a_out && k < b.nb; ++k)
      co.adust[static_cast<size_t>(c) * b.nb + k] = D.adh(c, k);
  }
  percell_host_cells_ = static_cast<int>(failed.size());
  if (!failed.empty()) {  // the host from the untouched inputs
    const int nf = static_cast<int>(failed.size());
    int Bp = 0;
#ifdef KROME_PHOTO_HPP_INCLUDED
    Bp = o.photo != nullptr ? o.photo->B : 0;
#endif
    const int nb = b.nb;
    std::vector<double> nn(static_cast<size_t>(nf) * NSPEC), TT(nf), dd(nf),
        ph(co.Phi != nullptr ? static_cast<size_t>(nf) * Bp : 0),
        td(co.Td != nullptr ? static_cast<size_t>(nf) * nb : 0),
        ad(co.adust != nullptr ? static_cast<size_t>(nf) * nb : 0);
    std::vector<Params> pp(pstride != 0 ? nf : 1, P[0]);
    for (int k = 0; k < nf; ++k) {
      const int c = failed[k];
      for (int i = 0; i < NSPEC; ++i)
        nn[static_cast<size_t>(k) * NSPEC + i] =
            n[static_cast<size_t>(c) * NSPEC + i];
      TT[k] = Tgas[c];
      dd[k] = dt[c];
      if (pstride != 0) pp[k] = P[c];
      for (int q = 0; !ph.empty() && q < Bp; ++q)
        ph[static_cast<size_t>(k) * Bp + q] =
            co.Phi[static_cast<size_t>(c) * Bp + q];
      for (int q = 0; !td.empty() && q < nb; ++q)
        td[static_cast<size_t>(k) * nb + q] =
            co.Td[static_cast<size_t>(c) * nb + q];
      for (int q = 0; !ad.empty() && q < nb; ++q)
        ad[static_cast<size_t>(k) * nb + q] =
            co.adust[static_cast<size_t>(c) * nb + q];
    }
    std::vector<int> fk(nf, 0);
    std::vector<SolveStats> sk(nf);
    CellsOpts ck = co;
    ck.path = PATH_PERCELL;
    ck.flags = fk.data();
    ck.stats = sk.data();
    ck.Phi = co.Phi != nullptr ? ph.data() : nullptr;
    ck.Td = co.Td != nullptr ? td.data() : nullptr;
    ck.adust = co.adust != nullptr ? ad.data() : nullptr;
    solve_cells_percell(nf, nn.data(), TT.data(), dd.data(), pp.data(),
                        pstride != 0 ? 1 : 0, o, ck);
    for (int k = 0; k < nf; ++k) {
      const int c = failed[k];
      fl[c] = fk[k];
      if (co.stats != nullptr) co.stats[c] = sk[k];
      if (fk[k] != 0) continue;
      for (int i = 0; i < NSPEC; ++i)
        n[static_cast<size_t>(c) * NSPEC + i] =
            nn[static_cast<size_t>(k) * NSPEC + i];
      Tgas[c] = TT[k];
      for (int q = 0; !ph.empty() && q < Bp; ++q)
        co.Phi[static_cast<size_t>(c) * Bp + q] =
            ph[static_cast<size_t>(k) * Bp + q];
      for (int q = 0; !td.empty() && q < nb; ++q)
        co.Td[static_cast<size_t>(c) * nb + q] =
            td[static_cast<size_t>(k) * nb + q];
      for (int q = 0; !ad.empty() && q < nb; ++q)
        co.adust[static_cast<size_t>(c) * nb + q] =
            ad[static_cast<size_t>(k) * nb + q];
    }
  }
  int first = 0;
  for (int c = 0; c < ncell; ++c) {
    if (co.flags != nullptr) co.flags[c] = fl[c];
    if (first == 0) first = fl[c];
  }
  return first;
}

// the layout of a call (or chunk) of ncell cells: CellsOpts::min_interleave
// when >= 0, else on CUDA interleaved above one wave of the contiguous
// kernel, elsewhere PERCELL_IL_ (contiguous unless forced); contiguous
// when the interleaved kernel's thread-local arrays are too small
inline bool percell_interleave_(int ncell, const CellsOpts& co,
                                const Batch& b) {
  if (b.Bev > BCAP_ || b.nb > NBCAP_) return false;
  if (co.min_interleave >= 0) return ncell >= co.min_interleave;
  const long w = percell_wave_();
  return w > 0 ? ncell > w : PERCELL_IL_;
}

// the per-cell kernel in the layout percell_interleave_ picks
inline int solve_percell_(DeviceBuffers_& D, int ncell, double* n,
                          double* Tgas, const double* dt, const Params* P,
                          int pstride, const SolverOpts& o,
                          const CellsOpts& co) {
  const bool il = percell_interleave_(ncell, co, D.b);
  return il ? solve_percell_layout_<true>(D, ncell, n, Tgas, dt, P, pstride,
                                          o, co)
            : solve_percell_layout_<false>(D, ncell, n, Tgas, dt, P,
                                           pstride, o, co);
}

// one chunk of cells on the device with the buffers D (the checks of
// solve_cells_kokkos done)
inline int solve_chunk_(DeviceBuffers_& D, int ncell, double* n,
                        double* Tgas, const double* dt, const Params* P,
                        int pstride, const SolverOpts& o,
                        const CellsOpts& co) {
  percell_host_cells_ = 0;
  int B = 0;
  const int sflag = setup_batch_(D, B, ncell, n, Tgas, dt, P, pstride, o, co);
  if (sflag != 0) return sflag;
  return solve_percell_(D, ncell, n, Tgas, dt, P, pstride, o, co);
}

// the device buffers of a call: the cache (CellsOpts::device_cache) when
// no other host thread holds it (lk, locked here, keeps it until the call
// returns; the cache is freed before Kokkos finalizes), else `local`,
// allocated for this call alone
inline DeviceBuffers_* acquire_buffers_(const CellsOpts& co,
                                        DeviceBuffers_& local,
                                        std::unique_lock<std::mutex>& lk) {
  if (!co.device_cache || !lk.try_lock()) return &local;
  if (device_cache_ == nullptr) device_cache_ = new DeviceBuffers_();
  if (!device_cache_hook_) {
    Kokkos::push_finalize_hook([] {
      std::lock_guard<std::mutex> g(device_cache_mutex_);
      if (device_cache_ != nullptr) *device_cache_ = DeviceBuffers_();
      device_cache_hook_ = false;
    });
    device_cache_hook_ = true;
  }
  return device_cache_;
}

// the default memory budget of a call: this fraction of the device's free
// memory (plus the buffers the call reuses)
constexpr double DEVICE_MEMORY_FRACTION_ = 0.8;

// the cells of one chunk: ncell, or as many as the memory budget holds
// (CellsOpts::device_memory bytes; 0: DEVICE_MEMORY_FRACTION_ of the free
// memory, where the runtime reports it; < 0: unbounded); 0 when not even
// one cell fits
inline int chunk_cells_(int ncell, const DeviceBuffers_& D,
                        const SolverOpts& o, const CellsOpts& co) {
  long long budget = co.device_memory;
  if (budget == 0) {
    const long long f = device_free_bytes_();
    budget = f < 0 ? -1
                   : static_cast<long long>(
                         DEVICE_MEMORY_FRACTION_ * static_cast<double>(
                             f + static_cast<long long>(buffers_bytes_(D))));
  }
  const long long per = percell_bytes_(o, co);
  if (budget > 0 && per > 0)
    return static_cast<int>(std::min<long long>(ncell, budget / per));
  return ncell;
}

// integrate ncell cells over their own dt[c] on the device, one thread
// per cell; P: one Params (pstride 0) or ncell of them (pstride 1). Calls
// of more cells than fit in one go are integrated in chunks on the same
// buffers: at most the cells that CellsOpts::device_memory holds (by
// default 80 % of the device's free memory, the buffers this call can
// reuse included; CUDA and HIP; chunk_cells_). An allocation that still
// fails (out of memory) frees the buffers and halves the chunk; when not
// even one cell fits, the remaining cells run on the host with PATH_AUTO,
// and get
// FLAG_DEVICE_MEMORY (-111) otherwise. Returns 0, the host's input flags
// (-102, -103, -105, -106, -107; every cell gets it), -109 (Kokkos not
// initialized), -111, or the flag of the first cell that failed on the
// device and again on the host (the others written).
inline int solve_cells_kokkos(int ncell, double* n, double* Tgas,
                              const double* dt, const Params* P,
                              int pstride, const SolverOpts& o,
                              const CellsOpts& co) {
  if (ncell <= 0) return 0;
  if (!Kokkos::is_initialized()) return FLAG_NO_KOKKOS;
  auto fail_all = [&](int f) {
    if (co.flags != nullptr)
      for (int c = 0; c < ncell; ++c) co.flags[c] = f;
    return f;
  };
#ifdef KROME_GAS_TOTAL
  if (o.conserve != 0) return fail_all(FLAG_CONSERVE_GAS_TOTAL);
#endif
  for (int c = 0; c < ncell; ++c)
    if (params_missing(P[static_cast<size_t>(c) * pstride]) >= 0)
      return fail_all(FLAG_PARAM_UNSET);
  // the buffers: the cache, unless off or in use by another caller
  DeviceBuffers_ local;
  std::unique_lock<std::mutex> lk(device_cache_mutex_, std::defer_lock);
  DeviceBuffers_* D = acquire_buffers_(co, local, lk);
  int chunk = chunk_cells_(ncell, *D, o, co);
  int B = 0, nb = 0;
#ifdef KROME_PHOTO_HPP_INCLUDED
  B = o.photo != nullptr ? o.photo->B : 0;
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  nb = o.dust != nullptr ? o.dust->nb : 0;
#endif
  // with radiation.evolve a single cell without a Phi row evolves opts.Phi
  CellsOpts c1 = co;
  if (EvolveConfig::rad && ncell == 1 && co.Phi == nullptr) c1.Phi = o.Phi;
  int first = 0, hosted = 0;
  for (int c0 = 0; c0 < ncell;) {
    const int m = std::min(chunk, ncell - c0);
    CellsOpts cc = c1;
    if (c1.Phi != nullptr) cc.Phi = c1.Phi + static_cast<size_t>(c0) * B;
    if (co.Td != nullptr) cc.Td = co.Td + static_cast<size_t>(c0) * nb;
    if (co.adust != nullptr)
      cc.adust = co.adust + static_cast<size_t>(c0) * nb;
    if (co.flags != nullptr) cc.flags = co.flags + c0;
    if (co.stats != nullptr) cc.stats = co.stats + c0;
    double* nc = n + static_cast<size_t>(c0) * NSPEC;
    const Params* Pc = P + static_cast<size_t>(c0) * pstride;
    int f = 0;
    bool oom = m < 1;
    if (!oom) {
      // an allocation (or launch) out of memory throws before any output
      // is written: the chunk's inputs are untouched
      try {
        f = solve_chunk_(*D, m, nc, Tgas + c0, dt + c0, Pc, pstride, o, cc);
      } catch (const std::exception& e) {
        if (!out_of_memory_(e)) throw;
        oom = true;
      }
    }
    if (oom) {
      *D = DeviceBuffers_();  // whatever the attempt allocated
      if (m > 1) {            // half the chunk, from the same cell
        chunk = m / 2;
        continue;
      }
      // not even one cell fits: the remaining cells
      const int r = ncell - c0;
      if (co.path == PATH_AUTO) {
        cc.path = PATH_PERCELL;
        f = solve_cells_percell(r, nc, Tgas + c0, dt + c0, Pc, pstride, o,
                                cc);
        hosted += r;
      } else {
        if (co.flags != nullptr)
          for (int c = c0; c < ncell; ++c) co.flags[c] = FLAG_DEVICE_MEMORY;
        f = FLAG_DEVICE_MEMORY;
      }
      if (first == 0) first = f;
      break;
    }
    hosted += percell_host_cells_;
    if (input_flag_(f)) return fail_all(f);  // the inputs, for every cell
    if (first == 0) first = f;
    c0 += m;
  }
  percell_host_cells_ = hosted;
  return first;
}

// ---- device-resident data (level A; gpu-percell-design.md 14.3) ---------
//
// The cells' arrays already in the device's memory (the Kokkos memory
// space): the call reads and writes them there, and only the cells that
// fail on the device travel to the host (re-solved by solve(), flagged).
// KROME-X's own conventions: cgs, Tgas, the network's species and
// parameters; the host code gathers its cells into these arrays first.

// the species arrays: one array per species, n[i * ld + c] (a Fortran
// x(ncell, nspec), the default), or one row per cell, n[c * NSPEC + i]
constexpr int LAYOUT_SPECIES_MAJOR = 0;
constexpr int LAYOUT_CELL_MAJOR = 1;

// the caller's device arrays of a call
struct DeviceCells {
  double* n = nullptr;          // number densities [cm^-3], or the mass
                                // fractions when rho is given
  const double* rho = nullptr;  // [ncell] gas density [g cm^-3] (x input)
  int layout = LAYOUT_SPECIES_MAJOR;
  long ld = 0;                  // species-major: the stride between the
                                // species' arrays (0: ncell)
  double* Tgas = nullptr;       // [ncell]
  const double* dt = nullptr;   // [ncell] [s]
  const double* par = nullptr;  // parameter rows (param_names order;
                                // nullptr: the defaults)
  int par_stride = 0;           // doubles between rows: 0 (one row) or NPAR
  double* Phi = nullptr;        // [ncell][photo bins], optional
  double* Td = nullptr;         // [ncell][dust bins], optional
  double* adust = nullptr;      // [ncell][dust bins], optional
  int* flags = nullptr;         // [ncell], optional: 0 solved on the device,
                                // FLAG_HOST_RECOVERED, or < 0 (failed)
};

// a call's outcome, on the host
struct DeviceSummary {
  long device = 0;  // cells solved on the device
  long host = 0;    // failed on the device, recovered on the host
  long failed = 0;  // failed on both (left untouched)
};

KROME_HD inline long dev_idx_(int layout, long ld, long c, int i) {
  return layout == LAYOUT_CELL_MAJOR ? c * NSPEC + i : i * ld + c;
}

// one chunk of m cells (a: the chunk's arrays, already offset); the flags
// of the failed cells and the summary updated; returns 0, an input flag
// for every cell (-102, -103, -106, -107), or the lowest failing cell's
// flag
inline int solve_chunk_device_(DeviceBuffers_& D, int m, const DeviceCells& a,
                               const SolverOpts& o, const CellsOpts& co,
                               DeviceSummary& sum) {
  CellsOpts cs = co;
  cs.Phi = a.Phi;
  Shape_ sh;
  const int sflag = batch_shape_(o, cs, sh);
  if (sflag != 0) return sflag;
  const int B = sh.B, nb = sh.nb;
  set_shape_(D, sh, m);
  Batch& b = D.b;
  upload_tables_(D, o, B, nb);
  // the rows' defaults: the options' field, the population's Td and radii
  DView phi0, td0, a0, xin("krome_xin", m);
  std::vector<double> hphi, htd, ha;
#ifdef KROME_PHOTO_HPP_INCLUDED
  if (B > 0 && o.Phi != nullptr) hphi.assign(o.Phi, o.Phi + B);
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  if (nb > 0) {
    htd = o.dust->Td;
    ha = o.dust->a;
  }
#endif
  upload_(phi0, hphi, "krome_phi0");
  upload_(td0, htd, "krome_td0");
  upload_(a0, ha, "krome_a0");
  const long ld = a.ld > 0 ? a.ld : m;
  const int lay = a.layout;
  const DeviceCells ac = a;
  const int ps = a.par_stride;
  // gather: the cells' inputs and contexts into the buffers
  Kokkos::parallel_for(
      "krome_gather", Kokkos::RangePolicy<ExecSpace>(0, m),
      KOKKOS_LAMBDA(const int c) {
        b.dt(c) = ac.dt[c];
        b.Tfix(c) = ac.Tgas[c];
        double xs = 0.0;
        for (int i = 0; i < NSPEC; ++i) {
          double v = ac.n[dev_idx_(lay, ld, c, i)];
          if (ac.rho != nullptr) {
            if (species_mass_g[i] > 0.0) xs += v;
            v = species_mass_g[i] > 0.0 ? ac.rho[c] * v / species_mass_g[i]
                                        : 0.0;
          }
          b.n0(c, i) = v;
        }
        xin(c) = xs;
        for (int q = 0; q < B; ++q)
          b.phi(c, q) = ac.Phi != nullptr
              ? ac.Phi[static_cast<long>(c) * B + q] : phi0(q);
        for (int k = 0; k < nb; ++k) {
          b.tdr(c, k) = ac.Td != nullptr
              ? ac.Td[static_cast<long>(c) * nb + k] : td0(k);
          b.adr(c, k) = ac.adust != nullptr
              ? ac.adust[static_cast<long>(c) * nb + k] : a0(k);
        }
        ThermCtx& t = b.tc(c);
        t = ThermCtx();
        if (ac.par != nullptr)
          params_from_row(ac.par + static_cast<long>(c) * ps, t.p);
        t.p.hk = StateHooks();
        for (int r = 0; r <= NPHOTO; ++r) t.shield[r] = 1.0;
        for (int q = 0; q <= NLW; ++q) t.khat[q] = 0.0;
      });
  setup_ctx_(b, B, m);
  percell_prepare_(D, o);
  const bool il = percell_interleave_(m, co, b);
  const int blk = b.blk;
  const DView y = D.y;
  {  // the states, in the kernel's layout: element i of cell c at
     // y[c cstride + i S] (interleaved: S = m cells, cstride 1;
     // contiguous: S 1, cstride blk)
    const long S = il ? m : 1;
    const long cstride = il ? 1 : blk;
    Kokkos::parallel_for(
        "krome_y0", Kokkos::RangePolicy<ExecSpace>(0, m),
        KOKKOS_LAMBDA(const int c) {
          for (int i = 0; i < NSPEC; ++i) y(c * cstride + i * S) = b.n0(c, i);
          if (b.iT >= 0) y(c * cstride + b.iT * S) = b.Tfix(c);
        });
  }
  if (il)
    percell_launch_<true>(D, m, o);
  else
    percell_launch_<false>(D, m, o);
  const IView1 flags = D.flags;
  // scatter: the solved cells back into the caller's arrays
  {
    const long S = il ? m : 1;
    const long cstride = il ? 1 : blk;
    const bool clip = o.clip_negative, renorm = o.renormalize_x;
    const int cons = o.conserve;
    Kokkos::parallel_for(
        "krome_scatter", Kokkos::RangePolicy<ExecSpace>(0, m),
        KOKKOS_LAMBDA(const int c) {
          const int f = flags(c);
          if (ac.flags != nullptr) ac.flags[c] = f;
          if (f != 0 || !(b.dt(c) > 0.0)) return;
          double nc[NSPEC];
          for (int i = 0; i < NSPEC; ++i) {
            double v = y(c * cstride + i * S);
            if (clip && v < 0.0) v = 0.0;
            nc[i] = v;
          }
          conserve_cell_(&b.n0(c, 0), nc, b.tc(c).p, cons);
          if (ac.rho != nullptr) {
            const double r = ac.rho[c];
            double xs = 0.0;
            for (int i = 0; i < NSPEC; ++i) {
              nc[i] = species_mass_g[i] * nc[i] / r;
              if (species_mass_g[i] > 0.0) xs += nc[i];
            }
            if (renorm && xs > 0.0)
              for (int i = 0; i < NSPEC; ++i) nc[i] *= xin(c) / xs;
          }
          for (int i = 0; i < NSPEC; ++i) ac.n[dev_idx_(lay, ld, c, i)] = nc[i];
          if (b.iT >= 0) ac.Tgas[c] = y(c * cstride + b.iT * S);
          if (ac.Phi != nullptr)
            for (int q = 0; q < b.Bev; ++q)
              ac.Phi[static_cast<long>(c) * B + q] = b.phi(c, q);
          if (ac.Td != nullptr)
            for (int k = 0; k < b.nb; ++k)
              ac.Td[static_cast<long>(c) * nb + k] = b.tdr(c, k);
          if (ac.adust != nullptr && b.A > 0)
            for (int k = 0; k < b.nb; ++k)
              ac.adust[static_cast<long>(c) * nb + k] = b.adr(c, k);
        });
  }
  Kokkos::fence();
  // the failed cells: their inputs, still in the buffers, to the host
  auto fh = D.fh;
  Kokkos::deep_copy(rows1_(fh, m), rows1_(flags, m));
  std::vector<int> failed;
  for (int c = 0; c < m; ++c)
    if (fh(c) != 0) failed.push_back(c);
  sum.device += m - static_cast<long>(failed.size());
  if (failed.empty()) return 0;
  const int nf = static_cast<int>(failed.size());
  IView1 fidx("krome_failed", nf);
  {
    Kokkos::View<const int*, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>> hv(failed.data(),
                                                             nf);
    Kokkos::deep_copy(fidx, hv);
  }
  // a failed cell's inputs, one row of st: the columns
  const int IN_N = 0, IN_T = NSPEC, IN_DT = IN_T + 1, IN_RHO = IN_DT + 1,
            IN_PAR = IN_RHO + 1, IN_PHI = IN_PAR + NPAR, IN_TD = IN_PHI + B,
            IN_A = IN_TD + nb, IN_W = IN_A + nb;
  DView2 st("krome_failed_cells", nf, IN_W);
  Kokkos::parallel_for(
      "krome_pack_failed", Kokkos::RangePolicy<ExecSpace>(0, nf),
      KOKKOS_LAMBDA(const int k) {
        const int c = fidx(k);
        for (int i = 0; i < NSPEC; ++i) st(k, IN_N + i) = b.n0(c, i);
        st(k, IN_T) = b.Tfix(c);
        st(k, IN_DT) = b.dt(c);
        st(k, IN_RHO) = ac.rho != nullptr ? ac.rho[c] : 0.0;
        for (int i = 0; i < NPAR; ++i)
          st(k, IN_PAR + i) = ac.par != nullptr
              ? ac.par[static_cast<long>(c) * ps + i] : 0.0;
        for (int q = 0; q < B; ++q) st(k, IN_PHI + q) = b.phi(c, q);
        for (int q = 0; q < nb; ++q) st(k, IN_TD + q) = b.tdr(c, q);
        for (int q = 0; q < nb; ++q) st(k, IN_A + q) = b.adr(c, q);
      });
  auto sth = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), st);
  auto xinh = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), xin);
  std::vector<double> nn(static_cast<size_t>(nf) * NSPEC), TT(nf), dd(nf),
      ph(static_cast<size_t>(nf) * B), td(static_cast<size_t>(nf) * nb),
      ad(static_cast<size_t>(nf) * nb);
  std::vector<Params> pp(nf);
  for (int k = 0; k < nf; ++k) {
    const size_t kn = static_cast<size_t>(k) * NSPEC,
                 kB = static_cast<size_t>(k) * B,
                 kb = static_cast<size_t>(k) * nb;
    for (int i = 0; i < NSPEC; ++i) nn[kn + i] = sth(k, IN_N + i);
    TT[k] = sth(k, IN_T);
    dd[k] = sth(k, IN_DT);
    double row[NPAR > 0 ? NPAR : 1];
    for (int i = 0; i < NPAR; ++i) row[i] = sth(k, IN_PAR + i);
    if (a.par != nullptr) params_from_row(row, pp[k]);
    for (int q = 0; q < B; ++q) ph[kB + q] = sth(k, IN_PHI + q);
    for (int q = 0; q < nb; ++q) td[kb + q] = sth(k, IN_TD + q);
    for (int q = 0; q < nb; ++q) ad[kb + q] = sth(k, IN_A + q);
  }
  std::vector<int> fk(nf, 0);
  CellsOpts ck = co;
  ck.path = PATH_PERCELL;
  ck.flags = fk.data();
  ck.stats = nullptr;
  ck.Phi = B > 0 ? ph.data() : nullptr;
  ck.Td = nb > 0 ? td.data() : nullptr;
  ck.adust = nb > 0 ? ad.data() : nullptr;
  solve_cells_percell(nf, nn.data(), TT.data(), dd.data(), pp.data(), 1, o,
                      ck);
  // the results back to the device and into the caller's arrays: one row
  // of ot per failed cell, the columns
  const int OUT_FLAG = 0, OUT_N = 1, OUT_T = OUT_N + NSPEC,
            OUT_PHI = OUT_T + 1, OUT_TD = OUT_PHI + B, OUT_A = OUT_TD + nb,
            OUT_W = OUT_A + nb;
  DView2 ot("krome_failed_out", nf, OUT_W);
  auto oh = Kokkos::create_mirror_view(ot);
  int first = 0;
  for (int k = 0; k < nf; ++k) {
    const int f = fk[k];
    oh(k, OUT_FLAG) = f == 0 ? FLAG_HOST_RECOVERED : f;
    if (f == 0) {
      ++sum.host;
    } else {
      ++sum.failed;
      if (first == 0) first = f;
    }
    double* nk = nn.data() + static_cast<size_t>(k) * NSPEC;
    if (a.rho != nullptr && f == 0) {
      double xk[NSPEC];
      n_to_x_out_(nk, sth(k, IN_RHO), xinh(failed[k]), o, xk);
      for (int i = 0; i < NSPEC; ++i) nk[i] = xk[i];
    }
    const size_t kB = static_cast<size_t>(k) * B,
                 kb = static_cast<size_t>(k) * nb;
    for (int i = 0; i < NSPEC; ++i) oh(k, OUT_N + i) = nk[i];
    oh(k, OUT_T) = TT[k];
    for (int q = 0; q < B; ++q) oh(k, OUT_PHI + q) = ph[kB + q];
    for (int q = 0; q < nb; ++q) oh(k, OUT_TD + q) = td[kb + q];
    for (int q = 0; q < nb; ++q) oh(k, OUT_A + q) = ad[kb + q];
  }
  Kokkos::deep_copy(ot, oh);
  Kokkos::parallel_for(
      "krome_unpack_failed", Kokkos::RangePolicy<ExecSpace>(0, nf),
      KOKKOS_LAMBDA(const int k) {
        const int c = fidx(k);
        const int f = static_cast<int>(ot(k, OUT_FLAG));
        if (ac.flags != nullptr) ac.flags[c] = f;
        if (f != FLAG_HOST_RECOVERED) return;
        for (int i = 0; i < NSPEC; ++i)
          ac.n[dev_idx_(lay, ld, c, i)] = ot(k, OUT_N + i);
        if (b.iT >= 0) ac.Tgas[c] = ot(k, OUT_T);
        if (ac.Phi != nullptr && b.Bev > 0)
          for (int q = 0; q < B; ++q)
            ac.Phi[static_cast<long>(c) * B + q] = ot(k, OUT_PHI + q);
        if (ac.Td != nullptr)
          for (int q = 0; q < nb; ++q)
            ac.Td[static_cast<long>(c) * nb + q] = ot(k, OUT_TD + q);
        if (ac.adust != nullptr && b.A > 0)
          for (int q = 0; q < nb; ++q)
            ac.adust[static_cast<long>(c) * nb + q] = ot(k, OUT_A + q);
      });
  Kokkos::fence();
  percell_host_cells_ += nf;
  return first;
}

// the device's flags of cells [0, ncell) set to f
inline void device_fill_flags_(int* flags, int ncell, int f) {
  if (flags == nullptr || ncell <= 0) return;
  Kokkos::parallel_for(
      "krome_fill_flags", Kokkos::RangePolicy<ExecSpace>(0, ncell),
      KOKKOS_LAMBDA(const int c) { flags[c] = f; });
  Kokkos::fence();
}

// Integrate ncell cells whose arrays are in the device's memory (a), each
// over its own dt, as solve_cells_kokkos does with host arrays: the same
// kernel, layouts, buffer cache, memory budget and chunking (CellsOpts:
// device_cache, device_memory, min_interleave; its Phi, Td, adust, flags
// and stats are not used: a's rows and flags instead). The results are
// written in place, on the device (fenced before and after: the caller's
// kernels are complete, its next ones see the results); the cells that
// fail on the device are re-solved on the host from their inputs and
// flagged FLAG_HOST_RECOVERED, or keep the host's (negative) flag and
// their inputs. sum (optional, on
// the host): the counts. Returns 0, the inputs' flag for every cell (-102,
// -103, -105 a REQUIRED parameter unset in some row, -106, -107), -109
// (Kokkos not initialized), -111 (not even one cell fits in the memory),
// or the lowest failing cell's flag
inline int solve_cells_device(int ncell, const DeviceCells& a,
                              const SolverOpts& o, const CellsOpts& co,
                              DeviceSummary* sum = nullptr) {
  DeviceSummary local_sum;
  DeviceSummary& sm = sum != nullptr ? *sum : local_sum;
  sm = DeviceSummary();
  if (ncell <= 0) return 0;
  if (!Kokkos::is_initialized()) return FLAG_NO_KOKKOS;
  Kokkos::fence();  // the caller's kernels on the arrays are complete
  auto fail_all = [&](int f) {
    device_fill_flags_(a.flags, ncell, f);
    sm.failed = ncell;
    return f;
  };
#ifdef KROME_GAS_TOTAL
  if (o.conserve != 0) return fail_all(FLAG_CONSERVE_GAS_TOTAL);
#endif
  {  // a REQUIRED parameter unset in any row: nothing is integrated
    const DeviceCells ac = a;
    const int ps = a.par_stride;
    int missing = 0;
    Kokkos::parallel_reduce(
        "krome_params_check", Kokkos::RangePolicy<ExecSpace>(0, ncell),
        KOKKOS_LAMBDA(const int c, int& acc) {
          Params p;
          if (ac.par != nullptr)
            params_from_row(ac.par + static_cast<long>(c) * ps, p);
          if (params_missing(p) >= 0) acc += 1;
        },
        missing);
    if (missing > 0) return fail_all(FLAG_PARAM_UNSET);
  }
  DeviceBuffers_ local;
  std::unique_lock<std::mutex> lk(device_cache_mutex_, std::defer_lock);
  DeviceBuffers_* D = acquire_buffers_(co, local, lk);
  CellsOpts cs = co;
  cs.Phi = a.Phi;
  int chunk = chunk_cells_(ncell, *D, o, cs);
  int B = 0, nb = 0;
#ifdef KROME_PHOTO_HPP_INCLUDED
  B = o.photo != nullptr ? o.photo->B : 0;
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
  nb = o.dust != nullptr ? o.dust->nb : 0;
#endif
  percell_host_cells_ = 0;
  int first = 0;
  const long ld = a.ld > 0 ? a.ld : ncell;
  for (int c0 = 0; c0 < ncell;) {
    const int m = std::min(chunk, ncell - c0);
    DeviceCells ac = a;  // the chunk's arrays
    ac.ld = ld;
    ac.n = a.layout == LAYOUT_CELL_MAJOR
        ? a.n + static_cast<size_t>(c0) * NSPEC : a.n + c0;
    if (a.rho != nullptr) ac.rho = a.rho + c0;
    ac.Tgas = a.Tgas + c0;
    ac.dt = a.dt + c0;
    if (a.par != nullptr)
      ac.par = a.par + static_cast<size_t>(c0) * a.par_stride;
    if (a.Phi != nullptr) ac.Phi = a.Phi + static_cast<size_t>(c0) * B;
    if (a.Td != nullptr) ac.Td = a.Td + static_cast<size_t>(c0) * nb;
    if (a.adust != nullptr) ac.adust = a.adust + static_cast<size_t>(c0) * nb;
    if (a.flags != nullptr) ac.flags = a.flags + c0;
    int f = 0;
    bool oom = m < 1;
    if (!oom) {
      try {
        f = solve_chunk_device_(*D, m, ac, o, cs, sm);
      } catch (const std::exception& e) {
        if (!out_of_memory_(e)) throw;
        oom = true;
      }
    }
    if (oom) {
      *D = DeviceBuffers_();
      if (m > 1) {
        chunk = m / 2;
        continue;
      }
      device_fill_flags_(ac.flags, ncell - c0, FLAG_DEVICE_MEMORY);
      sm.failed += ncell - c0;
      if (first == 0) first = FLAG_DEVICE_MEMORY;
      break;
    }
    if (input_flag_(f)) return fail_all(f);
    if (first == 0) first = f;
    c0 += m;
  }
  return first;
}

// consistent advection (conserve_to_n_; conserve_cells) on cells whose
// arrays are in the device's memory: n (or the mass fractions x, mass true;
// then ref per gram) in either layout, the reference rows ref (ref_stride
// doubles between rows: 0 or NELEM), the parameter rows par (par_stride:
// 0 or NPAR; nullptr: the defaults), flags[ncell] (optional). Returns the
// flag of the lowest failing cell (negative), else 0
inline int conserve_cells_device(int ncell, double* v, bool mass, int layout,
                                 long ld, const double* ref, int ref_stride,
                                 const double* par, int par_stride,
                                 int* flags) {
  if (ncell <= 0) return 0;
  if (!Kokkos::is_initialized()) return FLAG_NO_KOKKOS;
  Kokkos::fence();
  const long l = ld > 0 ? ld : ncell;
  IView1 own;  // the flags, when the caller wants none
  int* fp = flags;
  if (fp == nullptr) {
    own = IView1("krome_conserve_flags", ncell);
    fp = own.data();
  }
  int lowest = ncell;  // the lowest failing cell
  Kokkos::parallel_reduce(
      "krome_conserve_cells", Kokkos::RangePolicy<ExecSpace>(0, ncell),
      KOKKOS_LAMBDA(const int c, int& acc) {
        double x[NSPEC];
        for (int i = 0; i < NSPEC; ++i) x[i] = v[dev_idx_(layout, l, c, i)];
        Params p;
        if (par != nullptr)
          params_from_row(par + static_cast<long>(c) * par_stride, p);
        const double* r = ref + static_cast<long>(c) * ref_stride;
        const int f = mass ? conserve_to_x_(x, r, p) : conserve_to_(x, r, p);
        for (int i = 0; i < NSPEC; ++i) v[dev_idx_(layout, l, c, i)] = x[i];
        fp[c] = f;
        if (f < 0 && c < acc) acc = c;
      },
      Kokkos::Min<int>(lowest));
  Kokkos::fence();
  if (lowest >= ncell) return 0;
  int f = 0;
  Kokkos::deep_copy(f, Kokkos::View<int, MemSpace, Kokkos::MemoryTraits<
                           Kokkos::Unmanaged>>(fp + lowest));
  return f;
}

// registration with the solve_cells dispatcher (solve_cells.hpp)
inline const bool registered_ =
    (device_cells_ = &solve_cells_kokkos, true);

}  // namespace kk
}  // namespace krome
