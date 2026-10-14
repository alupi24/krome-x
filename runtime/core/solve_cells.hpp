// KROME-X cells frontend (Step 3): integrate many independent resolution
// elements through one call,
//
//   solve_cells(ncell, n[ncell][NSPEC], Tgas[ncell], dt[ncell], params,
//               opts, cells_opts)
//
// params is either one Params broadcast to every cell, or an array of
// ncell Params: in a hydro coupling most "parameters" are cell-dependent
// (crate, d2g, Av, tff, total_Z, opratio, Tdust, ...), so both backends
// read cell c's parameters as params[c] when the per-cell overload is
// used.
//
// Two backends behind one dispatcher:
//
//  PER-CELL   an (OpenMP-parallel) loop over cells, each running the
//             ordinary solve() with fully independent step control (the
//             modern analog of legacy KROME's per-cell DLSODES).
//  DEVICE     (solve_kokkos.hpp, -DKROME_USE_KOKKOS) all cells in one call
//             to the Kokkos execution space, one thread per cell running
//             the same integrator (core/cvode/bdf_cell.hpp).

// The dispatcher (PATH_AUTO) routes calls of CellsOpts::min_device cells
// and more to the device backend when one is compiled in, the others to
// PER-CELL.

// Implementation header: include krome.hpp (the public API), which
// declares solve_cells(), solve_cells_x() and CellsOpts.
#pragma once
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace krome {


// ---------------------------------------------------------------------------
// per-cell backend

inline int solve_cells_percell(int ncell, double* n, double* Tgas,
                               const double* dt, const Params* P,
                               int pstride, const SolverOpts& o,
                               const CellsOpts& co) {
  std::vector<int> fl(static_cast<size_t>(std::max(ncell, 0)), 0);
#ifdef KROME_PHOTO_HPP_INCLUDED
  const int B = o.photo ? o.photo->B : 0;
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int c = 0; c < ncell; ++c) {
    SolverOpts oc = o;
#ifdef KROME_PHOTO_HPP_INCLUDED
    // cell c's own field slice: solve() only reads a static field and
    // updates an evolved one (radiation.evolve) in place. Without rows, an
    // evolved field is opts.Phi for a single cell, and a private copy of
    // it per cell otherwise (the cells must not update one field together)
    std::vector<double> phic;
    if (co.Phi != nullptr && o.photo != nullptr) {
      oc.Phi = co.Phi + static_cast<size_t>(c) * B;
    } else if (EvolveConfig::rad && o.photo != nullptr && o.Phi != nullptr
               && ncell > 1) {
      phic.assign(o.Phi, o.Phi + B);
      oc.Phi = phic.data();
    }
#endif
#ifdef KROME_DUST_BINS_HPP_INCLUDED
    // binned dust carries per-cell state (Td): give each cell its own
    // copy, seeded from the caller's bins
    DustBins dustc;
    double* tdc = nullptr;  // cell c's row of co.Td
    double* adc = nullptr;  // cell c's row of co.adust (dust.evolve)
    if (o.dust != nullptr) {
      dustc = *o.dust;
      oc.dust = &dustc;
      if (co.Td != nullptr) {
        tdc = co.Td + static_cast<size_t>(c) * dustc.nb;
        for (int b = 0; b < dustc.nb; ++b) dustc.Td[b] = tdc[b];
      }
      if (co.adust != nullptr) {
        adc = co.adust + static_cast<size_t>(c) * dustc.nb;
        for (int b = 0; b < dustc.nb; ++b) dustc.a[b] = adc[b];
        dustc.set_radii(nullptr);  // the area weights of these radii
      }
    }
#endif
    report_cell_ = c;  // named in an error report
    report_ncell_ = ncell;
    fl[c] = dt[c] > 0.0
        ? solve(n + static_cast<size_t>(c) * NSPEC, Tgas[c],
                P[static_cast<size_t>(c) * pstride], dt[c], oc,
                co.stats != nullptr ? co.stats + c : nullptr)
        : 0;
    if (co.stats != nullptr && !(dt[c] > 0.0)) co.stats[c] = SolveStats();
    report_cell_ = -1;
#ifdef KROME_DUST_BINS_HPP_INCLUDED
    if (tdc != nullptr && fl[c] == 0 && dt[c] > 0.0)
      for (int b = 0; b < dustc.nb; ++b) tdc[b] = dustc.Td[b];
    if (adc != nullptr && fl[c] == 0 && dt[c] > 0.0)
      for (int b = 0; b < dustc.nb; ++b) adc[b] = dustc.a[b];
#endif
  }
  // the returned flag is that of the first failing cell (deterministic
  // whatever the thread schedule)
  int first = 0;
  for (int c = 0; c < ncell; ++c) {
    if (co.flags != nullptr) co.flags[c] = fl[c];
    if (first == 0) first = fl[c];
  }
  return first;
}

// ---------------------------------------------------------------------------
// dispatcher

// the device backend, when one is compiled in (solve_kokkos.hpp registers
// it): same arguments as solve_cells_percell, 0 or a flag
using DeviceCellsFn = int (*)(int, double*, double*, const double*,
                              const Params*, int, const SolverOpts&,
                              const CellsOpts&);
inline DeviceCellsFn device_cells_ = nullptr;

inline int solve_cells_(int ncell, double* n, double* Tgas,
                        const double* dt, const Params* P, int pstride,
                        const SolverOpts& o, const CellsOpts& co) {
  // a REQUIRED parameter unset in any cell: nothing is integrated
  bool missing = false;
  for (int c = 0; c < ncell; ++c) {
    const bool m = params_missing(P[static_cast<size_t>(c) * pstride]) >= 0;
    if (co.flags != nullptr) co.flags[c] = m ? -105 : 0;
    missing |= m;
  }
  if (missing) return -105;
  BatchPath path = co.path;
  const bool automatic = path == PATH_AUTO;
  if (automatic)
    path = device_cells_ != nullptr && ncell >= co.min_device ? PATH_DEVICE
                                                              : PATH_PERCELL;
  if (path == PATH_DEVICE) {
    if (device_cells_ == nullptr) {  // no device backend in this build
      if (co.flags != nullptr)
        for (int c = 0; c < ncell; ++c) co.flags[c] = -104;
      return -104;
    }
    const int f = device_cells_(ncell, n, Tgas, dt, P, pstride, o, co);
    // -104: a configuration the device backend does not support; solved
    // per cell on the host instead. -109 (Kokkos not initialized): with
    // PATH_AUTO the cells run on the host too, while an explicit
    // PATH_DEVICE reports it. Otherwise final (the cells that failed on
    // the device were already re-solved on the host)
    if (f != -104 && !(f == -109 && automatic)) return f;
  }
  return solve_cells_percell(ncell, n, Tgas, dt, P, pstride, o, co);
}

// one Params broadcast to every cell
inline int solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                       const Params& p, const SolverOpts& o,
                       const CellsOpts& co) {
  return solve_cells_(ncell, n, Tgas, dt, &p, 0, o, co);
}

// per-cell parameters: params[ncell]
inline int solve_cells(int ncell, double* n, double* Tgas, const double* dt,
                       const Params* params, const SolverOpts& o,
                       const CellsOpts& co) {
  return solve_cells_(ncell, n, Tgas, dt, params, 1, o, co);
}

// mass fractions x[ncell][NSPEC] with the per-cell gas density rho[ncell]
// [g cm^-3] (see solve_x): converted on the way in and, for the cells that
// succeed, out; a failed cell keeps its input x
inline int solve_cells_x_(int ncell, double* x, const double* rho,
                          double* Tgas, const double* dt, const Params* P,
                          int pstride, const SolverOpts& o,
                          const CellsOpts& co) {
  if (ncell <= 0) return 0;
  std::vector<double> n(static_cast<size_t>(ncell) * NSPEC), xin(ncell);
  for (int c = 0; c < ncell; ++c) {
    xin[c] = x_sum_(x + static_cast<size_t>(c) * NSPEC);
    x_to_n_(x + static_cast<size_t>(c) * NSPEC, rho[c],
            n.data() + static_cast<size_t>(c) * NSPEC);
  }
  std::vector<int> fl(ncell, 0);
  CellsOpts cc = co;
  cc.flags = co.flags != nullptr ? co.flags : fl.data();
  const int r = solve_cells_(ncell, n.data(), Tgas, dt, P, pstride, o, cc);
  for (int c = 0; c < ncell; ++c)
    if (cc.flags[c] == 0) {
      double* xc = x + static_cast<size_t>(c) * NSPEC;
      n_to_x_out_(n.data() + static_cast<size_t>(c) * NSPEC, rho[c],
                  xin[c], o, xc);
    }
  return r;
}
inline int solve_cells_x(int ncell, double* x, const double* rho,
                         double* Tgas, const double* dt, const Params& p,
                         const SolverOpts& o,
                         const CellsOpts& co) {
  return solve_cells_x_(ncell, x, rho, Tgas, dt, &p, 0, o, co);
}
inline int solve_cells_x(int ncell, double* x, const double* rho,
                         double* Tgas, const double* dt,
                         const Params* params, const SolverOpts& o,
                         const CellsOpts& co) {
  return solve_cells_x_(ncell, x, rho, Tgas, dt, params, 1, o, co);
}

} // namespace krome
