// KROME-X: a per-cell BDF integrator, ported from CVODE (SUNDIALS 7.9.0)
// for the per-cell device path (docs/gpu-percell-design.md).
//
// One cell, one integrator, no allocation, no library: the state lives in a
// caller-provided workspace (work_doubles(N) doubles, N ints) and every
// function is KROME_HD, so the same code runs per cell on the host and in a
// device kernel (one thread per cell). The algorithm is CVODE's, transcribed
// function by function (the CVODE names are kept in the comments): the
// fixed-leading-coefficient BDF in Nordsieck form, orders 1-5, the initial
// step (cvHin), the step and order selection, the local error test, the
// modified Newton iteration with its convergence test (the SUNDIALS Newton
// module, at most 3 iterations), the Jacobian / Newton-matrix reuse
// heuristics of the linear-solver interface (msbp 20, msbj 51, dgmax 0.3 /
// 0.2) with the gamma-ratio scaling of the correction, CVODE's dense
// difference-quotient Jacobian, SUNDIALS' dense LU with partial pivoting,
// and the interpolated output at tout (CV_NORMAL). The vector operations
// reproduce the special cases of the serial N_Vector (N_VLinearSum's
// branches), so with the same rhs and Jacobian the integration repeats
// CVODE's operations in CVODE's order.
//
// Not ported (unused by KROME): Adams methods, stability-limit detection,
// root finding, tstop, constraints, projection, sensitivities, iterative
// linear solvers, the fixed-point nonlinear solver, the fused vector
// kernels.
//
// The problem is a class with
//   KROME_HD int rhs(double t, const double* y, double* f);
//   KROME_HD int jac(double t, const double* y, const double* fy,
//                    double* J, double* t1, double* t2, double* t3);
// returning 0, > 0 (recoverable failure: CVODE retries with a smaller
// step) or < 0 (unrecoverable); J is N x N, column-major (element (i, j)
// at J[i + j N], zeroed on entry), t1..t3 are N-vector scratch. jac is not
// called with Opts::dq_jac.
//
// Opts::sparse replaces the dense LU of the Newton matrix by a sparse one
// with static pivots (KROME's addition, the legacy DLSODES way: a
// fill-reducing order fixed in advance, no pivoting): the leading ns
// equations (the species) form the sparse block, factorized row by row
// in the order and with the fill of the generated tables
// (Codegen.sparse_lu_tables); the remaining m = N - ns equations (Tgas,
// ...) are a dense border, eliminated through the Schur complement
// (dense LU with partial pivoting). A zero pivot is a recoverable setup
// failure (a smaller step brings the matrix closer to the identity).
//
// ---------------------------------------------------------------------
// Derived from SUNDIALS 7.9.0: src/cvode/cvode.c, cvode_ls.c, cvode_nls.c
// (Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban, and Dan
// Shumaker @ LLNL), src/sunnonlinsol/newton, src/sunlinsol/dense,
// src/sundials/sundials_dense.c, src/nvector/serial.
//
// SUNDIALS Copyright Start
// Copyright (c) 2025-2026, Lawrence Livermore National Security,
// University of Maryland Baltimore County, and the SUNDIALS contributors.
// Copyright (c) 2013-2025, Lawrence Livermore National Security
// and Southern Methodist University.
// Copyright (c) 2002-2013, Lawrence Livermore National Security.
// All rights reserved.
//
// See the LICENSE and NOTICE of SUNDIALS, reproduced in
// runtime/core/cvode/LICENSE.sundials.
//
// SPDX-License-Identifier: BSD-3-Clause
// SUNDIALS Copyright End
// ---------------------------------------------------------------------
#pragma once
#include <cfloat>
#include <cmath>

#ifndef KROME_HD  // as the generated network defines it
#if defined(__CUDACC__) || defined(__HIPCC__)
#define KROME_HD __host__ __device__
#else
#define KROME_HD
#endif
#endif

// LAPACK's dense LU (host builds with -DKROME_USE_LAPACK; never in device
// code): the calls of SUNDIALS' LapackDense linear solver, dgetrf in the
// setup and dgetrs in the solve, on the same column-major matrix
#if defined(KROME_USE_LAPACK) && !defined(__CUDA_ARCH__) \
    && !defined(__HIP_DEVICE_COMPILE__)
#define KROME_BDF_LAPACK_ 1
extern "C" {
void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv,
             int* info);
void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a,
             const int* lda, const int* ipiv, double* b, const int* ldb,
             int* info);
}
#endif

// benchmarking (bench/gpu_lu): -DKROME_BDF_NEWTON_HOOK=name passes every
// Newton matrix the host integrator is about to factorize (column-major,
// n x n) to `void name(const double* A, int n)`, defined by the caller;
// off by default
#ifdef KROME_BDF_NEWTON_HOOK
void KROME_BDF_NEWTON_HOOK(const double* A, int n);
#endif

namespace krome {
namespace bdf {

constexpr int QMAX = 5;  // BDF_Q_MAX (CVODE's default maxord)
constexpr int L_MAX = QMAX + 1;

// return flags: CVODE's values (cvode.h)
enum : int {
  SUCCESS = 0,
  TOO_MUCH_WORK = -1,
  TOO_MUCH_ACC = -2,
  ERR_FAILURE = -3,
  CONV_FAILURE = -4,
  LSETUP_FAIL = -6,
  LSOLVE_FAIL = -7,
  RHSFUNC_FAIL = -8,
  FIRST_RHSFUNC_ERR = -9,
  REPTD_RHSFUNC_ERR = -10,
  UNREC_RHSFUNC_ERR = -11,
  NLS_FAIL = -16,
  ILL_INPUT = -22,
  TOO_CLOSE = -27
};

// the sparse block's tables (Codegen.sparse_lu_tables), in pivot order:
// perm[k] the equation of pivot k, row k at ci[rp[k] .. rp[k+1]) (columns
// ascending; L left of the diagonal, U from it on), di[k] its diagonal
struct Sparse {
  int ns = 0, nnz = 0;
  const int *perm = nullptr, *rp = nullptr, *ci = nullptr, *di = nullptr;
};

struct Opts {
  double rtol = 1e-4;
  const double* atol = nullptr;    // N absolute tolerances
  // N relative tolerances, or nullptr (rtol for all): the weights are then
  // 1 / fma(rtol_i, |y_i|, atol_i), KROME's ewt_cb_ (CVodeWFtolerances)
  const double* rtol_v = nullptr;
  long max_steps = 500;   // per call (mxstep): 0 the default 500, < 0 none
  int max_err_fails = 7;  // within one step (maxnef)
  int max_conv_fails = 10;  // within one step (maxncf)
  double h0 = 0.0;        // first step (hin); 0: cvHin's estimate
  bool dq_jac = false;    // CVODE's difference-quotient Jacobian
  // LAPACK's dgetrf / dgetrs for the Newton systems (SUNDIALS' LapackDense);
  // only in host builds with KROME_USE_LAPACK and unit stride, else ignored
  bool lapack = false;
  // the static-pivot sparse LU (sp: the tables) instead of the dense one
  bool sparse = false;
  Sparse sp;
  // a linear invariant c . y of the rhs (c . f = 0) kept exactly by the
  // Newton iteration (KROME's addition; KROME: the charge, component the
  // electrons). With a Jacobian that satisfies c . J = 0 an exact Newton
  // correction has c . delta = -c . G (G the residual) for any gamma, so
  // c . y stays at its predicted value; the difference-quotient Jacobian
  // breaks this at roundoff over tiny increments (trace species), and the
  // error, conserved from then on, can outgrow a species the invariant
  // ties (the electrons of a collapsing cloud). With inv_row >= 0 every
  // correction's component inv_row is set so that c . delta = -c . G
  // over the first inv_n components (the Jacobian itself is untouched:
  // the fixed point is the same), while component inv_row (a density)
  // carries at least 1e-6 of sum_{i != r} |c_i y_i| and stays
  // non-negative: below that its value from the invariant would be lost
  // to cancellation (and a state whose invariant it cannot carry, KROME:
  // an input with more net charge than electrons, would need it
  // negative); there the correction stays the solver's own. inv_row < 0:
  // off (CVODE's iteration)
  int inv_row = -1, inv_n = 0;
  const int* inv_c = nullptr;
};

// CVODE's counters (CVodeGetNum*): steps, rhs evaluations (without the
// difference-quotient Jacobian's, nfeDQ), Jacobian evaluations, Newton
// matrix setups, Newton iterations and convergence failures (nni, nnf),
// step attempts failed in the nonlinear solve (ncfn) and in the error
// test (netf)
struct Stats {
  long nst = 0, nfe = 0, nfeDQ = 0, nje = 0, nsetups = 0, nni = 0, nnf = 0,
       ncfn = 0, netf = 0;
  KROME_HD void add(const Stats& s) {
    nst += s.nst;
    nfe += s.nfe;
    nfeDQ += s.nfeDQ;
    nje += s.nje;
    nsetups += s.nsetups;
    nni += s.nni;
    nnf += s.nnf;
    ncfn += s.ncfn;
    netf += s.netf;
  }
};

// the integrator's state when integrate() returned (the CVodeGet* values
// of the last CVode call: the error report's): current time, last / next
// step and order, the tolerance scale factor, that call's counters; and,
// when the arrays are given, the error weights and the estimated local
// errors (N values each)
struct Diag {
  double tcur = 0.0, hlast = 0.0, hcur = 0.0, tolsf = 1.0;
  int qlast = 0, qcur = 0;
  Stats last;
  double* ewt = nullptr;
  double* ele = nullptr;
};

// the Newton matrix's storage: N x N, or with Opts::sparse the factors
// (nnz), a row accumulator (ns), the border's columns and rows (2 ns m)
// and its Schur complement (m x m), m = N - ns
KROME_HD constexpr long newton_doubles(int N, bool sparse = false,
                                       int ns = 0, int nnz = 0) {
  return sparse ? nnz + ns + 2L * ns * (N - ns) + (long)(N - ns) * (N - ns)
                : (long)N * N;
}
// the workspace of an N-equation cell: the Nordsieck array, 9 vectors, the
// saved Jacobian (N x N) and the Newton matrix (newton_doubles)
KROME_HD constexpr long work_doubles(int N, bool sparse = false, int ns = 0,
                                     int nnz = 0) {
  return (L_MAX + 9) * (long)N + (long)N * N
         + newton_doubles(N, sparse, ns, nnz);
}
KROME_HD constexpr long work_ints(int N) { return N; }

namespace detail_ {

// internal flags (cvode_impl.h, cvode_ls_impl.h; the Newton module's
// SUN_NLS_CONTINUE / SUN_NLS_CONV_RECVR)
enum : int {
  DO_ERROR_TEST = 2,
  PREDICT_AGAIN = 3,
  TRY_AGAIN = 5,
  FIRST_CALL = 6,
  PREV_CONV_FAIL = 7,
  PREV_ERR_FAIL = 9,
  RHSFUNC_RECVR = 10,
  NLS_CONTINUE = 901,
  NLS_CONV_RECVR = 902,
  LS_LUFACT_FAIL = 1  // a positive (recoverable) linear-solver setup flag
};
enum : int { NO_FAILURES = 0, FAIL_BAD_J = 1, FAIL_OTHER = 2 };

// cvode.c / cvode_impl.h / cvode_ls_impl.h / cvode_nls.c constants
constexpr double FUZZ_FACTOR = 100.0;
constexpr double HLB_FACTOR = 100.0;
constexpr double HUB_FACTOR = 0.1;
constexpr double H_BIAS = 0.5;
constexpr int MAX_ITERS = 4;
constexpr double CORTES = 0.1;  // nlscoef
constexpr double ETA_MIN_FX = 0.0;
constexpr double ETA_MAX_FX = 1.5;
constexpr double ETA_MAX_FS = 10000.0;
constexpr double ETA_MAX_ES = 10.0;
constexpr double ETA_MAX_GS = 10.0;
constexpr double ETA_MIN = 0.1;
constexpr double ETA_MAX_EF = 0.2;
constexpr double ETA_MIN_EF = 0.1;
constexpr double ETA_CF = 0.25;
constexpr long SMALL_NST = 10;
constexpr int SMALL_NEF = 2;
constexpr double ONEPSM = 1.000001;
constexpr double ADDON = 1e-6;
constexpr double BIAS1 = 6.0;
constexpr double BIAS2 = 6.0;
constexpr double BIAS3 = 10.0;
constexpr int LONG_WAIT = 10;
constexpr int MXNEF1 = 3;
constexpr long MSBP = 20;              // MSBP_DEFAULT
constexpr double DGMAX_LSETUP = 0.3;   // DGMAX_LSETUP_DEFAULT
constexpr long MSBJ = 51;              // CVLS_MSBJ
constexpr double DGMAX_JBAD = 0.2;     // CVLS_DGMAX
constexpr double MIN_INC_MULT = 1000.0;
constexpr int NLS_MAXCOR = 3;
constexpr double CRDOWN = 0.3;
constexpr double RDIV = 2.0;

// sundials_math.h
KROME_HD inline double smax_(double a, double b) { return a > b ? a : b; }
KROME_HD inline double smin_(double a, double b) { return a < b ? a : b; }
KROME_HD inline double ssqrt_(double x) { return x <= 0.0 ? 0.0 : sqrt(x); }
KROME_HD inline double spowi_(double base, int e) {
  double prod = 1.0;
  const int ex = e < 0 ? -e : e;
  for (int i = 1; i <= ex; i++) prod *= base;
  if (e < 0) prod = 1.0 / prod;
  return prod;
}

// element i of a vector of stride S: unit stride (U, the host layout) or
// interleaved cells (element i of cell c at c + i S, the device layout)
template <bool U>
KROME_HD inline long ix_(long i, long S) {
  return U ? i : i * S;
}

// the serial N_Vector operations, with N_VLinearSum's special cases
// (each rounds differently); every vector of a call has stride S
template <bool U>
KROME_HD inline void axpy_(int N, double a, const double* x, double* y,
                           long S) {
  if (a == 1.0) {
    for (int i = 0; i < N; i++) y[ix_<U>(i, S)] += x[ix_<U>(i, S)];
    return;
  }
  if (a == -1.0) {
    for (int i = 0; i < N; i++) y[ix_<U>(i, S)] -= x[ix_<U>(i, S)];
    return;
  }
  for (int i = 0; i < N; i++) y[ix_<U>(i, S)] += a * x[ix_<U>(i, S)];
}
template <bool U>
KROME_HD inline void lin_sum_(int N, double a, const double* x, double b,
                              const double* y, double* z, long S) {
  if (b == 1.0 && z == y) {
    axpy_<U>(N, a, x, z, S);
    return;
  }
  if (a == 1.0 && z == x) {
    axpy_<U>(N, b, y, z, S);
    return;
  }
  if (a == 1.0 && b == 1.0) {
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = x[k] + y[k];
    }
    return;
  }
  bool test;
  if ((test = (a == 1.0 && b == -1.0)) || (a == -1.0 && b == 1.0)) {
    const double* v1 = test ? y : x;
    const double* v2 = test ? x : y;
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = v2[k] - v1[k];
    }
    return;
  }
  if ((test = (a == 1.0)) || b == 1.0) {
    const double c = test ? b : a;
    const double* v1 = test ? y : x;
    const double* v2 = test ? x : y;
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = (c * v1[k]) + v2[k];
    }
    return;
  }
  if ((test = (a == -1.0)) || b == -1.0) {
    const double c = test ? b : a;
    const double* v1 = test ? y : x;
    const double* v2 = test ? x : y;
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = (c * v1[k]) - v2[k];
    }
    return;
  }
  if (a == b) {
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = a * (x[k] + y[k]);
    }
    return;
  }
  if (a == -b) {
    for (int i = 0; i < N; i++) {
      const long k = ix_<U>(i, S);
      z[k] = a * (x[k] - y[k]);
    }
    return;
  }
  for (int i = 0; i < N; i++) {
    const long k = ix_<U>(i, S);
    z[k] = (a * x[k]) + (b * y[k]);
  }
}
template <bool U>
KROME_HD inline void scale_(int N, double c, const double* x, double* z,
                            long S) {
  for (int i = 0; i < N; i++) z[ix_<U>(i, S)] = c * x[ix_<U>(i, S)];
}
template <bool U>
KROME_HD inline void const_(int N, double c, double* z, long S) {
  for (int i = 0; i < N; i++) z[ix_<U>(i, S)] = c;
}
template <bool U>
KROME_HD inline double wrms_(int N, const double* x, const double* w,
                             long S) {
  double sum = 0.0;
  for (int i = 0; i < N; i++) {
    const double p = x[ix_<U>(i, S)] * w[ix_<U>(i, S)];
    sum += p * p;
  }
  return ssqrt_(sum / N);
}
template <bool U>
KROME_HD inline double max_norm_(int N, const double* x, long S) {
  double mx = 0.0;
  for (int i = 0; i < N; i++)
    if (fabs(x[ix_<U>(i, S)]) > mx) mx = fabs(x[ix_<U>(i, S)]);
  return mx;
}

// SUNDlsMat_denseGETRF / GETRS: column-major a (element (i, k) at
// a[(i + k n) S]), partial pivoting (pivot k at p[k S]); returns 0 or
// k + 1 at a zero pivot
template <bool U>
KROME_HD inline int getrf_(double* a, int n, int* p, long S) {
  for (int k = 0; k < n; k++) {
    double* col_k = a + ix_<U>((long)k * n, S);
    int l = k;
    for (int i = k + 1; i < n; i++)
      if (fabs(col_k[ix_<U>(i, S)]) > fabs(col_k[ix_<U>(l, S)])) l = i;
    p[ix_<U>(k, S)] = l;
    if (col_k[ix_<U>(l, S)] == 0.0) return k + 1;
    if (l != k) {
      for (int i = 0; i < n; i++) {
        double* ci = a + ix_<U>((long)i * n, S);
        const double temp = ci[ix_<U>(l, S)];
        ci[ix_<U>(l, S)] = ci[ix_<U>(k, S)];
        ci[ix_<U>(k, S)] = temp;
      }
    }
    const double mult = 1.0 / col_k[ix_<U>(k, S)];
    for (int i = k + 1; i < n; i++) col_k[ix_<U>(i, S)] *= mult;
    for (int j = k + 1; j < n; j++) {
      double* col_j = a + ix_<U>((long)j * n, S);
      const double a_kj = col_j[ix_<U>(k, S)];
      if (a_kj != 0.0)
        for (int i = k + 1; i < n; i++)
          col_j[ix_<U>(i, S)] -= a_kj * col_k[ix_<U>(i, S)];
    }
  }
  return 0;
}
template <bool U>
KROME_HD inline void getrs_(const double* a, int n, const int* p, double* b,
                            long S) {
  for (int k = 0; k < n; k++) {
    const int pk = p[ix_<U>(k, S)];
    if (pk != k) {
      const double tmp = b[ix_<U>(k, S)];
      b[ix_<U>(k, S)] = b[ix_<U>(pk, S)];
      b[ix_<U>(pk, S)] = tmp;
    }
  }
  for (int k = 0; k < n - 1; k++) {
    const double* col_k = a + ix_<U>((long)k * n, S);
    for (int i = k + 1; i < n; i++)
      b[ix_<U>(i, S)] -= col_k[ix_<U>(i, S)] * b[ix_<U>(k, S)];
  }
  for (int k = n - 1; k > 0; k--) {
    const double* col_k = a + ix_<U>((long)k * n, S);
    b[ix_<U>(k, S)] /= col_k[ix_<U>(k, S)];
    for (int i = 0; i < k; i++)
      b[ix_<U>(i, S)] -= col_k[ix_<U>(i, S)] * b[ix_<U>(k, S)];
  }
  b[0] /= a[0];
}

}  // namespace detail_

// The integrator of one cell (CVODE's CVodeMem, BDF, dense direct solver,
// Newton). reinit() starts an integration at (t0, y0); solve() is one
// CVode(tout, CV_NORMAL) call from there.
// U = false: every vector and matrix of the cell has stride S in its
// workspace (the interleaved device layout: element i of cell c at
// c + i S, coalesced across threads), the problem's rhs / jac see the same
// strides; U = true: unit stride (S unused), the host layout.
template <class Prob, bool U = true>
class Cell {

 public:
  Stats st;  // the counters since the last reinit()

  KROME_HD Cell(Prob& prob, int n, double* work, int* iwork, const Opts& opts,
                long stride = 1)
      : P(prob), N(n), o(opts), S(stride) {
    const long vN = ix((long)N);  // one vector, in elements of the work
    double* w = work;
    for (int j = 0; j < L_MAX; j++, w += vN) zn[j] = w;
    ewt = w;
    w += vN;
    y = w;
    w += vN;
    acor = w;
    w += vN;
    tempv = w;
    w += vN;
    ftemp = w;
    w += vN;
    delta = w;
    w += vN;
    vtemp1 = w;
    w += vN;
    vtemp2 = w;
    w += vN;
    vtemp3 = w;
    w += vN;
    A = w;
    w += ix(newton_doubles(N, o.sparse, o.sp.ns, o.sp.nnz));
    savedJ = w;
    piv = iwork;
    mxstep = o.max_steps == 0 ? 500 : o.max_steps;
    atolmin0 = false;
    if (o.atol != nullptr) {  // CVodeSVtolerances: N_VMin(abstol) == 0
      double mn = o.atol[0];
      for (int i = 1; i < N; i++)
        if (o.atol[i] < mn) mn = o.atol[i];
      atolmin0 = mn == 0.0;
    }
  }

  // CVodeReInit (+ the counters reset by cvLsInitialize / cvNlsInit),
  // with the first step hin (0: estimated)
  KROME_HD void reinit(double t0, const double* y0, double h_init) {
    tn = t0;
    for (int i = 0; i < N; i++) zn[0][ix(i)] = y0[ix(i)];
    q = 1;
    L = 2;
    qwait = L;
    etamax = detail_::ETA_MAX_FS;
    qu = 0;
    hu = 0.0;
    tolsf = 1.0;
    st = Stats();
    nstlp = 0;
    nscon = 0;
    h0u = 0.0;
    next_h = 0.0;
    next_q = 0;
    hin = h_init;
    nstlj = 0;
    acnrmcur = false;
    nls_jcur = false;
    jcur = false;
    saved_tq5 = 0.0;
    for (int i = 0; i <= L_MAX; i++) tau[i] = 0.0;
    eta = 1.0;
  }

  // CVode(tout, CV_NORMAL), once after reinit(): integrates to tout and
  // writes y(tout) (interpolated) to yout; on failure yout is the last
  // accepted state at *tret, unless the failure precedes the first step
  // (then yout is untouched)
  KROME_HD int solve(double tout, double* yout, double* tret) {
    using namespace detail_;
    int istate = SUCCESS;
    *tret = tn;
    {  // the first-call block (nst == 0)
      // cvInitialSetup
      const double tdist = fabs(tout - tn);
      const double tround = uround * smax_(fabs(tn), fabs(tout));
      if (tdist == 0.0 || tdist < 2.0 * tround) return TOO_CLOSE;
      if (ewt_set_(zn[0], ewt) != 0) return ILL_INPUT;
      const int retval = P.rhs(tn, zn[0], zn[1]);
      st.nfe++;
      if (retval < 0) return RHSFUNC_FAIL;
      if (retval > 0) return FIRST_RHSFUNC_ERR;
      h = hin;
      if (h != 0.0 && (tout - tn) * h < 0.0) return ILL_INPUT;
      if (h == 0.0) {
        const int hflag = hin_(tout);
        if (hflag != SUCCESS) return hflag;
      }
      const double rh = fabs(h) * hmax_inv;
      if (rh > 1.0) h /= rh;
      if (fabs(h) < hmin) h *= hmin / fabs(h);
      hscale = h;
      h0u = h;
      hprime = h;
      scale_<U>(N, h, zn[1], zn[1], S);
    }
    long nstloc = 0;
    for (;;) {
      next_h = h;
      next_q = q;
      if (st.nst > 0 && ewt_set_(zn[0], ewt) != 0) {
        istate = ILL_INPUT;
        *tret = tn;
        copy_(zn[0], yout);
        break;
      }
      if (mxstep > 0 && nstloc >= mxstep) {
        istate = TOO_MUCH_WORK;
        *tret = tn;
        copy_(zn[0], yout);
        break;
      }
      const double nrm = wrms_<U>(N, zn[0], ewt, S);
      tolsf = uround * nrm;
      if (tolsf > 1.0) {
        istate = TOO_MUCH_ACC;
        *tret = tn;
        copy_(zn[0], yout);
        tolsf *= 2.0;
        break;
      }
      tolsf = 1.0;
      const int kflag = step_();
      if (kflag != SUCCESS) {
        istate = kflag;
        *tret = tn;
        copy_(zn[0], yout);
        break;
      }
      nstloc++;
      if ((tn - tout) * h >= 0.0) {
        istate = SUCCESS;
        *tret = tout;
        get_dky_(tout, yout);
        next_q = qprime;
        next_h = hprime;
        break;
      }
    }
    return istate;
  }

  // CVodeGetCurrentStep / GetLastStep / GetCurrentOrder / GetLastOrder /
  // GetCurrentTime
  KROME_HD double current_step() const { return next_h; }
  KROME_HD double last_step() const { return hu; }
  KROME_HD int current_order() const { return next_q; }
  KROME_HD int last_order() const { return qu; }
  KROME_HD double current_time() const { return tn; }
  KROME_HD double tol_scale_factor() const { return tolsf; }
  // the error weights and the estimated local errors (CVodeGetErrWeights,
  // CVodeGetEstLocalErrors) into w and e (unit stride, N values)
  KROME_HD void get_weights_errors(double* w, double* e) const {
    for (int i = 0; i < N; i++) {
      if (w != nullptr) w[i] = ewt[ix(i)];
      if (e != nullptr) e[i] = acor[ix(i)];
    }
  }

 private:
  Prob& P;
  const int N;
  const Opts o;
  const long S;  // the vectors' stride (U = false)
  const double uround = DBL_EPSILON;
  const double hmin = 0.0, hmax_inv = 0.0;  // CVODE's defaults
  long mxstep;
  bool atolmin0;
  double* zn[L_MAX];
  double *ewt, *y, *acor, *tempv, *ftemp, *delta, *vtemp1, *vtemp2, *vtemp3;
  double *A, *savedJ;
  int* piv;
  // CVodeMem
  double tn = 0.0, h = 0.0, hin = 0.0, hprime = 0.0, hscale = 0.0, h0u = 0.0;
  double next_h = 0.0, hu = 0.0, eta = 1.0, etamax = 0.0, etaq = 0.0;
  double etaqm1 = 0.0, etaqp1 = 0.0, tolsf = 1.0, saved_tq5 = 0.0;
  double tau[L_MAX + 1] = {}, tq[6] = {}, l[L_MAX] = {};
  double rl1 = 0.0, gamma = 0.0, gammap = 0.0, gamrat = 0.0;
  double crate = 0.0, delp = 0.0, delnrm = 0.0, acnrm = 0.0;
  int q = 1, qprime = 1, next_q = 0, qu = 0, L = 2, qwait = 2;
  long nscon = 0, nstlp = 0, nstlj = 0;
  int convfail = 0;
  bool acnrmcur = false, jcur = false, nls_jcur = false;

  KROME_HD long ix(long i) const { return detail_::ix_<U>(i, S); }

  KROME_HD void copy_(const double* x, double* z) const {
    for (int i = 0; i < N; i++) z[ix(i)] = x[ix(i)];
  }

  // cvEwtSetSV (rtol scalar, atol vector) or KROME's ewt_cb_ (rtol_v)
  KROME_HD int ewt_set_(const double* ycur, double* weight) {
    // (the tolerances o.atol / o.rtol_v: unit stride, shared by the cells)
    if (o.rtol_v != nullptr) {
      for (int i = 0; i < N; i++) {
        const double d = fma(o.rtol_v[i], fabs(ycur[ix(i)]), o.atol[i]);
        if (!(d > 0.0)) return -1;
        weight[ix(i)] = 1.0 / d;
      }
      return 0;
    }
    // tempv = rtol |y| + atol (N_VLinearSum(rtol, |y|, 1, atol): VLin1,
    // (rtol |y|) + atol; rtol = 1, axpy: |y| + atol, the same value), then
    // N_VMin, N_VInv
    for (int i = 0; i < N; i++) weight[ix(i)] = fabs(ycur[ix(i)]);
    for (int i = 0; i < N; i++)
      weight[ix(i)] = (o.rtol * weight[ix(i)]) + o.atol[i];
    if (atolmin0) {
      double mn = weight[0];
      for (int i = 1; i < N; i++)
        if (weight[ix(i)] < mn) mn = weight[ix(i)];
      if (mn <= 0.0) return -1;
    }
    for (int i = 0; i < N; i++) weight[ix(i)] = 1.0 / weight[ix(i)];
    return 0;
  }

  // ---- initial step (cvHin, cvUpperBoundH0, cvYddNorm) ----------------

  KROME_HD int hin_(double tout) {
    using namespace detail_;
    const double tdiff = tout - tn;
    const int sign = (tdiff > 0.0) ? 1 : -1;
    const double tdist = fabs(tdiff);
    const double tround = uround * smax_(fabs(tn), fabs(tout));
    const double hlb = HLB_FACTOR * tround;
    const double hub = upper_bound_h0_(tdist);
    double hg = ssqrt_(hlb * hub);
    if (hub < hlb) {
      h = sign == -1 ? -hg : hg;
      return SUCCESS;
    }
    double hs = hg, hnew = hg, yddnrm = 0.0;
    for (int count1 = 1; count1 <= MAX_ITERS; count1++) {
      bool hgOK = false;
      for (int count2 = 1; count2 <= MAX_ITERS; count2++) {
        const double hgs = hg * sign;
        const int retval = ydd_norm_(hgs, &yddnrm);
        if (retval < 0) return RHSFUNC_FAIL;
        if (retval == SUCCESS) {
          hgOK = true;
          break;
        }
        hg *= 0.2;
      }
      if (!hgOK) {
        if (count1 <= 2) return REPTD_RHSFUNC_ERR;
        hnew = hs;
        break;
      }
      hs = hg;
      hnew = (yddnrm * hub * hub > 2.0) ? ssqrt_(2.0 / yddnrm)
                                         : ssqrt_(hg * hub);
      if (count1 == MAX_ITERS) break;
      const double hrat = hnew / hg;
      if ((hrat > 0.5) && (hrat < 2.0)) break;
      if ((count1 > 1) && (hrat > 2.0)) {
        hnew = hg;
        break;
      }
      hg = hnew;
    }
    double h0 = H_BIAS * hnew;
    if (h0 < hlb) h0 = hlb;
    if (h0 > hub) h0 = hub;
    if (sign == -1) h0 = -h0;
    h = h0;
    return SUCCESS;
  }

  KROME_HD double upper_bound_h0_(double tdist) {
    using namespace detail_;
    double* temp1 = tempv;
    double* temp2 = acor;
    for (int i = 0; i < N; i++) temp2[ix(i)] = fabs(zn[0][ix(i)]);
    ewt_set_(zn[0], temp1);
    for (int i = 0; i < N; i++) temp1[ix(i)] = 1.0 / temp1[ix(i)];
    lin_sum_<U>(N, HUB_FACTOR, temp2, 1.0, temp1, temp1, S);
    for (int i = 0; i < N; i++) temp2[ix(i)] = fabs(zn[1][ix(i)]);
    for (int i = 0; i < N; i++) temp1[ix(i)] = temp2[ix(i)] / temp1[ix(i)];
    const double hub_inv = max_norm_<U>(N, temp1, S);
    double hub = HUB_FACTOR * tdist;
    if (hub * hub_inv > 1.0) hub = 1.0 / hub_inv;
    return hub;
  }

  KROME_HD int ydd_norm_(double hg, double* yddnrm) {
    using namespace detail_;
    lin_sum_<U>(N, hg, zn[1], 1.0, zn[0], y, S);
    const int retval = P.rhs(tn + hg, y, tempv);
    st.nfe++;
    if (retval < 0) return RHSFUNC_FAIL;
    if (retval > 0) return RHSFUNC_RECVR;
    lin_sum_<U>(N, 1.0 / hg, tempv, -1.0 / hg, zn[1], tempv, S);
    *yddnrm = wrms_<U>(N, tempv, ewt, S);
    return SUCCESS;
  }

  // ---- one step (cvStep) ----------------------------------------------

  KROME_HD int step_() {
    using namespace detail_;
    int ncf = 0, nef = 0;
    double dsm = 0.0;
    if (st.nst > 0 && hprime != h) adjust_params_();
    const double saved_t = tn;
    int nflag = FIRST_CALL;
    for (;;) {
      predict_();
      set_();
      nflag = nls_(nflag);
      const int kflag = handle_nflag_(&nflag, saved_t, &ncf);
      if (kflag == PREDICT_AGAIN) continue;
      if (kflag != DO_ERROR_TEST) return kflag;
      const int eflag = do_error_test_(&nflag, saved_t, &nef, &dsm);
      if (eflag == TRY_AGAIN) continue;
      if (eflag != SUCCESS) return eflag;
      break;
    }
    complete_step_();
    prepare_next_step_(dsm);
    etamax = (st.nst <= SMALL_NST) ? ETA_MAX_ES : ETA_MAX_GS;
    scale_<U>(N, tq[2], acor, acor, S);
    return SUCCESS;
  }

  // cvAdjustParams, cvAdjustOrder, cvIncreaseBDF, cvDecreaseBDF, cvRescale
  KROME_HD void adjust_params_() {
    if (qprime != q) {
      adjust_order_(qprime - q);
      q = qprime;
      L = q + 1;
      qwait = L;
    }
    rescale_();
  }

  KROME_HD void adjust_order_(int deltaq) {
    if ((q == 2) && (deltaq != 1)) return;
    if (deltaq == 1) increase_bdf_();
    else if (deltaq == -1) decrease_bdf_();
  }

  KROME_HD void increase_bdf_() {
    for (int i = 0; i <= QMAX; i++) l[i] = 0.0;
    double alpha1, prod, xiold, alpha0;
    l[2] = alpha1 = prod = xiold = 1.0;
    alpha0 = -1.0;
    double hsum = hscale;
    if (q > 1) {
      for (int j = 1; j < q; j++) {
        hsum += tau[j + 1];
        const double xi = hsum / hscale;
        prod *= xi;
        alpha0 -= 1.0 / (j + 1);
        alpha1 += 1.0 / xi;
        for (int i = j + 2; i >= 2; i--) l[i] = l[i] * xiold + l[i - 1];
        xiold = xi;
      }
    }
    const double A1 = (-alpha0 - alpha1) / prod;
    detail_::scale_<U>(N, A1, zn[QMAX], zn[L], S);  // zn[indx_acor] = zn[qmax]
    if (q > 1)
      for (int j = 2; j <= q; j++) detail_::axpy_<U>(N, l[j], zn[L], zn[j], S);
  }

  KROME_HD void decrease_bdf_() {
    for (int i = 0; i <= QMAX; i++) l[i] = 0.0;
    l[2] = 1.0;
    double hsum = 0.0;
    for (int j = 1; j <= q - 2; j++) {
      hsum += tau[j];
      const double xi = hsum / hscale;
      for (int i = j + 2; i >= 2; i--) l[i] = l[i] * xi + l[i - 1];
    }
    if (q > 2)
      for (int j = 2; j < q; j++) detail_::axpy_<U>(N, -l[j], zn[q], zn[j], S);
  }

  KROME_HD void rescale_() {
    double c = eta;
    for (int j = 1; j <= q; j++) {
      detail_::scale_<U>(N, c, zn[j], zn[j], S);
      c = eta * c;
    }
    h = hscale * eta;
    next_h = h;
    hscale = h;
    nscon = 0;
  }

  // cvPredict (no tstop), cvRestore
  KROME_HD void predict_() {
    tn += h;
    for (int k = 1; k <= q; k++)
      for (int j = q; j >= k; j--) detail_::axpy_<U>(N, 1.0, zn[j], zn[j - 1], S);
  }

  KROME_HD void restore_(double saved_t) {
    tn = saved_t;
    for (int k = 1; k <= q; k++)
      for (int j = q; j >= k; j--) detail_::axpy_<U>(N, -1.0, zn[j], zn[j - 1], S);
  }

  // cvSet with cvSetBDF and cvSetTqBDF
  KROME_HD void set_() {
    double xi_inv, xistar_inv, alpha0, alpha0_hat;
    l[0] = l[1] = xi_inv = xistar_inv = 1.0;
    for (int i = 2; i <= q; i++) l[i] = 0.0;
    alpha0 = alpha0_hat = -1.0;
    double hsum = h;
    if (q > 1) {
      for (int j = 2; j < q; j++) {
        hsum += tau[j - 1];
        xi_inv = h / hsum;
        alpha0 -= 1.0 / j;
        for (int i = j; i >= 1; i--) l[i] += l[i - 1] * xi_inv;
      }
      alpha0 -= 1.0 / q;
      xistar_inv = -l[1] - alpha0;
      hsum += tau[q - 1];
      xi_inv = h / hsum;
      alpha0_hat = -l[1] - xi_inv;
      for (int i = q; i >= 1; i--) l[i] += l[i - 1] * xistar_inv;
    }
    set_tq_bdf_(hsum, alpha0, alpha0_hat, xi_inv, xistar_inv);
    rl1 = 1.0 / l[1];
    gamma = h * rl1;
    if (st.nst == 0) gammap = gamma;
    gamrat = (st.nst > 0) ? gamma / gammap : 1.0;
  }

  KROME_HD void set_tq_bdf_(double hsum, double alpha0, double alpha0_hat,
                            double xi_inv, double xistar_inv) {
    const double A1 = 1.0 - alpha0_hat + alpha0;
    const double A2 = 1.0 + q * A1;
    tq[2] = fabs(A1 / (alpha0 * A2));
    tq[5] = fabs(A2 * xistar_inv / (l[q] * xi_inv));
    if (qwait == 1) {
      if (q > 1) {
        const double C = xistar_inv / l[q];
        const double A3 = alpha0 + 1.0 / q;
        const double A4 = alpha0_hat + xi_inv;
        const double Cpinv = (1.0 - A4 + A3) / A3;
        tq[1] = fabs(C * Cpinv);
      } else {
        tq[1] = 1.0;
      }
      hsum += tau[q];
      xi_inv = h / hsum;
      const double A5 = alpha0 - (1.0 / (q + 1));
      const double A6 = alpha0_hat - xi_inv;
      const double Cppinv = (1.0 - A6 + A5) / A2;
      tq[3] = fabs(Cppinv / (xi_inv * (q + 2) * A5));
    }
    tq[4] = detail_::CORTES / tq[2];
  }

  // ---- the nonlinear solve (cvNls + SUNNonlinSolSolve_Newton) ---------

  KROME_HD int nls_(int nflag) {
    using namespace detail_;
    convfail = ((nflag == FIRST_CALL) || (nflag == PREV_ERR_FAIL))
                   ? NO_FAILURES : FAIL_OTHER;
    bool callSetup = (nflag == PREV_CONV_FAIL) || (nflag == PREV_ERR_FAIL)
                     || (st.nst == 0) || (st.nst >= nstlp + MSBP)
                     || (fabs(gamrat - 1.0) > DGMAX_LSETUP);
    const_<U>(N, 0.0, acor, S);
    const double tol = tq[4];
    bool jbad = false;
    long niters = 0, nconvfails = 0;
    int retval;
    for (;;) {
      int curiter = 0;
      retval = residual_();
      if (retval != SUCCESS) break;
      if (callSetup) {
        retval = nls_lsetup_(jbad, &nls_jcur);
        if (retval != SUCCESS) break;
      }
      for (;;) {
        niters++;
        const bool inv = o.dq_jac && o.inv_row >= 0;
        const double cg = inv ? inv_dot_(delta) : 0.0;  // c . G
        scale_<U>(N, -1.0, delta, delta, S);
        retval = lsolve_(delta);
        if (retval != SUCCESS) break;
        if (inv) inv_set_(delta, -cg);
        lin_sum_<U>(N, 1.0, acor, 1.0, delta, acor, S);
        retval = conv_test_(curiter, tol);
        curiter++;
        if (retval == SUCCESS) {
          nls_jcur = false;
          break;
        }
        if (retval != NLS_CONTINUE) break;
        if (curiter >= NLS_MAXCOR) {
          retval = NLS_CONV_RECVR;
          break;
        }
        retval = residual_();
        if (retval != SUCCESS) break;
      }
      if (retval == SUCCESS) break;
      if ((retval > 0) && !nls_jcur) {
        nconvfails++;
        callSetup = true;
        jbad = true;
        const_<U>(N, 0.0, acor, S);
        continue;
      }
      break;
    }
    if (retval != SUCCESS) nconvfails++;
    st.nni += niters;
    st.nnf += nconvfails;
    if (retval != SUCCESS) return retval;
    lin_sum_<U>(N, 1.0, zn[0], 1.0, acor, y, S);
    if (!acnrmcur) acnrm = wrms_<U>(N, acor, ewt, S);
    jcur = false;
    return SUCCESS;
  }

  // cvNlsResidual: y = zn[0] + acor, delta = rl1 zn[1] + acor - gamma f(y)
  KROME_HD int residual_() {
    using namespace detail_;
    lin_sum_<U>(N, 1.0, zn[0], 1.0, acor, y, S);
    const int retval = P.rhs(tn, y, ftemp);
    st.nfe++;
    if (retval < 0) return RHSFUNC_FAIL;
    if (retval > 0) return RHSFUNC_RECVR;
    lin_sum_<U>(N, rl1, zn[1], 1.0, acor, delta, S);
    lin_sum_<U>(N, -gamma, ftemp, 1.0, delta, delta, S);
    return SUCCESS;
  }

  // cvNlsConvTest
  KROME_HD int conv_test_(int m, double tol) {
    using namespace detail_;
    delnrm = wrms_<U>(N, delta, ewt, S);
    if (m > 0) crate = smax_(CRDOWN * crate, delnrm / delp);
    const double dcon = delnrm * smin_(1.0, crate) / tol;
    if (dcon <= 1.0) {
      acnrm = (m == 0) ? delnrm : wrms_<U>(N, acor, ewt, S);
      acnrmcur = true;
      return SUCCESS;
    }
    if ((m >= 1) && (delnrm > RDIV * delp)) return NLS_CONV_RECVR;
    delp = delnrm;
    return NLS_CONTINUE;
  }

  // cvNlsLSetup
  KROME_HD int nls_lsetup_(bool jbad, bool* jcur_out) {
    using namespace detail_;
    if (jbad) convfail = FAIL_BAD_J;
    const int retval = ls_setup_();
    st.nsetups++;
    *jcur_out = jcur;
    gamrat = 1.0;
    gammap = gamma;
    crate = 1.0;
    delnrm = 0.0;
    nstlp = st.nst;
    if (retval < 0) return LSETUP_FAIL;
    if (retval > 0) return NLS_CONV_RECVR;
    return SUCCESS;
  }

  // cvLsSetup + cvLsLinSys (A = I - gamma J, J reused or re-evaluated at
  // the predicted y and f(y)) + SUNLinSolSetup_Dense
  KROME_HD int ls_setup_() {
    using namespace detail_;
    const long NN = (long)N * N;
    const double dgamma = fabs((gamma / gammap) - 1.0);
    const bool jbad_ls = (st.nst == 0) || (st.nst >= nstlj + MSBJ)
                         || ((convfail == FAIL_BAD_J) && (dgamma < DGMAX_JBAD))
                         || (convfail == FAIL_OTHER);
    // CVODE evaluates J into A, copies it to savedJ and scales A in place
    // (zero, fill, copy, scale: four passes over the matrix). Here J goes
    // straight into savedJ and A = I - gamma savedJ is one pass, with the
    // same operations (A = savedJ * c, then + 1 on the diagonal), so the
    // result is bitwise CVODE's. A failed evaluation leaves savedJ partly
    // written where CVODE kept the previous J, but that J is never reused:
    // the retry after a setup failure has convfail FAIL_OTHER, so its setup
    // evaluates J again.
    int retval = 0;
    if (!jbad_ls) {
      jcur = false;
    } else {
      jcur = true;
      if (!o.dq_jac)  // the DQ Jacobian writes every element
        for (long k = 0; k < NN; k++) savedJ[ix(k)] = 0.0;
      const int jr = o.dq_jac
          ? dq_jac_(tn, y, ftemp, savedJ, vtemp1)
          : P.jac(tn, y, ftemp, savedJ, vtemp1, vtemp2, vtemp3);
      if (jr < 0) retval = -1;
      else if (jr > 0) retval = 1;
    }
    // SUNMatScaleAddI_Dense(-gamma, A), A = savedJ (the sparse LU gathers
    // its entries from savedJ instead: sp_setup_)
    if (retval == 0 && !o.sparse) {
      const double c = -gamma;
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          const long k = i + (long)j * N;
          A[ix(k)] = savedJ[ix(k)] * c;
        }
        A[ix(j + (long)j * N)] += 1.0;
      }
    }
    if (jcur) {
      st.nje++;
      nstlj = st.nst;
    }
    if (retval != 0) return retval;
    if (o.sparse) return sp_setup_(-gamma);
#if defined(KROME_BDF_NEWTON_HOOK) && !defined(__CUDA_ARCH__) \
    && !defined(__HIP_DEVICE_COMPILE__)
    if constexpr (U) ::KROME_BDF_NEWTON_HOOK(A, N);
#endif
#ifdef KROME_BDF_LAPACK_
    if constexpr (U)
      if (o.lapack) {  // SUNLinSolSetup_LapackDense
        const int n = N;
        int ier = 0;
        dgetrf_(&n, &n, A, &n, piv, &ier);
        if (ier > 0) return LS_LUFACT_FAIL;
        return ier < 0 ? -1 : 0;
      }
#endif
    return getrf_<U>(A, N, piv, S) > 0 ? LS_LUFACT_FAIL : 0;
  }

  // cvLsSolve (dense): x = A^-1 b, scaled by 2 / (1 + gamrat) when gamma
  // changed since the setup
  KROME_HD int lsolve_(double* b) {
    if (o.sparse) {
      sp_solve_(b);
      if (gamrat != 1.0) detail_::scale_<U>(N, 2.0 / (1.0 + gamrat), b, b, S);
      return SUCCESS;
    }
#ifdef KROME_BDF_LAPACK_
    if constexpr (U)
      if (o.lapack) {  // SUNLinSolSolve_LapackDense
        const int n = N, one = 1;
        int ier = 0;
        dgetrs_("N", &n, &one, A, &n, piv, b, &n, &ier);
        if (gamrat != 1.0) detail_::scale_<U>(N, 2.0 / (1.0 + gamrat), b, b, S);
        return ier < 0 ? LSOLVE_FAIL : SUCCESS;
      }
#endif
    detail_::getrs_<U>(A, N, piv, b, S);
    if (gamrat != 1.0) detail_::scale_<U>(N, 2.0 / (1.0 + gamrat), b, b, S);
    return SUCCESS;
  }

  // ---- the static-pivot sparse LU (Opts::sparse) ----------------------
  // A's storage: the factors' values (nnz, the tables' positions), the row
  // accumulator w (ns), the border B' = L^-1 B (ns x m, column c at c ns),
  // C' = C U^-1 (m x ns, row r at r ns) and the Schur complement
  // Z = D - C' B' (m x m, column-major; its pivots in piv)

  // A = c savedJ + I (c = -gamma) gathered and factorized: 0, or
  // LS_LUFACT_FAIL at a zero or non-finite pivot
  KROME_HD int sp_setup_(double c) {
    using namespace detail_;
    const Sparse& sp = o.sp;
    const int ns = sp.ns, m = N - ns;
    double* val = A;
    double* w = A + ix(sp.nnz);
    double* Bt = w + ix(ns);
    double* Ct = Bt + ix((long)ns * m);
    double* Z = Ct + ix((long)m * ns);
    for (int k = 0; k < ns; k++) {
      const long r = sp.perm[k];
      for (int p = sp.rp[k]; p < sp.rp[k + 1]; p++)
        val[ix(p)] = savedJ[ix(r + (long)sp.perm[sp.ci[p]] * N)] * c;
      val[ix(sp.di[k])] += 1.0;
    }
    for (int i = 0; i < ns; i++) {  // row i of L and U (the IKJ form)
      const int r0 = sp.rp[i], r1 = sp.rp[i + 1], d = sp.di[i];
      for (int p = r0; p < r1; p++) w[ix(sp.ci[p])] = val[ix(p)];
      for (int p = r0; p < d; p++) {
        const int k = sp.ci[p];
        const double lik = w[ix(k)] / val[ix(sp.di[k])];
        w[ix(k)] = lik;
        for (int q = sp.di[k] + 1; q < sp.rp[k + 1]; q++)
          w[ix(sp.ci[q])] -= lik * val[ix(q)];
      }
      const double piv = w[ix(i)];
      if (piv == 0.0 || !(fabs(piv) <= DBL_MAX)) return LS_LUFACT_FAIL;
      for (int p = r0; p < r1; p++) val[ix(p)] = w[ix(sp.ci[p])];
    }
    if (m == 0) return 0;
    for (int b = 0; b < m; b++) {  // B' = L^-1 B
      double* x = Bt + ix((long)b * ns);
      for (int k = 0; k < ns; k++)
        x[ix(k)] = savedJ[ix(sp.perm[k] + (long)(ns + b) * N)] * c;
      sp_lower_(x);
    }
    for (int r = 0; r < m; r++) {  // C' = C U^-1 (x U = row r of C)
      double* x = Ct + ix((long)r * ns);
      for (int k = 0; k < ns; k++)
        x[ix(k)] = savedJ[ix(ns + r + (long)sp.perm[k] * N)] * c;
      for (int k = 0; k < ns; k++) {
        const double xk = x[ix(k)] / val[ix(sp.di[k])];
        x[ix(k)] = xk;
        if (xk != 0.0)
          for (int q = sp.di[k] + 1; q < sp.rp[k + 1]; q++)
            x[ix(sp.ci[q])] -= xk * val[ix(q)];
      }
    }
    for (int b = 0; b < m; b++)  // Z = D - C' B'
      for (int r = 0; r < m; r++) {
        double z = savedJ[ix(ns + r + (long)(ns + b) * N)] * c;
        if (r == b) z += 1.0;
        for (int k = 0; k < ns; k++)
          z -= Ct[ix((long)r * ns + k)] * Bt[ix((long)b * ns + k)];
        Z[ix(r + (long)b * m)] = z;
      }
    return getrf_<U>(Z, m, piv, S) > 0 ? LS_LUFACT_FAIL : 0;
  }

  // x = L^-1 x (unit lower triangle; x in pivot order)
  KROME_HD void sp_lower_(double* x) const {
    const Sparse& sp = o.sp;
    const double* val = A;
    for (int i = 0; i < sp.ns; i++) {
      double s = x[ix(i)];
      for (int p = sp.rp[i]; p < sp.di[i]; p++) s -= val[ix(p)] * x[ix(sp.ci[p])];
      x[ix(i)] = s;
    }
  }

  // b = A^-1 b with the factors of sp_setup_
  KROME_HD void sp_solve_(double* b) {
    const Sparse& sp = o.sp;
    const int ns = sp.ns, m = N - ns;
    const double* val = A;
    double* x = A + ix(sp.nnz);
    const double* Bt = x + ix(ns);
    const double* Ct = Bt + ix((long)ns * m);
    double* Z = const_cast<double*>(Ct) + ix((long)m * ns);
    for (int k = 0; k < ns; k++) x[ix(k)] = b[ix(sp.perm[k])];
    sp_lower_(x);
    if (m > 0) {
      for (int r = 0; r < m; r++) {
        double s = b[ix(ns + r)];
        for (int k = 0; k < ns; k++) s -= Ct[ix((long)r * ns + k)] * x[ix(k)];
        b[ix(ns + r)] = s;
      }
      detail_::getrs_<U>(Z, m, piv, b + ix(ns), S);
      for (int bb = 0; bb < m; bb++) {
        const double v = b[ix(ns + bb)];
        if (v != 0.0)
          for (int k = 0; k < ns; k++) x[ix(k)] -= Bt[ix((long)bb * ns + k)] * v;
      }
    }
    for (int i = ns - 1; i >= 0; i--) {
      double s = x[ix(i)];
      for (int p = sp.di[i] + 1; p < sp.rp[i + 1]; p++)
        s -= val[ix(p)] * x[ix(sp.ci[p])];
      x[ix(i)] = s / val[ix(sp.di[i])];
    }
    for (int k = 0; k < ns; k++) b[ix(sp.perm[k])] = x[ix(k)];
  }

  // Opts::inv_row: c . x over the invariant's components
  KROME_HD double inv_dot_(const double* x) const {
    double s = 0.0;
    for (int i = 0; i < o.inv_n; i++)
      if (o.inv_c[i] != 0) s += o.inv_c[i] * x[ix(i)];
    return s;
  }
  // component inv_row of the correction x set so that c . x = target,
  // where the iterate y = zn[0] + acor carries enough of the invariant
  // there (Opts::inv_row) and stays non-negative
  KROME_HD void inv_set_(double* x, double target) const {
    const int r = o.inv_row;
    double s = 0.0, sa = 0.0;
    for (int i = 0; i < o.inv_n; i++)
      if (i != r && o.inv_c[i] != 0) {
        s += o.inv_c[i] * x[ix(i)];
        sa += fabs(o.inv_c[i] * (zn[0][ix(i)] + acor[ix(i)]));
      }
    const double yr = zn[0][ix(r)] + acor[ix(r)];
    const double xr = (target - s) / o.inv_c[r];
    if (fabs(o.inv_c[r] * yr) >= 1e-6 * sa && yr + xr >= 0.0) x[ix(r)] = xr;
  }

  // cvLsDenseDQJac (no constraints): column j = (f(y + inc e_j) - fy) / inc
  KROME_HD int dq_jac_(double t, double* yv, const double* fy, double* Jac,
                       double* tmp1) {
    using namespace detail_;
    double* ftmp = tmp1;
    const double srur = ssqrt_(uround);
    const double fnorm = wrms_<U>(N, fy, ewt, S);
    const double minInc = (fnorm != 0.0)
        ? (MIN_INC_MULT * fabs(h) * uround * (double)N * fnorm) : 1.0;
    int retval = 0;
    for (int j = 0; j < N; j++) {
      const double yjsaved = yv[ix(j)];
      const double inc = smax_(srur * fabs(yjsaved), minInc / ewt[ix(j)]);
      yv[ix(j)] += inc;
      retval = P.rhs(t, yv, ftmp);
      st.nfeDQ++;
      if (retval != 0) break;
      yv[ix(j)] = yjsaved;
      const double inc_inv = 1.0 / inc;
      lin_sum_<U>(N, inc_inv, ftmp, -inc_inv, fy, Jac + ix((long)j * N), S);
    }
    return retval;
  }

  // ---- after the nonlinear solve (cvHandleNFlag, cvDoErrorTest) -------

  KROME_HD int handle_nflag_(int* nflagPtr, double saved_t, int* ncfPtr) {
    using namespace detail_;
    const int nflag = *nflagPtr;
    if (nflag == SUCCESS) return DO_ERROR_TEST;
    st.ncfn++;
    restore_(saved_t);
    if (nflag < 0) {
      if (nflag == LSETUP_FAIL) return LSETUP_FAIL;
      if (nflag == LSOLVE_FAIL) return LSOLVE_FAIL;
      if (nflag == RHSFUNC_FAIL) return RHSFUNC_FAIL;
      return NLS_FAIL;
    }
    (*ncfPtr)++;
    etamax = 1.0;
    if ((fabs(h) <= hmin * ONEPSM) || (*ncfPtr == o.max_conv_fails)) {
      if (nflag == NLS_CONV_RECVR) return CONV_FAILURE;
      if (nflag == RHSFUNC_RECVR) return REPTD_RHSFUNC_ERR;
    }
    eta = smax_(ETA_CF, hmin / fabs(h));
    *nflagPtr = PREV_CONV_FAIL;
    rescale_();
    return PREDICT_AGAIN;
  }

  KROME_HD int do_error_test_(int* nflagPtr, double saved_t, int* nefPtr,
                              double* dsmPtr) {
    using namespace detail_;
    const double dsm = acnrm * tq[2];
    *dsmPtr = dsm;
    if (dsm <= 1.0) return SUCCESS;
    (*nefPtr)++;
    st.netf++;
    *nflagPtr = PREV_ERR_FAIL;
    restore_(saved_t);
    if ((fabs(h) <= hmin * ONEPSM) || (*nefPtr == o.max_err_fails))
      return ERR_FAILURE;
    etamax = 1.0;
    if (*nefPtr <= MXNEF1) {
      eta = 1.0 / (pow(BIAS2 * dsm, 1.0 / L) + ADDON);
      eta = smax_(ETA_MIN_EF, smax_(eta, hmin / fabs(h)));
      if (*nefPtr >= SMALL_NEF) eta = smin_(eta, ETA_MAX_EF);
      rescale_();
      return TRY_AGAIN;
    }
    if (q > 1) {
      eta = smax_(ETA_MIN_EF, hmin / fabs(h));
      adjust_order_(-1);
      L = q;
      q--;
      qwait = L;
      rescale_();
      return TRY_AGAIN;
    }
    // at order 1: restart from y' = f(tn, y)
    eta = smax_(ETA_MIN_EF, hmin / fabs(h));
    h *= eta;
    next_h = h;
    hscale = h;
    qwait = LONG_WAIT;
    nscon = 0;
    const int retval = P.rhs(tn, zn[0], tempv);
    st.nfe++;
    if (retval < 0) return RHSFUNC_FAIL;
    if (retval > 0) return UNREC_RHSFUNC_ERR;
    scale_<U>(N, h, tempv, zn[1], S);
    return TRY_AGAIN;
  }

  // ---- after a successful step ----------------------------------------

  // cvCompleteStep
  KROME_HD void complete_step_() {
    st.nst++;
    nscon++;
    hu = h;
    qu = q;
    for (int i = q; i >= 2; i--) tau[i] = tau[i - 1];
    if ((q == 1) && (st.nst > 1)) tau[2] = tau[1];
    tau[1] = h;
    for (int j = 0; j <= q; j++) detail_::axpy_<U>(N, l[j], acor, zn[j], S);
    qwait--;
    if ((qwait == 1) && (q != QMAX)) {
      copy_(acor, zn[QMAX]);
      saved_tq5 = tq[5];
    }
  }

  // cvPrepareNextStep
  KROME_HD void prepare_next_step_(double dsm) {
    using namespace detail_;
    if (etamax == 1.0) {
      qwait = qwait > 2 ? qwait : 2;
      qprime = q;
      hprime = h;
      eta = 1.0;
      return;
    }
    etaq = 1.0 / (pow(BIAS2 * dsm, 1.0 / L) + ADDON);
    if (qwait != 0) {
      eta = etaq;
      qprime = q;
      set_eta_();
    } else {
      qwait = 2;
      etaqm1 = compute_etaqm1_();
      etaqp1 = compute_etaqp1_();
      choose_eta_();
      set_eta_();
    }
  }

  // cvSetEta
  KROME_HD void set_eta_() {
    using namespace detail_;
    if ((eta > ETA_MIN_FX) && (eta < ETA_MAX_FX)) {
      eta = 1.0;
      hprime = h;
      return;
    }
    if (eta >= ETA_MAX_FX) {
      eta = smin_(eta, etamax);
      eta /= smax_(1.0, fabs(h) * hmax_inv * eta);
    } else {
      eta = smax_(eta, ETA_MIN);
      eta = smax_(eta, hmin / fabs(h));
    }
    hprime = h * eta;
    if (qprime < q) nscon = 0;
  }

  // cvComputeEtaqm1, cvComputeEtaqp1, cvChooseEta
  KROME_HD double compute_etaqm1_() {
    using namespace detail_;
    etaqm1 = 0.0;
    if (q > 1) {
      const double ddn = wrms_<U>(N, zn[q], ewt, S) * tq[1];
      etaqm1 = 1.0 / (pow(BIAS1 * ddn, 1.0 / q) + ADDON);
    }
    return etaqm1;
  }

  KROME_HD double compute_etaqp1_() {
    using namespace detail_;
    etaqp1 = 0.0;
    if (q != QMAX) {
      if (saved_tq5 == 0.0) return etaqp1;
      const double cquot = (tq[5] / saved_tq5) * spowi_(h / tau[2], L);
      lin_sum_<U>(N, -cquot, zn[QMAX], 1.0, acor, tempv, S);
      const double dup = wrms_<U>(N, tempv, ewt, S) * tq[3];
      etaqp1 = 1.0 / (pow(BIAS3 * dup, 1.0 / (L + 1)) + ADDON);
    }
    return etaqp1;
  }

  KROME_HD void choose_eta_() {
    using namespace detail_;
    const double etam = smax_(etaqm1, smax_(etaq, etaqp1));
    if ((etam > ETA_MIN_FX) && (etam < ETA_MAX_FX)) {
      eta = 1.0;
      qprime = q;
      return;
    }
    if (etam == etaq) {
      eta = etaq;
      qprime = q;
    } else if (etam == etaqm1) {
      eta = etaqm1;
      qprime = q - 1;
    } else {
      eta = etaqp1;
      qprime = q + 1;
      copy_(acor, zn[QMAX]);  // Delta_n for the order increase
    }
  }

  // CVodeGetDky(t, 0): the interpolating polynomial at t (N_VLinearCombination's
  // generic form: a scale and axpys, highest order first)
  KROME_HD void get_dky_(double t, double* dky) {
    const double s = (t - tn) / h;
    for (int j = q; j >= 0; j--) {
      double c = 1.0;
      for (int i = 0; i < j; i++) c *= s;
      if (j == q) detail_::scale_<U>(N, c, zn[j], dky, S);
      else detail_::axpy_<U>(N, c, zn[j], dky, S);
    }
  }
};

// KROME's use of CVODE (solve_cell.hpp integrate_): one CVode(dt) call from
// t = 0 with y; error-test, convergence and repeated-rhs failures restart
// at the reached state from a tenth of the step CVODE had reduced to, up
// to max_restarts times. y: in, the initial state; out, y(dt) on success
// or the last accepted state on failure. Returns 0 or a CVODE flag.
//
// U = false: y, work and iwork are a cell's slices of interleaved arrays
// (element i at stride S; see Cell).
//
// diag (optional): the integrator's final state (Diag).
template <class Prob, bool U = true>
KROME_HD inline int integrate(Prob& prob, int N, double* y, double dt,
                              const Opts& o, int max_restarts, double* work,
                              int* iwork, Stats* stats, int* nrestarts,
                              long S = 1, Diag* diag = nullptr) {
  Cell<Prob, U> cv(prob, N, work, iwork, o, S);
  cv.reinit(0.0, y, o.h0);
  double t = 0.0;
  int restarts = 0;
  auto finish = [&](int flag) {
    if (nrestarts != nullptr) *nrestarts = restarts;
    if (diag != nullptr) {
      diag->tcur = cv.current_time();
      diag->hlast = cv.last_step();
      diag->hcur = cv.current_step();
      diag->tolsf = cv.tol_scale_factor();
      diag->qlast = cv.last_order();
      diag->qcur = cv.current_order();
      diag->last = cv.st;
      cv.get_weights_errors(diag->ewt, diag->ele);
    }
    return flag;
  };
  for (;;) {
    const int flag = cv.solve(dt, y, &t);
    if (stats != nullptr) stats->add(cv.st);
    if (flag >= 0) return finish(0);
    const bool restartable = flag == ERR_FAILURE || flag == CONV_FAILURE
                             || flag == REPTD_RHSFUNC_ERR;
    if (!restartable || restarts >= max_restarts) return finish(flag);
    ++restarts;
    const double hcur = cv.current_step();
    cv.reinit(t, y, hcur > 0.0 ? 0.1 * hcur : 0.0);
  }
}

}  // namespace bdf
}  // namespace krome
