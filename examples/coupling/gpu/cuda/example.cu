/*****************************************************************
 KROME-X in a GPU-resident host code (CUDA C++): the hydro state stays
 on the GPU for the whole run, one array per variable (u[v * ncell + c]:
 density, internal energy density, the passive scalars rho x of each
 species). Each step the chemistry
  1. gathers the cells into KROME-X's arrays (a kernel): mass fractions
     x[i * ncell + c] (species-major: KROME-X's default layout), rho, and
     Tgas from the internal energy (cgs: a real code converts its units
     here);
  2. restores each cell's element composition after the advection
     (consistent advection, on the device);
  3. integrates every cell over dt on the device (KROME-X, one call);
  4. scatters the results back (a kernel): rho x, and the internal energy
     from the new Tgas and composition.
 The cells: n_H from 1 to 1e6 cm-3 at 100 K, the free-fall time of each
 density; one is the quick start's (n_H = 1e4 cm-3), whose history goes
 to abundances.dat as in the CPU examples (README.md), from 1 yr to
 1e6 yr in steps 10 % longer each time. The same as the CUDA Fortran
 example, with a flat cell array instead of octs.
*****************************************************************/
#include <cuda_runtime.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.h"       /* the C interface: krome_solve_cells_x_device */
#include "krome_net.h"   /* KROME_NSPEC, KROME_IDX_*, KROME_PAR_*, ... */

#define CK(x)                                                         \
  do {                                                                \
    cudaError_t e_ = (x);                                             \
    if (e_ != cudaSuccess) {                                          \
      printf("CUDA error: %s (%s:%d)\n", cudaGetErrorString(e_),      \
             __FILE__, __LINE__);                                     \
      return 1;                                                       \
    }                                                                 \
  } while (0)

static const int NS = KROME_NSPEC;
enum { IRHO = 0, IEINT = 1, IPASS = 2, NVAR = IPASS + KROME_NSPEC };
static const double KB = KROME_CONST_boltzmann_erg;  /* erg / K */

/* u -> x, rho, Tgas = (gamma - 1) rho e / (kB n_tot) */
__global__ void gather(int ncell, const double* u, const double* m,
                       double gam, double* x, double* rho, double* T) {
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= ncell) return;
  const double r = u[IRHO * ncell + c];
  double ntot = 0.0;
  for (int i = 0; i < NS; ++i) {
    const double rx = u[(IPASS + i) * ncell + c];
    x[i * ncell + c] = rx / r;
    ntot += rx / m[i];
  }
  rho[c] = r;
  T[c] = (gam - 1.0) * u[IEINT * ncell + c] / (KB * ntot);
}

/* back into u; a failed cell (flag < 0) keeps its state */
__global__ void scatter(int ncell, double* u, const double* m, double gam,
                        const double* x, const double* rho,
                        const double* T, const int* flags) {
  const int c = blockIdx.x * blockDim.x + threadIdx.x;
  if (c >= ncell || flags[c] < 0) return;
  double ntot = 0.0;
  for (int i = 0; i < NS; ++i) {
    const double rx = rho[c] * x[i * ncell + c];
    u[(IPASS + i) * ncell + c] = rx;
    ntot += rx / m[i];
  }
  u[IEINT * ncell + c] = KB * ntot * T[c] / (gam - 1.0);
}

int main(void) {
  const double spy = 3.15e7;            /* seconds per year */
  const int ncell = 32768, iref = ncell / 2;   /* iref: the quick start */
  /* one GPU per MPI rank: krome_device_initialize_id(rank % ngpu) */
  krome_device_initialize_id(0);
  krome_opts* h = krome_opts_new();
  krome_opts_set(h, "evolve_T", 1);     /* Tgas evolved */

  /* the cells (on the host only to set them up) */
  std::vector<double> m(NS), u((size_t)NVAR * ncell),
      par((size_t)KROME_NPAR * ncell);
  for (int i = 0; i < NS; ++i) m[i] = krome_species_mass(i);
  for (int c = 0; c < ncell; ++c)
    krome_params_default(&par[(size_t)c * KROME_NPAR]);
  const double gam = par[KROME_PAR_gamma_ad];   /* the network's gamma */
  for (int c = 0; c < ncell; ++c) {
    const double nH = c == iref ? 1e4 : pow(10.0, 6.0 * c / (ncell - 1));
    double n[KROME_NSPEC] = {0.0}, rho = 0.0, ntot = 0.0;
    n[KROME_IDX_H] = nH;
    n[KROME_IDX_Hj] = n[KROME_IDX_E] = 1e-4 * nH;
    n[KROME_IDX_He] = 0.0775 * nH;
    for (int i = 0; i < NS; ++i) {
      rho += m[i] * n[i];
      ntot += n[i];
      u[(size_t)(IPASS + i) * ncell + c] = m[i] * n[i];
    }
    u[(size_t)IRHO * ncell + c] = rho;
    u[(size_t)IEINT * ncell + c] = KB * ntot * 100.0 / (gam - 1.0);
    par[(size_t)c * KROME_NPAR + KROME_PAR_tff] = 1e14 * sqrt(1e4 / nH);
  }
  /* the reference element composition per gram (the same for every
     cell here), for the consistent advection */
  double xref[KROME_NSPEC], ref[KROME_NELEM];
  for (int i = 0; i < NS; ++i)
    xref[i] = u[(size_t)(IPASS + i) * ncell] / u[IRHO * ncell];
  krome_get_element_totals_x(1, xref, ref);

  /* the host code's device data, and KROME-X's arguments */
  double *u_d, *m_d, *x_d, *rho_d, *T_d, *dt_d, *par_d, *ref_d;
  int* flags_d;
  CK(cudaMalloc(&u_d, u.size() * sizeof(double)));
  CK(cudaMalloc(&m_d, NS * sizeof(double)));
  CK(cudaMalloc(&x_d, (size_t)NS * ncell * sizeof(double)));
  CK(cudaMalloc(&rho_d, ncell * sizeof(double)));
  CK(cudaMalloc(&T_d, ncell * sizeof(double)));
  CK(cudaMalloc(&dt_d, ncell * sizeof(double)));
  CK(cudaMalloc(&par_d, par.size() * sizeof(double)));
  CK(cudaMalloc(&ref_d, KROME_NELEM * sizeof(double)));
  CK(cudaMalloc(&flags_d, ncell * sizeof(int)));
  CK(cudaMemcpy(u_d, u.data(), u.size() * sizeof(double),
                cudaMemcpyHostToDevice));
  CK(cudaMemcpy(m_d, m.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
  CK(cudaMemcpy(par_d, par.data(), par.size() * sizeof(double),
                cudaMemcpyHostToDevice));
  CK(cudaMemcpy(ref_d, ref, sizeof(ref), cudaMemcpyHostToDevice));

  FILE* out = fopen("abundances.dat", "w");
  fprintf(out, "#t_yr Tgas %s\n", krome_names_header());
  printf("%5s%11s%11s\n", "step", "t(yr)", "Tgas(K)");
  const int nt = 128, nb = (ncell + nt - 1) / nt;
  std::vector<double> dtv(ncell);
  double t = 0.0, dt = spy;
  long nfailed = 0;
  int step = 0;
  const auto t0 = std::chrono::steady_clock::now();
  while (t < 1e6 * spy) {
    /* ... the hydro step (advection of rho, rho e, rho x) would be here */
    gather<<<nb, nt>>>(ncell, u_d, m_d, gam, x_d, rho_d, T_d);
    krome_conserve_cells_x_device(ncell, x_d, KROME_LAYOUT_SPECIES_MAJOR,
                                  ncell, ref_d, 0, par_d, KROME_NPAR, NULL);
    for (int c = 0; c < ncell; ++c) dtv[c] = dt;
    CK(cudaMemcpy(dt_d, dtv.data(), ncell * sizeof(double),
                  cudaMemcpyHostToDevice));
    krome_device_summary s;
    krome_solve_cells_x_device(ncell, x_d, rho_d, KROME_LAYOUT_SPECIES_MAJOR,
                               ncell, T_d, dt_d, par_d, KROME_NPAR, NULL,
                               flags_d, &s, h);
    /* s: cells solved on the device, recovered on the host, failed. A
       failed cell (flags_d < 0) is left as it was: retry it (a smaller dt,
       a larger max_steps), or accept it unchanged for this step */
    nfailed += s.failed;
    scatter<<<nb, nt>>>(ncell, u_d, m_d, gam, x_d, rho_d, T_d, flags_d);

    t += dt;
    dt *= 1.1;
    ++step;
    /* the quick start's cell, to the host for the output */
    double uc[NVAR], n[KROME_NSPEC], ntot = 0.0;
    for (int v = 0; v < NVAR; ++v)
      CK(cudaMemcpy(&uc[v], u_d + (size_t)v * ncell + iref, sizeof(double),
                    cudaMemcpyDeviceToHost));
    for (int i = 0; i < NS; ++i) {
      n[i] = uc[IPASS + i] / m[i];
      ntot += n[i];
    }
    const double Tref = (gam - 1.0) * uc[IEINT] / (KB * ntot);
    fprintf(out, "%17.8E %17.8E", t / spy, Tref);
    for (int i = 0; i < NS; ++i) fprintf(out, " %17.8E", n[i] / 1e4);
    fprintf(out, "\n");
    if (step % 20 == 0) printf("%5d%11.3E%11.3E\n", step, t / spy, Tref);
  }
  CK(cudaDeviceSynchronize());
  const double wall = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  fclose(out);
  printf("wrote abundances.dat; %d steps of %d cells in %.2f s\n", step,
         ncell, wall);
  printf("cells failed (left unchanged for that step): %ld\n", nfailed);
  cudaFree(u_d); cudaFree(m_d); cudaFree(x_d); cudaFree(rho_d);
  cudaFree(T_d); cudaFree(dt_d); cudaFree(par_d); cudaFree(ref_d);
  cudaFree(flags_d);
  krome_opts_free(h);
  krome_device_finalize();
  return 0;
}
