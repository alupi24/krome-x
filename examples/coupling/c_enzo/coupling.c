/* KROME-X coupling example: an ENZO-like C grid code.
 *
 * The grid keeps one baryon field per quantity, in code units: the gas
 * density, the specific internal energy, and one mass-density field per
 * species (ENZO's HIDensity, HIIDensity, ... fields). KROME wants the
 * number densities of a cell contiguous, n[ncell][nspec] [cm^-3], so the
 * fields are gathered once per step into that buffer (unit conversion
 * included), the WHOLE grid is integrated with ONE krome_solve_cells call
 * (every cell over the same hydro dt; the library parallelizes it with
 * OpenMP when built with -fopenmp), and scattered back.
 *
 * Operator splitting per hydro step:
 *   1. "hydro": a uniform compression rho -> 1.02 rho, the species fields
 *      advected with it, adiabatic heating e -> e 1.02^(gamma-1)
 *   2. chemistry + cooling: Tgas from e with KROME's mu and gamma,
 *      krome_solve_cells over dt, then e from the new Tgas with the new mu
 *      and gamma (H2 formation changes both)
 * A cell whose solve fails (flags[c] != 0, the cell left untouched by
 * KROME) is retried with two, then four substeps.
 *
 * Output: result.dat, per cell: index, initial nH [cm^-3], final Tgas [K],
 * n(H2)/nH, n(e-)/nH (the same problem as the Fortran, C++ and Python
 * examples). */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "krome.h"
#include "krome_net.h"

enum { NCELL = 16, NSTEP = 20 };
#define NS KROME_NSPEC
static const double KB = 1.380649e-16, MP = 1.67262192e-24;
static const double DT = 3e11, FCOMP = 1.02;

/* code units (ENZO's DensityUnits, VelocityUnits): cgs = code * units */
static const double DENSITY_UNITS = 1.67262192e-24;  /* g cm^-3 */
static const double VELOCITY_UNITS = 1e5;            /* cm s^-1 */

/* the grid: one field per quantity (ENZO BaryonField) */
typedef struct {
  double density[NCELL];      /* gas density [code] */
  double eint[NCELL];         /* specific internal energy [code] */
  double species[NS][NCELL];  /* species mass densities [code] */
} Grid;

static double mass[NS];  /* species masses [g] */
static double par[KROME_NPAR];
static krome_opts* h;

/* gather the species fields of the grid into n[NCELL][NS] [cm^-3] */
static void gather(const Grid* g, double* n) {
  for (int c = 0; c < NCELL; ++c)
    for (int s = 0; s < NS; ++s)
      n[c * NS + s] = g->species[s][c] * DENSITY_UNITS / mass[s];
}

static void scatter(const double* n, Grid* g) {
  for (int c = 0; c < NCELL; ++c)
    for (int s = 0; s < NS; ++s)
      g->species[s][c] = n[c * NS + s] * mass[s] / DENSITY_UNITS;
}

static void mu_gamma(int nc, const double* n, const double* T, double* mu,
                     double* gam) {
  if (krome_get_mu(nc, n, mu) || krome_get_gamma(nc, n, T, par, 0, gam)) {
    fprintf(stderr, "KROME query failed\n");
    exit(1);
  }
}

/* Tgas [K] of every cell from its specific internal energy; gamma may
 * depend on Tgas (codegen thermal.gamma): two fixed-point passes from 1e4 K */
static void temperature(const Grid* g, const double* n, double* T) {
  double mu[NCELL], gam[NCELL];
  for (int c = 0; c < NCELL; ++c) T[c] = 1e4;
  for (int pass = 0; pass < 2; ++pass) {
    mu_gamma(NCELL, n, T, mu, gam);
    for (int c = 0; c < NCELL; ++c)
      T[c] = g->eint[c] * VELOCITY_UNITS * VELOCITY_UNITS *
             (gam[c] - 1.0) * mu[c] * MP / KB;
  }
}

/* a failed cell: 2, then 4 substeps of dt (each from the last good state) */
static int subcycle(double* nc, double* Tc) {
  for (int nsub = 2; nsub <= 4; nsub += 2) {
    double ns[NS], Ts = *Tc, d = DT / nsub;
    int ierr = 0;
    for (int s = 0; s < NS; ++s) ns[s] = nc[s];
    for (int j = 0; j < nsub && !ierr; ++j)
      ierr = krome_solve_cells(1, ns, &Ts, &d, par, 0, NULL, NULL, h);
    if (!ierr) {
      for (int s = 0; s < NS; ++s) nc[s] = ns[s];
      *Tc = Ts;
      return 0;
    }
  }
  return -1;
}

int main(void) {
  static Grid g;
  double n[NCELL * NS], T[NCELL], dt[NCELL], mu[NCELL], gam[NCELL];
  double nH0[NCELL], nH[NCELL];
  int flags[NCELL];

  for (int s = 0; s < NS; ++s) mass[s] = krome_species_mass(s);
  h = krome_opts_new();                /* solver options, once */
  krome_opts_set(h, "evolve_T", 1.0);  /* thermal evolution on */
  krome_params_default(par);           /* one row, broadcast (stride 0) */

  /* initial state: nH from 1 to 1e4 cm^-3, 1e4 K, 1% ionized */
  for (int c = 0; c < NCELL; ++c) {
    double* nc = n + c * NS;
    nH0[c] = pow(10.0, 4.0 * c / (NCELL - 1));
    for (int s = 0; s < NS; ++s) nc[s] = 1e-20 * nH0[c];
    nc[KROME_IDX_H] = 0.99 * nH0[c];
    nc[KROME_IDX_Hj] = 0.01 * nH0[c];
    nc[KROME_IDX_H2] = 1e-6 * nH0[c];
    nc[KROME_IDX_He] = 0.0775 * nH0[c];
    nc[KROME_IDX_E] = 0.01 * nH0[c];
    T[c] = 1e4;
  }
  mu_gamma(NCELL, n, T, mu, gam);
  for (int c = 0; c < NCELL; ++c) {
    double rho = 0.0;
    for (int s = 0; s < NS; ++s) rho += mass[s] * n[c * NS + s];
    g.density[c] = rho / DENSITY_UNITS;
    g.eint[c] = KB * T[c] / ((gam[c] - 1.0) * mu[c] * MP) /
                (VELOCITY_UNITS * VELOCITY_UNITS);
  }
  scatter(n, &g);

  for (int step = 0; step < NSTEP; ++step) {
    /* 1. hydro: uniform compression */
    gather(&g, n);
    temperature(&g, n, T);
    mu_gamma(NCELL, n, T, mu, gam);
    for (int c = 0; c < NCELL; ++c) {
      g.density[c] *= FCOMP;
      for (int s = 0; s < NS; ++s) g.species[s][c] *= FCOMP;
      g.eint[c] *= pow(FCOMP, gam[c] - 1.0);
    }

    /* 2. chemistry and cooling: the whole grid in one call */
    gather(&g, n);
    temperature(&g, n, T);
    for (int c = 0; c < NCELL; ++c) dt[c] = DT;
    if (krome_solve_cells(NCELL, n, T, dt, par, 0, NULL, flags, h))
      for (int c = 0; c < NCELL; ++c)
        if (flags[c] && subcycle(n + c * NS, T + c))
          fprintf(stderr, "cell %d failed after subcycling: %s\n", c,
                  krome_flag_message(flags[c]));
    mu_gamma(NCELL, n, T, mu, gam);
    for (int c = 0; c < NCELL; ++c)
      g.eint[c] = KB * T[c] / ((gam[c] - 1.0) * mu[c] * MP) /
                  (VELOCITY_UNITS * VELOCITY_UNITS);
    scatter(n, &g);
  }

  gather(&g, n);
  temperature(&g, n, T);
  krome_get_Hnuclei(NCELL, n, nH);
  FILE* f = fopen("result.dat", "w");
  for (int c = 0; c < NCELL; ++c)
    fprintf(f, "%4d %19.11E %19.11E %19.11E %19.11E\n", c, nH0[c], T[c],
            n[c * NS + KROME_IDX_H2] / nH[c], n[c * NS + KROME_IDX_E] / nH[c]);
  fclose(f);
  krome_opts_free(h);
  printf("wrote result.dat\n");
  return 0;
}
