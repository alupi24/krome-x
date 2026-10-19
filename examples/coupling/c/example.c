/* KROME-X from C (the C interface, krome.h): the cell of the quick start
 * (README), primordial gas at nH = 1e4 cm^-3 and 100 K, followed in time
 * to 1e6 yr; the C++, Fortran and Python examples solve the same cell
 * (README.md) */
#include <stdio.h>
#include "krome.h"
#include "krome_net.h"

int main(void) {
  const double spy = 3.15e7, nH = 1e4;
  double n[KROME_NSPEC] = {0.0}, par[KROME_NPAR], Tgas = 100.0;
  double t = 0.0, dt = spy;
  int i, flag;
  FILE* f;

  /* solver options: Tgas evolved with the chemistry */
  krome_opts* h = krome_opts_new();
  krome_opts_set(h, "evolve_T", 1.0);

  /* parameters: the network defaults, and the free-fall time of the
   * COMPRESS heating */
  krome_params_default(par);
  par[KROME_PAR_tff] = 1e14;

  /* initial densities, cm^-3 */
  n[KROME_IDX_H] = nH;
  n[KROME_IDX_Hj] = 1.0;
  n[KROME_IDX_E] = 1.0;
  n[KROME_IDX_He] = 775.0;

  f = fopen("abundances.dat", "w");
  fprintf(f, "#t_yr Tgas");
  for (i = 0; i < KROME_NSPEC; ++i) fprintf(f, " %s", krome_species_name(i));
  fprintf(f, "\n");
  /* steps from 1 yr, 10% longer each, to 1e6 yr: one cell per call */
  while (t < 1e6 * spy) {
    flag = krome_solve_cells(1, n, &Tgas, &dt, par, 0, NULL, NULL, h);
    if (flag != 0) {
      printf("KROME failed: %s\n", krome_flag_message(flag));
      return 1;
    }
    t += dt;
    dt *= 1.1;
    fprintf(f, "%17.8E%17.8E", t / spy, Tgas);
    for (i = 0; i < KROME_NSPEC; ++i) fprintf(f, "%17.8E", n[i] / nH);
    fprintf(f, "\n");
  }
  fclose(f);
  krome_opts_free(h);
  printf("wrote abundances.dat\n");
  return 0;
}
