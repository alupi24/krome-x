// In this test we have a semi-infinite slab in a box with size rmin-rmax
// defined by ngrid grid points. The light is attenuated via the opacity
// determined by the photochemical xsecs present in the network and the
// dust opacity loaded from file. The emission is 10 times the Draine flux
// in the energy interval [5,13.6] eV. (KROME-X port of legacy tests/rtchem)
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int ngrid = 300;  // grid points
  const double fav = 1.0 / 6.289e-22;  // Av/NH conversion factor
  const double ntot = 1e3;             // total density
  const double spy = 365.0 * 24.0 * 3600.0;  // s
  const double rmin = fav / ntot * 1e-4;  // min radius: cm
  const double rmax = fav / ntot * 3.0;   // box radius: cm
  const double tend = 1e6 * spy;
  const double d2g = 1e-2;  // dust-to-gas mass ratio of the dust opacity

  // gas temperature in K (all cells), fixed
  std::vector<double> Tgas(ngrid, 5e1);

  // runtime parameters: the cosmic-ray rate (legacy user_crate)
  Params p;
  p.crate = 5e-17;

  // init abundances
  std::vector<double> xall(static_cast<size_t>(ngrid) * NSPEC, 0.0);
  for (int i = 0; i < ngrid; ++i) {
    double* x = &xall[static_cast<size_t>(i) * NSPEC];
    x[IDX_H] = ntot / 3.0;
    x[IDX_H2] = ntot / 3.0 * 2.0;
    x[IDX_C] = 1e-4 * ntot;
    x[IDX_O] = 3e-4 * ntot;
    x[IDX_He] = 1e-1 * ntot;
  }

  // initial time step
  double dt = .1 * spy;
  double t = 0.0;

  // the position of the grid points
  std::vector<double> rx(ngrid);
  for (int i = 0; i < ngrid; ++i)
    rx[i] = std::pow(1e1, i * (std::log10(rmax) - std::log10(rmin))
                              / (ngrid - 1) + std::log10(rmin));

  // photon bins (legacy -photoBins=10): 10 linear bins in [5, 13.6] eV
  // with the Draine flux scaled by 10 (the stored, unattenuated field)
  const int B = 10;
  Radiation rad;
  rad.Elo.resize(B);
  rad.Ehi.resize(B);
  rad.Phi.resize(B);
  spectra::bins_lin(B, 5.0, 13.6, rad.Elo.data(), rad.Ehi.data());
  spectra::draine(B, rad.Elo.data(), rad.Ehi.data(), rad.Phi.data());
  for (double& v : rad.Phi) v *= 1e1;
  PhotoModule photo;
  photo.init(rad);

  // energy (bin midpoint, for the plot) and the stored intensity J
  // [eV cm^-2 s^-1 Hz^-1 sr^-1]
  std::vector<double> energy(B), flux(B);
  for (int j = 0; j < B; ++j) {
    energy[j] = 0.5 * (rad.Elo[j] + rad.Ehi[j]);
    flux[j] = spectra::phi_to_J(rad.Phi[j], rad.Elo[j], rad.Ehi[j]);
  }

  // load energy (micron), opacity (cm2/g) table from file, and average it
  // on the bins; kdust: the dust opacity per gas mass, d2g kappa
  std::vector<double> tE, tk;
  {
    std::FILE* ft = std::fopen("opacityDraineR31.dat", "r");
    if (ft == nullptr) {
      std::printf("ERROR: missing opacityDraineR31.dat\n");
      return 1;
    }
    double wl, k;
    while (std::fscanf(ft, "%lf %lf", &wl, &k) == 2) {
      tE.push_back(1.23984193 / wl);  // micron -> eV
      tk.push_back(k);
    }
    std::fclose(ft);
  }
  if (tE.size() > 1 && tE[1] < tE[0]) {  // ascending energies
    std::vector<double> e2(tE.rbegin(), tE.rend()), k2(tk.rbegin(), tk.rend());
    tE.swap(e2);
    tk.swap(k2);
  }
  std::vector<double> kdust(B);
  if (spectra::dust_kappa_bins(static_cast<int>(tE.size()), tE.data(),
                               tk.data(), B, rad.Elo.data(), rad.Ehi.data(),
                               kdust.data()) != 0) {
    std::printf("ERROR: the opacity table does not cover the bins\n");
    return 1;
  }
  for (double& v : kdust) v *= d2g;

  // solver options: defaults (isothermal), the photon bins and the cell's
  // field
  std::vector<double> phi(B), op(B), kabs(B);
  SolverOpts o;
  o.photo = &photo;
  o.Phi = phi.data();

  // the composition that attenuates cell i: as in legacy, the state left
  // by the previous solve (cell i-1; for the first cell the last cell of
  // the previous time step, zero at the first step)
  std::vector<double> x(NSPEC, 0.0);

  std::FILE* f66 = std::fopen("profiles.dat", "w");
  std::FILE* f77 = std::fopen("abundances_center.dat", "w");
  std::FILE* f55 = std::fopen("profiles_final.dat", "w");
  const std::string hdr = names_header();
  std::fprintf(f66, "#time Av Tgas %s\n", hdr.c_str());
  std::fprintf(f77, "#time Tgas %s\n", hdr.c_str());
  std::fprintf(f55, "#time Av Tgas %s\n", hdr.c_str());

  // loop on time
  for (;;) {
    // a file with the opacity (replaced every time step)
    std::FILE* f45 = std::fopen("opacity.dat", "w");
    dt *= 1.2;  // increase time-step size
    t += dt;    // increase total time
    std::printf("%12.4f %%\n", t / tend * 1e2);
    double rold = 0.0;
    for (double& v : op) v = 1.0;  // initial opacity is 1 (thin)
    double Av = 0.0;
    // loop on grid points
    for (int i = 0; i < ngrid; ++i) {
      double* xi = &xall[static_cast<size_t>(i) * NSPEC];
      const double r = rx[i];
      const double dr = r - rold;  // grid spacing
      // compute Av
      double sx = 0.0;
      for (int k = 0; k < NSPEC; ++k) sx += xi[k];
      Av += dr * sx / fav;
      // attenuation by photochemistry and dust across the cell
      krome_test::check(get_absorption(x.data(), Tgas[i], p, o, kdust.data(),
                                       kabs.data()), "get_absorption");
      for (int j = 0; j < B; ++j) op[j] *= std::exp(-kabs[j] * dr);
      // energy-dependent flux
      for (int j = 0; j < B; ++j)
        std::fprintf(f45, "%17.8E %17.8E %17.8E\n", Av, energy[j],
                     op[j] * flux[j]);
      std::fprintf(f45, "\n");
      // the flux of the cell: the stored flux scaled by the opacity, over
      // 4 (a semi-infinite slab: krome assumes 4 pi)
      for (int j = 0; j < B; ++j) phi[j] = rad.Phi[j] * op[j] / 4.0;
      // call KROME to do chemistry
      krome_test::check(solve(xi, Tgas[i], p, dt, o), "solve");
      for (int k = 0; k < NSPEC; ++k) x[k] = xi[k];
      rold = r;
    }
    std::fclose(f45);

    // the time evolution of the ngrid/2 cell
    krome_test::write_row(f77, {t / spy, Tgas[ngrid / 2 - 1]},
                          &xall[static_cast<size_t>(ngrid / 2 - 1) * NSPEC]);

    // full evolution (all cells, all time-steps)
    rold = 0.0;
    Av = 0.0;
    for (int i = 0; i < ngrid; ++i) {
      const double* xi = &xall[static_cast<size_t>(i) * NSPEC];
      double sx = 0.0;
      for (int k = 0; k < NSPEC; ++k) sx += xi[k];
      Av += (rx[i] - rold) * sx / fav;
      krome_test::write_row(f66, {t / spy, Av, Tgas[i]}, xi);
      rold = rx[i];
    }
    std::fprintf(f66, "\n");

    if (t > tend) break;  // exit after tend
  }

  // final conditions
  double Av = 0.0, rold = 0.0;
  for (int i = 0; i < ngrid; ++i) {
    const double* xi = &xall[static_cast<size_t>(i) * NSPEC];
    double sx = 0.0;
    for (int k = 0; k < NSPEC; ++k) sx += xi[k];
    Av += (rx[i] - rold) * sx / fav;
    krome_test::write_row(f55, {t / spy, Av, Tgas[i]}, xi);
    rold = rx[i];
  }
  std::fclose(f66);
  std::fclose(f77);
  std::fclose(f55);

  std::printf("done!\nTo plot: python3 plot.py\n");
  return 0;
}
