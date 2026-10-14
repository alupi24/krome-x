// This is a simple one-zone collapse test following the chemical and
// thermal evolution of a primordial cloud irradiated by both a UV flux
// below the Lyman limit (13.6 eV) and a X-ray flux between 2-10 keV.
// For further details we refer to Inayoshi, MNRAS, 416, 2748-2759, 2011.
// Note that the results are different compared to Inayoshi+2011 as the
// cross-sections employed here for H and He photoionization are from
// Verner+1996, while Inayoshi adopted older outdated cross-sections.
// The dynamics is described by the Larson-Penston-type similar solution.
// For a direct benchmark see Fig. 3 of Latif+2014.
// (KROME-X port of legacy tests/collapseUV_Xrays)
//
// KROME-X: the X-rays are photon bins. The legacy tabulated X-ray rates and
// heating (data/ratexH.dat, ...) were the H and He photoionization by
// J = J21_X 1e-21 (E / 1 keV)^-1.5 erg s^-1 cm^-2 Hz^-1 sr^-1 between 2
// and 10 keV, attenuated by the H and He Jeans columns (num2col), with the
// secondary ionization and heating of the fast photoelectrons. Here the
// driver attenuates the field by exp(-tau) of the same columns
// (get_opacity) before each step; the photoionization, its heating
// (thermal.heating = ["PHOTO"]) and the secondary processes (radiation.secondary_ionization)
// are computed by KROME.
#include <cmath>
#include <cstdio>
#include <vector>
#include "krome.hpp"
#include "krome_test_utils.hpp"
using namespace krome;

int main() {
  const int rstep = 500000;
  const double j21xs = 1e-1;  // J21 for X-rays

  // X-ray photon bins, 2-10 keV, and the unattenuated field: the J21
  // power law of the spectra is normalized at 13.6 eV, (13.6 / E)^alpha,
  // so J21_X at 1 keV is J21 = J21_X (1000 / 13.6)^1.5
  const int nbins = 100;
  Radiation rad = Radiation::log_bins(nbins, 2e3, 1e4);
  PhotoModule photo;
  photo.init(rad);
  std::vector<double> phi0(nbins), phi(nbins), tau(nbins);
  spectra::j21(nbins, rad.Elo.data(), rad.Ehi.data(),
               j21xs * std::pow(1e3 / 13.6, 1.5), 1.5, phi0.data());
  std::printf(" J21 Xray = %g\n", j21xs);

  // solver options: thermal evolution on; legacy -conserve (element and
  // charge conservation) is the runtime option conserve = 2
  SolverOpts o;
  o.evolve_T = true;
  o.conserve = 2;
  o.photo = &photo;
  o.Phi = phi.data();

  // preset J21 values for UV flux
  const double j21s[] = {0.0, 1e1, 1e3, 3e3};

  std::FILE* f22 = std::fopen("abundances.dat", "w");

  for (double j21 : j21s) {
    // INITIAL CONDITIONS: legacy krome_set_zredshift(15) set only the
    // redshift (Compton cooling); the CMB temperature (CIE cooling, lower
    // temperature cap) stayed at the legacy default 2.73 K
    Params p;
    p.zredshift = 15.0;
    p.Tcmb = 2.73;
    const double ntot = 1e-1;  // total density in 1/cm3
    double Tgas = 1.6e2;       // temperature in kelvin

    // KROME J21 parameter (legacy krome_set_user_J21)
    p.J21 = j21;

    // species initialization in 1/cm3
    double x[NSPEC];
    for (double& v : x) v = 1e-40;

    x[IDX_H] = 0.9999 * ntot;   // H
    x[IDX_H2] = 2e-6 * ntot;    // H2
    x[IDX_Hj] = 2e-3 * ntot;    // H+
    x[IDX_He] = 0.0775 * ntot;  // He
    x[IDX_E] = get_electrons_neutral(x);

    double dd = ntot;

    std::printf(" solving for J21= %g\n%5s%11s%11s\n", j21, "step",
                "n(cm-3)", "Tgas(K)");

    // output file header
    std::fprintf(f22, " #j21 ntot Tgas\n");

    // loop over the hydro time-step
    for (int i = 1; i <= rstep; ++i) {
      const double dd1 = dd;

      // free-fall time (legacy user_tff)
      const double tff = get_free_fall_time(x);
      p.tff = tff;
      const double dtH = 0.01 * tff;  // TIME-STEP
      dd += dd / tff * dtH;           // UPDATE DENSITY

      // rescale species number density to current density
      for (double& v : x) v *= dd / dd1;

      // run until 1e8 number density
      if (dd > 1e8) break;

      // the X-ray field attenuated by the H and He Jeans columns
      krome_test::check(get_opacity(x, Tgas, p, o, 0.0, nullptr,
                                    tau.data()), "get_opacity");
      for (int b = 0; b < nbins; ++b) phi[b] = phi0[b] * std::exp(-tau[b]);

      // solve the chemistry
      krome_test::check(solve(x, Tgas, p, dtH, o), "solve");

      std::fprintf(f22, "%17.8E%17.8E%17.8E\n", j21, dd, Tgas);
      if (i % 100 == 0) std::printf("%5d%11.3E%11.3E\n", i, dd, Tgas);
    }
    std::fprintf(f22, "\n");
  }
  std::fclose(f22);

  std::printf("To plot: python3 plot.py\nThat's all! have a nice day!\n");
  return 0;
}
