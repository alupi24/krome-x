# Changelog

## Unreleased (after 2.0.0)

- **`n2x` takes the gas density** (interface change). Mass fractions are
  now x = m n / rho with the caller's rho, as in `x2n` and in the original
  KROME, since the network's species need not carry all the gas mass:
  C++ `n2x(n, rho, x)`, C `krome_n2x(ncell, n, rho, x)`, Fortran
  `krome_n2x(n, rho, x)`, Python `net.n2x(n, rho)`.
- **`DH` leaves out reactions whose energy is not the gas's.** Radiative
  channels (two reactants, one product), cosmic-ray and radiation-driven
  rates, and, with the `ATOMIC` cooling, the H/He collisional ionizations
  it already counts are no longer in the reaction-enthalpy terms; the
  generator lists them among its warnings, and a reaction's `; dh=0|1`
  overrides the choice.
- **In-bin photon spectrum** (changes results with coarse bins). Every bin
  average (photo cross sections, photoheating, `J_LW`, `h2pd()` lines,
  the Habing flux, grain Qabs) now assumes the same photon spectrum inside
  the bin, dN/dE ~ E^-alpha, set when the bins are built (C++
  `Radiation::alpha`, C `krome_opts_set_photo_alpha`, Fortran
  `krome_set_photo_alpha`, Python `opts.set_photo_alpha`). The default
  alpha = 1 (constant photons per ln E, J constant: legacy KROME) replaces
  2.0.0's flat photon spectrum (alpha = 0) of the photo() averages and
  the grain Qabs taken at the bin midpoint; alpha = 2 is the flat
  nu F_nu of He, Wibking & Krumholz (2024). The two Solomon forms
  (`J_LW` rates and photo()) now agree in bins straddling the band edges.
