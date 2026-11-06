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
