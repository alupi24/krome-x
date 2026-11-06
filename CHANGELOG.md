# Changelog

## Unreleased (after 2.0.0)

- **Background photon field.** With `radiation.evolve`, a fixed field per
  bin (`Phi_bg`, e.g. a UV background) can be added to the evolved one:
  it acts wherever the field acts (rates, heating, grains) but is never
  absorbed nor returned. One array for every cell (C++
  `SolverOpts::Phi_bg`, C `krome_opts_set_photo_background`, Fortran
  `krome_set_photo_background`, Python `opts.set_photo_background`) or
  one row per cell (`CellsOpts::Phi_bg`, C `krome_cells.Phi_bg`, Fortran
  `krome_solve(..., phi_bg=)`, Python `net.solve(..., Phi_bg=)`). In C,
  `krome_cells` has a fourth member, `Phi_bg`, at the end: initialize
  the struct (`{0}` or `{Phi, Td, adust}`) so that it is NULL when unused.
- **`n2x` takes the gas density** (interface change). Mass fractions are
  now x = m n / rho with the caller's rho, as in `x2n` and in the original
  KROME, since the network's species need not carry all the gas mass:
  C++ `n2x(n, rho, x)`, C `krome_n2x(ncell, n, rho, x)`, Fortran
  `krome_n2x(n, rho, x)`, Python `net.n2x(n, rho)`.
