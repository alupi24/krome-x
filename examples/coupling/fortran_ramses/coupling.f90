! KROME-X coupling example: a RAMSES-like Fortran grid code.
!
! The hydro state is RAMSES' conservative array uold(ncell, nvar): density,
! momentum, total energy, then the species as passive scalars rho*X_s (mass
! fractions X_s). Cells are processed in vectors of nvector cells, as
! RAMSES' cooling_fine does; each vector is packed into KROME's layout
! x(nspec, nvector) (one cell's species contiguous), integrated with ONE
! krome_solve_x call (every cell over the same hydro dt), and unpacked.
!
! Operator splitting per hydro step:
!   1. "hydro": a uniform compression rho -> 1.02 rho, the passive scalars
!      advected with it, adiabatic heating e -> e 1.02^(gamma-1)
!   2. chemistry + cooling: Tgas from the specific internal energy with
!      KROME's mu and gamma, krome_solve_x over dt, then the internal energy
!      from the new Tgas with the new mu and gamma (H2 formation changes
!      both)
! A cell whose solve fails (flag /= 0, the cell left untouched by KROME) is
! retried with two half steps, then four quarter steps (subcycling).
!
! Output: result.dat, per cell: index, initial nH [cm^-3], final Tgas [K],
! n(H2)/nH, n(e-)/nH (the same problem as the C, C++ and Python examples).
program coupling_ramses
  use krome
  use krome_net
  implicit none
  integer, parameter :: ncell = 16, nstep = 20, nvector = 8
  integer, parameter :: irho = 1, imom = 2, iene = 3, ivar0 = 3
  integer, parameter :: nvar = ivar0 + krome_nspec
  real(8), parameter :: kb = 1.380649d-16, mp = 1.67262192d-24
  real(8), parameter :: dt = 3d11, fcomp = 1.02d0
  real(8) :: uold(ncell, nvar), nH0(ncell), mass(krome_nspec)
  real(8) :: x(krome_nspec, nvector), n(krome_nspec, nvector)
  real(8) :: rho(nvector), Tgas(nvector), eint(nvector), dtv(nvector)
  real(8) :: mu(nvector), gam(nvector), par(krome_npar, 1)
  integer(c_int) :: flags(nvector)
  type(c_ptr) :: h
  integer :: i, c, k, s, ic, nv, istep
  real(8) :: ek

  do s = 1, krome_nspec
    mass(s) = krome_species_mass_f(s)
  end do

  ! solver options, once: thermal evolution on
  h = krome_opts_new()
  call krome_set(h, "evolve_T", 1d0)
  call krome_params_default(par)  ! one row, broadcast to every cell

  ! initial state: nH from 1 to 1e4 cm^-3, 1e4 K, 1% ionized
  do i = 1, ncell
    nH0(i) = 10d0**(4d0 * (i - 1) / (ncell - 1))
    n(:, 1) = 1d-20 * nH0(i)
    n(krome_idx_H, 1) = 0.99d0 * nH0(i)
    n(krome_idx_Hj, 1) = 0.01d0 * nH0(i)
    n(krome_idx_H2, 1) = 1d-6 * nH0(i)
    n(krome_idx_He, 1) = 0.0775d0 * nH0(i)
    n(krome_idx_E, 1) = 0.01d0 * nH0(i)
    uold(i, irho) = sum(mass * n(:, 1))
    uold(i, imom) = 0d0
    uold(i, ivar0 + 1:nvar) = mass * n(:, 1)  ! rho X_s
    Tgas(1) = 1d4
    call mu_gamma(n(:, 1:1), Tgas(1:1), mu(1:1), gam(1:1))
    uold(i, iene) = uold(i, irho) * kb * Tgas(1) / ((gam(1) - 1d0) &
         * mu(1) * mp)
  end do

  do istep = 1, nstep
    ! 1. hydro: uniform compression (momentum stays 0)
    do i = 1, ncell
      n(:, 1) = uold(i, ivar0 + 1:nvar) / mass
      Tgas(1) = temperature(uold(i, :), n(:, 1:1), eint(1))
      call mu_gamma(n(:, 1:1), Tgas(1:1), mu(1:1), gam(1:1))
      uold(i, irho) = uold(i, irho) * fcomp
      uold(i, ivar0 + 1:nvar) = uold(i, ivar0 + 1:nvar) * fcomp
      ek = 0.5d0 * uold(i, imom)**2 / uold(i, irho)
      uold(i, iene) = ek + uold(i, irho) * eint(1) * fcomp**(gam(1) - 1d0)
    end do

    ! 2. chemistry and cooling, nvector cells at a time
    do ic = 1, ncell, nvector
      nv = min(nvector, ncell - ic + 1)
      do k = 1, nv
        c = ic + k - 1
        rho(k) = uold(c, irho)
        x(:, k) = uold(c, ivar0 + 1:nvar) / rho(k)   ! mass fractions
        n(:, k) = x(:, k) * rho(k) / mass
        Tgas(k) = temperature(uold(c, :), n(:, k:k), eint(k))
      end do
      dtv(1:nv) = dt
      if (krome_solve_x(x(:, 1:nv), rho(1:nv), Tgas(1:nv), dtv(1:nv), h, &
           par=par, flags=flags(1:nv)) /= 0) then
        do k = 1, nv
          if (flags(k) /= 0) call subcycle(x(:, k), rho(k), Tgas(k), dt, &
               flags(k))
        end do
      end if
      do k = 1, nv
        c = ic + k - 1
        n(:, k) = x(:, k) * rho(k) / mass
        call mu_gamma(n(:, k:k), Tgas(k:k), mu(k:k), gam(k:k))
        uold(c, ivar0 + 1:nvar) = x(:, k) * rho(k)
        ek = 0.5d0 * uold(c, imom)**2 / uold(c, irho)
        uold(c, iene) = ek + rho(k) * kb * Tgas(k) / ((gam(k) - 1d0) &
             * mu(k) * mp)
      end do
    end do
  end do

  open(10, file="result.dat", status="replace")
  do i = 1, ncell
    n(:, 1) = uold(i, ivar0 + 1:nvar) / mass
    Tgas(1) = temperature(uold(i, :), n(:, 1:1), eint(1))
    if (krome_get_Hnuclei(n(:, 1:1), mu(1:1)) /= 0) stop "Hnuclei"
    write(10, '(I4,4ES20.11)') i - 1, nH0(i), Tgas(1), &
         n(krome_idx_H2, 1) / mu(1), n(krome_idx_E, 1) / mu(1)
  end do
  close(10)
  call krome_opts_free(h)
  print *, "wrote result.dat"

contains

  ! Tgas [K] of a cell from its specific internal energy e = (E - Ek)/rho
  function temperature(u, nc, e) result(T)
    real(8), intent(in) :: u(:), nc(:, :)
    real(8), intent(out) :: e
    real(8) :: T, m(1), g(1), T1(1)
    e = (u(iene) - 0.5d0 * u(imom)**2 / u(irho)) / u(irho)
    ! gamma may depend on Tgas (codegen thermal.gamma): two fixed-point passes
    T1 = 1d4
    call mu_gamma(nc, T1, m, g)
    T1 = e * (g(1) - 1d0) * m(1) * mp / kb
    call mu_gamma(nc, T1, m, g)
    T = e * (g(1) - 1d0) * m(1) * mp / kb
  end function temperature

  subroutine mu_gamma(nc, T, m, g)
    real(8), intent(in) :: nc(:, :), T(:)
    real(8), intent(out) :: m(:), g(:)
    if (krome_get_mu(nc, m) /= 0) stop "krome_get_mu"
    if (krome_get_gamma(nc, T, g, par=par) /= 0) stop "krome_get_gamma"
  end subroutine mu_gamma

  ! a failed cell: 2, then 4 substeps of dt (each from the last good state)
  subroutine subcycle(xc, rc, Tc, dtc, flag)
    real(8), intent(inout) :: xc(:), Tc
    real(8), intent(in) :: rc, dtc
    integer(c_int), intent(inout) :: flag
    real(8) :: xs(krome_nspec, 1), Ts(1), r(1), d(1)
    integer :: nsub, j, ierr
    do nsub = 2, 4, 2
      xs(:, 1) = xc
      Ts(1) = Tc
      r(1) = rc
      d(1) = dtc / nsub
      do j = 1, nsub
        ierr = krome_solve_x(xs, r, Ts, d, h, par=par)
        if (ierr /= 0) exit
      end do
      if (ierr == 0) then
        xc = xs(:, 1)
        Tc = Ts(1)
        flag = 0
        return
      end if
    end do
    print *, "cell failed after subcycling: ", trim(krome_flag_message_f(flag))
  end subroutine subcycle

end program coupling_ramses
