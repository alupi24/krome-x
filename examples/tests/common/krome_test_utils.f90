! helpers of the Fortran test driver (python3 -m kromex --test fortran), on the KROME-X
! Fortran binding (modules krome + krome_net); the legacy krome_user
! conveniences (krome_names_header, krome_print_best_flux,
! krome_explore_flux, krome_print_info, krome_scale_Z) are part of the
! binding. One cell: n(krome_nspec).
module krome_test_utils
  use krome
  use krome_net
  implicit none
  private
  public :: krome_free_fall_time, krome_check

contains

  ! free-fall time sqrt(3 pi / (32 G rho)) [s] of one cell
  function krome_free_fall_time(n) result(tff)
    real(8), intent(in) :: n(:)
    real(8) :: tff, nn(krome_nspec, 1), out(1)
    nn(:, 1) = n
    if (krome_get_free_fall_time(nn, out) /= 0) stop "free-fall time"
    tff = out(1)
  end function krome_free_fall_time

  ! stop with the solver's message on a failed solve
  subroutine krome_check(ierr, where)
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: where
    if (ierr == 0) return
    print *, "KROME failed (", trim(where), "): flag", ierr, " ", &
         trim(krome_flag_message_f(ierr))
    stop 1
  end subroutine krome_check

end module krome_test_utils
