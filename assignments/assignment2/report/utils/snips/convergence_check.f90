subroutine convergence_check(phi, a, b, stat)
!------------------------------------------------------------------------------
! This subroutine check whether the iteration x = phi(x) converges in [a, b].
! It return the results in stat. If converges, stat is true; if not, false.
!
! Argments:
!          phi: a real(8) function
!          a: real(8), the left boundary of the interval
!          b: real(8), the right boundary of the interval
!          stat: logical, save the state of convergence
!------------------------------------------------------------------------------
    implicit none
    real(8), external :: phi
    real(8), intent(in) :: a, b
    logical, intent(out) :: stat

    real(8), dimension(EVAL_NUM) :: xs, phis, dphis
    real(8) :: eval_step
    integer :: i

    eval_step = abs(a - b) / EVAL_NUM

    do i = 1, EVAL_NUM
        xs(i) = a + eval_step * i
        phis(i) = phi(xs(i))
        dphis(i) = numerical_derivative(phi, xs(i))
    end do

    stat = .false.
    if (maxval(phis) <= maxval((/a, b/)) .and. minval(phis) >= minval((/a, b/))) then
        if (maxval(abs(dphis)) < 1) then
            stat = .true.
        end if
    end if

    return
end subroutine