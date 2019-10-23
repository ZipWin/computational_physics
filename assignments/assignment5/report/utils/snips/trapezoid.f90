subroutine trapezoid(f, a, b, n, I)
    implicit none
    real(8), external :: f
    real(8), intent(in) :: a, b
    integer, intent(in) :: n
    type(result), intent(out) :: I

    real(8) :: h, I_hh
    integer :: j

    h = (b - a) / n
    I%value = (f(a) + f(b)) * h/2.0d0
    do j = 1, n-1
        I%value = I%value + 2.0d0*f(a+j*h) * h/2.0d0
    end do

    h = (b - a) / (2*n)
    I_hh = (f(a) + f(b)) * h/2.0d0
    do j = 1, 2*n-1
        I_hh = I_hh + 2.0d0*f(a+j*h) * h/2.0d0
    end do

    I%error = 4.0d0/3.0d0 * (I_hh - I%value)

    return
end subroutine