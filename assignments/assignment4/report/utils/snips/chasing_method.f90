subroutine chasing_method(ndim, s, a, u, d, x)
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim), intent(in) :: a, d
    real(8), dimension(ndim-1), intent(in) :: s, u
    real(8), dimension(ndim), intent(out) :: x

    real(8), dimension(ndim) :: alpha, y
    real(8), dimension(ndim-1) :: g, b
    integer :: i

    do i = 1, ndim-1
        b(i) = u(i)
    end do
    alpha(1) = a(1)
    do i = 2, ndim
        g(i-1) = s(i-1)/alpha(i-1)
        alpha(i) = a(i) - g(i-1)*b(i-1)
    end do

    y(1) = d(1)
    do i = 2, ndim
        y(i) = d(i) - g(i-1)*y(i-1)
    end do
    x(ndim) = y(ndim)/alpha(ndim)
    do i = ndim-1, 1, -1
        x(i) = (y(i)-b(i)*x(i+1))/alpha(i)
    end do

    return
end subroutine