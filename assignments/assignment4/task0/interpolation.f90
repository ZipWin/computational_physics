module utils
    implicit none
    type result
        real(8) :: value, error
    end type

    contains
    recursive real(8) function diff_quot(points) result(dq)
        implicit none
        real(8), dimension(*, *), intent(in) :: points(:, :)

        integer :: ndim

        ndim = size(points, dim=1)

        if (ndim == 1) then
            dq = points(1, 2)
        else if (ndim >= 2) then
            dq = (diff_quot(points(2:ndim, :)) - diff_quot(points(1:ndim-1, :))) / (points(ndim, 1) - points(1, 1))
        else
            stop "ERROR: the size of list for differential quotient must greater then 2"
        end if

        return
    end function

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

end module


module interpolation
    use utils
    contains
    recursive type(result) function lagrange(sample, x, require_error) result(y)
        implicit none
        real(8), dimension(*, *), intent(in) :: sample(:, :)
        real(8), intent(in) :: x
        logical, optional :: require_error

        real(8) :: item
        type(result) :: L, R
        integer :: i, j, nos

        if (.not. present(require_error)) then
            require_error = .true.
        end if

        nos = size(sample, dim=1)

        y%value = 0.0d0
        do i = 1, nos
            item = 1.0d0
            do j = 1, i-1
                item = item * (x - sample(j, 1)) / (sample(i, 1) - sample(j, 1))
            end do
            do j = i+1, nos
                item = item * (x - sample(j, 1)) / (sample(i, 1) - sample(j, 1))
            end do
            y%value = y%value + item * sample(i, 2)
        end do
        if (require_error) then
            L = lagrange(sample(1:nos-2, :), x, require_error=.false.)
            R = lagrange(sample(3:nos  , :), x, require_error=.false.)
            y%error = (x-sample(nos, 1))*(x-sample(nos-1, 1))/&
            ((x-sample(nos, 1))*(x-sample(nos-1, 1))-(x-sample(1, 1))*(x-sample(2, 1)))*L%value&
             - (x-sample(1, 1))*(x-sample(2, 1))/&
            ((x-sample(nos, 1))*(x-sample(nos-1, 1))-(x-sample(1, 1))*(x-sample(2, 1)))*R%value&
            -(L%value + R%value)/2.0d0
        end if

        return
    end function

    recursive type(result) function newton(sample, x, require_error) result(y)
        implicit none
        real(8), dimension(*, *), intent(in) :: sample(:, :)
        real(8), intent(in) :: x
        logical, optional :: require_error

        real(8) :: item
        type(result) :: L, R
        integer :: i, j, nos

        if (.not. present(require_error)) then
            require_error = .true.
        end if

        nos = size(sample, dim=1)

        y%value = 0.0d0
        do i = 1, nos
            item = 1.0d0
            do j = 1, i-1
                item = item * (x - sample(j, 1))
            end do
            y%value = y%value + item * diff_quot(sample(1:i, :))
        end do
        if (require_error) then
            L = newton(sample(1:nos-2, :), x, require_error=.false.)
            R = newton(sample(3:nos  , :), x, require_error=.false.)
            y%error = (x-sample(nos, 1))*(x-sample(nos-1, 1))/&
            ((x-sample(nos, 1))*(x-sample(nos-1, 1))-(x-sample(1, 1))*(x-sample(2, 1)))*L%value&
             - (x-sample(1, 1))*(x-sample(2, 1))/&
            ((x-sample(nos, 1))*(x-sample(nos-1, 1))-(x-sample(1, 1))*(x-sample(2, 1)))*R%value&
            -(L%value + R%value)/2.0d0
        end if

        return
    end function

    type(result) function cubic_spline(sample, x) result(y)
        implicit none
        real(8), dimension(*, *), intent(in) :: sample(:, :)
        real(8), intent(in) :: x

        integer :: nos
        real(8), dimension(size(sample, dim=1)) :: m
        real(8), dimension(size(sample, dim=1)-1) :: h
        real(8), dimension(size(sample, dim=1)-2) :: a, b, q, dia
        real(8), dimension(size(sample, dim=1)-3) :: u
        real(8) :: L, R
        integer :: i

        nos = size(sample, dim=1)

        do i = 1, nos-1
            h(i) = sample(i+1, 1) - sample(i, 1)
        end do
        do i = 1, nos-2
            a(i) = h(i) / (h(i) + h(i+1))
            b(i) = 3.0d0 * ((1.0d0-a(i))*(sample(i+1, 2)-sample(i, 2))/h(i) + a(i)*(sample(i+2, 2)-sample(i+1, 2))/h(i+1))
        end do
        
        m(1) = 0.0d0
        m(nos) = 0.0d0
        ! solve m by chasing method
        dia(:) = 2.0d0
        call chasing_method(nos-2, 1.0d0-a(2:nos-2), dia, a(1:nos-3), b(:), m(2:nos-1))

        do i = 1, nos-1
            if (x >= sample(i, 1) .and. x <= sample(i+1, 1)) then
                exit
            end if
        end do

        y%value = sample(i, 2) * (h(i)+2*(x-sample(i, 1))) * (x-sample(i+1, 1))**2 / h(i)**3
        y%value = y%value + sample(i+1, 2) * (h(i)-2*(x-sample(i+1, 1))) * (x-sample(i, 1))**2 / h(i)**3
        y%value = y%value + m(i) * (x-sample(i, 1)) * (x-sample(i+1, 1))**2 / h(i)**2
        y%value = y%value + m(i+1) * (x-sample(i+1, 1)) * (x-sample(i, 1))**2 / h(i)**2

        y%error = 1.0d0/(1.0d0+x**2) - y%value

        return
    end function

end module


program main
    use utils
    use interpolation
    implicit none
    integer, parameter :: N = 15, NS = 9999
    real(8), dimension(N+1, 2) :: sample
    real(8) :: x
    type(result) :: y
    integer :: i

    open(file='samples.dat', unit=10, action='write')
    do i = 0, N
        sample(i+1, 1) = -5.0d0 + 10.0d0/N * i
        sample(i+1, 2) = 1.0d0 / (1.0d0 + sample(i+1, 1)**2)
        write (10, "(2f16.8)") sample(i+1, 1), sample(i+1, 2)
    end do
    close(unit=10)
    
    open(file='lagrange.dat', unit=10, action='write')
    do i = 0, NS
        x = -5.0d0 + 10.0d0/NS * i
        y = lagrange(sample, x, require_error=.true.)
        write (10, "(4f16.8)") x, 1.0d0/(1.0d0+x**2), y%value, y%error
    end do
    close(unit=10)

    open(file='newton.dat', unit=10, action='write')
    do i = 0, NS
        x = -5.0d0 + 10.0d0/NS * i
        y = newton(sample, x, require_error=.true.)
        write (10, "(4f16.8)") x, 1.0d0/(1.0d0+x**2), y%value, y%error
    end do
    close(unit=10)

    open(file='cubic_spline.dat', unit=10, action='write')
    do i = 0, NS
        x = -5.0d0 + 10.0d0/NS * i
        y = cubic_spline(sample, x)
        write (10, "(4f16.8)") x, 1.0d0/(1.0d0+x**2), y%value, y%error
    end do
    close(unit=10)

end program