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