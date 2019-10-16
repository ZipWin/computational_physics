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