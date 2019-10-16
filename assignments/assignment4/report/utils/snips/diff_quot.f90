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