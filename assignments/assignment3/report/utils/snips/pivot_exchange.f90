subroutine pivot_exchange(matrix, ndim1, ndim2)
    implicit none
    integer(8), intent(in) :: ndim1, ndim2
    real(8), dimension(ndim1, ndim2) :: matrix

    integer :: col, i
    real(8), dimension(ndim2) :: tmp_row

    do col = 1, ndim2
        do i = col, ndim1
            if (matrix(i, col) > matrix(col, col)) then
                tmp_row = matrix(col, :)
                matrix(col, :) = matrix(i, :)
                matrix(i, :) = tmp_row
            end if
        end do
    end do

    return
end subroutine