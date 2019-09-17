function input_matrix(ndim1, ndim2) result(matrix)
    ! Input matrix manually
    implicit none
    integer, intent(in) :: ndim1, ndim2
    integer, parameter :: ndim = 100
    integer :: i
    real*8, dimension(ndim, ndim) :: matrix

    do i = 1, ndim1
        print "(a17,i2,a8)", "Please input the ", i, 'th row: '
        read *, matrix(i, :ndim2)
    end do

    return
end function  input_matrix