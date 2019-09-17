program matrix_algebra
    implicit none
    real(8), dimension(5, 5) :: matrix
    real(8), dimension(5, 1) :: vector
    real(8), dimension(5, 1) :: vector_result
    real(8) :: scalar_result
    integer :: i
    integer :: ndim = 5

    matrix = reshape((/ 4, 2, 2, 5, 8, &
                        2, 5, 1, 3, 4, &
                        2, 1, 6, 2, 6, &
                        5, 3, 2, 1, 3, &
                        8, 4, 6, 3, 3 /), (/ 5, 5 /))
    vector = reshape((/ 2, 4, 5, 2, 1 /), (/ 5, 1 /))
    
    call calculate(matrix, vector, ndim, vector_result, scalar_result)

    print "(a14)", "The matrix A ="
    print "(5f8.2)", (matrix(i, :), i=1,5)
    print "(a21)", "The column vector V ="
    print "(f8.2)", (vector(i, 1), i=1,5)
    print "(a5)", "A*V ="
    print "(f8.2)", (vector_result(i, 1), i=1,5)
    print "(a9)", "V^T*A*V ="
    print "(f8.2)", scalar_result

end program

subroutine calculate(matrix, vector, ndim, vector_result, scalar_result)
!------------------------------------------------------------------------------------
! Given a matrix A and a vector V, this routine compute A*V and
! V^T*A*V, where * is the matrix multiplication.
!
! Arguments:
!           matrix: a square matrix with shape (ndim, ndim)
!           vector: a column vector with shape (ndim, 1)
!           ndim: the size of matrix and vector
!           vector_result: a column vector with shape (ndin, 1) to save the result A*V
!           scalar_result: a scalar to save the result V^T*A*V
!------------------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: ndim
    real(8), dimension(ndim, ndim), intent(in) :: matrix
    real(8), dimension(ndim, 1), intent(in) :: vector
    real(8), dimension(ndim, 1), intent(out) :: vector_result
    real(8), dimension(1, 1) :: fs
    real(8), intent(out) :: scalar_result

    vector_result = matmul(matrix, vector)
    fs = matmul(transpose(vector), vector_result)
    scalar_result = fs(1, 1)
    
    return
end subroutine