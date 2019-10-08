subroutine gauss_elimination(co_matrix, ndim, solution)
    implicit none
    integer(8), intent(in) :: ndim
    real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
    type(result), dimension(ndim), intent(out) :: solution

    real(8), dimension(ndim, ndim+1) :: A
    real(8), dimension(ndim) :: tmp_row, X
    integer(8) :: t, i

    A = co_matrix
    do t = 1, ndim
        ! pivot exchange
        do i = t, ndim
            if (abs(A(i, t)) > abs(A(t, t))) then
                tmp_row = A(t, :)
                A(t, :) = A(i, :)
                A(i, :) = tmp_row
            end if
        end do
        ! Gauss elimination
        do i = t+1, ndim
            A(i, :) = A(i, :) - A(i, t) / A(t, t) * A(t, :)
        end do
    end do

    print "('The Gauss-eliminated extented coefficients matrix is')"
    print "(10f7.2)", (A(i, :), i=1,9)

    X = 0.0d0
    do t = ndim, 1, -1
        X(t) = (A(t, ndim+1) - dot_product(A(t, t+1:ndim), X(t+1:ndim))) / A(t, t)
    end do

    solution%value = X
    solution%error = 0.0d0

    return
end subroutine