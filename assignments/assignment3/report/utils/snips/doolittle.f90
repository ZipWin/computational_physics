subroutine doolittle(co_matrix, ndim, solution)
    implicit none
    integer(8), intent(in) :: ndim
    real(8), dimension(ndim, ndim+1), intent(in) :: co_matrix
    type(result), dimension(ndim), intent(out) :: solution

    real(8), dimension(ndim, ndim+1) :: A
    real(8), dimension(ndim, ndim) :: L, U
    real(8), dimension(ndim) :: X, Y
    integer(8) :: i, j, k

    A = co_matrix
    L = 0.0d0; U = 0.0d0;
    
    call pivot_exchange(A, ndim, ndim+1)

    do k = 1, ndim
        do j = k, ndim
            U(k, j) = A(k, j) - dot_product(L(k, :k-1), U(:k-1, j))
        end do
        L(k, k) = 1.0d0
        do i = k+1, ndim
            L(i, k) = (A(i, k) - dot_product(L(i, :k-1), U(:k-1, k))) / U(k, k)
        end do
    end do

    do i = 1, ndim
        Y(i) = A(i, ndim+1) - dot_product(L(i, :i-1), Y(:i-1))
    end do
    do i = ndim, 1, -1
        X(i) = (Y(i) - dot_product(U(i, i+1:), X(i+1:))) / U(i, i)
    end do

    solution%value = X
    solution%error = 0.0d0

    return
end subroutine